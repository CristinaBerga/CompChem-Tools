# -*- coding: utf-8 -*-

####################################################################
# Script: NMR.py
# Author: Cristina Berga, https://github.com/CristinaBerga/CompChem-Tools/
# Creation Date: June 2025
####################################################################
#
# This Python script is a tool designed for the automated analysis
# of Nuclear Magnetic Resonance (NMR) calculations from Gaussian
# output files.
#
####################################################################
#
# Key Features:
#
# 1. Extraction of key data, including the total energy (SCF Done)
#    and isotropic magnetic shielding values (GIAO) for ¹H and ¹³C
#    nuclei.
#
# 2. Automatic calculation of chemical shifts using a TMS
#    (Tetramethylsilane) reference file.
#
# 3. Generation of a CSV file (NMR.csv) with all extracted
#    and calculated values, organized by molecule.
#
# 4. Creation of SVG (Scalable Vector Graphics) images of the
#    molecules, with the corresponding chemical shifts embedded
#    as labels for quick visualization.
#
####################################################################
#
# Usage:
#
# - Run the script and follow the prompts to select the nuclei
#   (C, H, or both).
#
# - You can specify the .out files directly on the command line
#   or provide the path to a folder containing them.
#
# - Ensure that the TMS reference file (`TMS.out`) is located in
#   the same directory as the other files.
#
# - The installation of Open Babel is required for format conversion.
#
####################################################################

import os
import sys
import re
import csv
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D
import matplotlib.pyplot as plt
import numpy as np

HARTREE_TO_KCAL_MOL = 627.509

def demanar_elements():
    resposta = ""
    while resposta not in ["C", "H", "CH", "HC"]:
        resposta = input("Which NMR do you want? C, H or both (CH)?").strip().upper()
    return set(resposta)

def extreure_dades_fitxer(nom_fitxer, elements_interessants):
    energia_scf = None
    shieldings = {}

    with open(nom_fitxer, 'r') as f:
        lines = f.readlines()

    normal_termination_indices = []
    for i, line in enumerate(lines):
        if "Normal termination of Gaussian" in line:
            normal_termination_indices.append(i)

    search_range_for_scf = lines
    
    if normal_termination_indices:
        first_normal_termination_idx = normal_termination_indices[0]
        search_range_for_scf = lines[:first_normal_termination_idx]
        
    for line in reversed(search_range_for_scf):
        if "SCF Done" in line:
            if "=" in line:
                energia_str = line.split("=")[1].strip().split()[0]
                try:
                    energia_scf = float(energia_str)
                except ValueError:
                    energia_scf = None
            break
    
    in_shielding_block = False
    for line in lines:
        if "SCF GIAO Magnetic shielding tensor" in line:
            in_shielding_block = True
            shieldings = {}
            continue
        if in_shielding_block:
            match = re.match(r"\s*(\d+)\s+(\w+)\s+Isotropic\s*=\s*([-+]?\d*\.\d+)", line)
            if match:
                atom_index = match.group(1)
                atom_type = match.group(2)
                shielding = float(match.group(3))
                if atom_type in elements_interessants:
                    atom_label = atom_index + atom_type
                    shieldings[atom_label] = shielding
            elif line.strip() == "" and "Isotropic" not in line:
                break

    return energia_scf, shieldings

def obtenir_referencia_tms(tms_fitxer):
    with open(tms_fitxer, 'r') as f:
        lines = f.readlines()

    tms_C = None
    tms_H = None
    for line in lines:
        if "SCF GIAO Magnetic shielding tensor" in line:
            break

    for line in lines[lines.index(line)+1:]:
        match = re.match(r"\s*(\d+)\s+(\w+)\s+Isotropic\s*=\s*([-+]?\d*\.\d+)", line)
        if match:
            atom_type = match.group(2)
            shielding = float(match.group(3))
            if atom_type == "C" and tms_C is None:
                tms_C = shielding
            elif atom_type == "H" and tms_H is None:
                tms_H = shielding
            if tms_C and tms_H:
                break
    return {"C": tms_C, "H": tms_H}

def extreure_xyz(nom_fitxer, xyz_fitxer):
    with open(nom_fitxer, 'r') as f:
        lines = f.readlines()

    geometria = []
    bloc_inici = None
    for i in range(len(lines) - 1, -1, -1):
        if "Standard orientation" in lines[i] or "Input orientation" in lines[i]:
            bloc_inici = i
            break
    if bloc_inici is None:
        return False

    bloc_dades_inici = bloc_inici + 5
    simbols = Chem.GetPeriodicTable()
    with open(xyz_fitxer, 'w') as out_xyz:
        count = 0
        geometria_xyz = []
        for i in range(bloc_dades_inici, len(lines)):
            if "-----" in lines[i] or lines[i].strip() == "":
                break
            parts = lines[i].split()
            if len(parts) >= 6:
                try:
                    atomic_number = int(parts[1])
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    geometria_xyz.append((simbols.GetElementSymbol(atomic_number), x, y, z))
                    count += 1
                except ValueError:
                    continue
        out_xyz.write(f"{count}\n\n")
        for sym, x, y, z in geometria_xyz:
            out_xyz.write(f"{sym:2} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n")
    return True

def convertir_xyz_a_mol(nom_base_path):
    xyz_file = f"{nom_base_path}.xyz"
    mol_file = f"{nom_base_path}.mol"
    if os.path.exists(xyz_file):
        subprocess.run(["obabel", "-ixyz", xyz_file, "-O", mol_file, "-f", "-k"], check=True)
        return mol_file
    return None

def genera_imatge(mol, desplaçaments_corr, fitxer_nom_path):
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)

    conf3D = None
    if mol.GetNumConformers() > 0:
        conf3D = mol.GetConformer(0)

    if conf3D:
        Chem.AssignStereochemistryFrom3D(mol)
    
    AllChem.Compute2DCoords(mol)
    conf2D = mol.GetConformer()

    # Adjust the position of hydrogens for label display
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            idx_H = atom.GetIdx()
            veïns = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
            if len(veïns) == 1:
                idx_X = veïns[0]
                pos_H = np.array(conf2D.GetAtomPosition(idx_H))
                pos_X = np.array(conf2D.GetAtomPosition(idx_X))
                direcció = pos_H - pos_X
                nova_pos = pos_X + 0.6 * direcció # Move the label away from the H
                conf2D.SetAtomPosition(idx_H, Point3D(*nova_pos))

    width, height = 1000, 1000
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()

    # Add the labels for the displacements
    atom_labels = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        key = f"{idx+1}{atom.GetSymbol()}"
        if key in desplaçaments_corr:
            ppm = desplaçaments_corr[key]
            atom_labels[idx] = f"{ppm:.2f}"

    for idx, label in atom_labels.items():
        opts.atomLabels[idx] = label

    # Ensure that stereo annotations are enabled
    opts.addStereoAnnotation = True
    opts.addAtomIndices = False

    # Draw the molecule
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()

    svg_text = drawer.GetDrawingText()
    with open(f"{fitxer_nom_path}.svg", "w") as f:
        f.write(svg_text)

def processar_fitxers(llista_fitxers, elements_interessants, output_dir):
    fitxers = sorted(llista_fitxers)
    tots_labels_potencials = set()
    dades = []

    tms_fitxer = next((f for f in fitxers if os.path.basename(f).startswith("TMS")), None)
    tms_refs = obtenir_referencia_tms(tms_fitxer) if tms_fitxer else {}

    for fitxer in fitxers:
        base_nom = os.path.splitext(os.path.basename(fitxer))[0]
        xyz_path = os.path.join(output_dir, f"{base_nom}.xyz")
        mol_path = os.path.join(output_dir, f"{base_nom}.mol")
        svg_path_base = os.path.join(output_dir, base_nom)

        extret = extreure_xyz(fitxer, xyz_path)
        if extret:
            convertir_xyz_a_mol(os.path.join(output_dir, base_nom))

        scf, shieldings = extreure_dades_fitxer(fitxer, elements_interessants)
        tots_labels_potencials.update(shieldings.keys())
        dades.append((os.path.basename(fitxer), scf, shieldings))

        desplaçaments_corr = {}
        for label, val in shieldings.items():
            tipus_atom = re.sub(r'\d+', '', label)
            tms_val = tms_refs.get(tipus_atom)
            if tms_val:
                desplaçaments_corr[label] = tms_val - val

        if os.path.exists(mol_path):
            mol = Chem.MolFromMolFile(mol_path, removeHs=False, sanitize=False)
                
            if mol:
                try:
                    genera_imatge(mol, desplaçaments_corr, svg_path_base)
                except Exception as e:
                    print(f"ERROR: Unable to generate SVG for {base_nom} with RDKit: {e}")
            else:
                print(f"WARNING: RDKit could not load the .mol file ({mol_path}). Skipping image generation.")
        else:
            print(f"WARNING: The .mol file ({mol_path}) does not exist. Skipping image generation.")

    dades_per_csv = []
    
    labels_amb_valors = set()
    for fitxer_original, scf, shieldings in dades:
        if os.path.basename(fitxer_original).startswith("TMS"):
            continue

        energia_kcal_mol = None
        if scf is not None:
            energia_kcal_mol = round(scf * HARTREE_TO_KCAL_MOL, 3)

        fila_corregida = {
            "Molecule": os.path.splitext(os.path.basename(fitxer_original))[0],
            "Energy (hartree)": scf,
            "Energy (kcal/mol)": energia_kcal_mol
        }
        
        desplaçaments_corr = {}
        for label, val in shieldings.items():
            tipus_atom = re.sub(r'\d+', '', label)
            tms_val = tms_refs.get(tipus_atom)
            if tms_val:
                try:
                    desplaçaments_corr[label] = round(tms_val - float(val), 3)
                except (TypeError, ValueError):
                    desplaçaments_corr[label] = ""
            else:
                desplaçaments_corr[label] = "" 

        for label in tots_labels_potencials:
            val_corregit = desplaçaments_corr.get(label, "")
            fila_corregida[label] = val_corregit
            if val_corregit != "":
                labels_amb_valors.add(label)
        
        dades_per_csv.append(fila_corregida)

    labels_finals_ordenats = sorted(list(labels_amb_valors), key=lambda x: int(re.findall(r'\d+', x)[0]))

    # Generate the CSV in the output folder
    csv_output_path = os.path.join(output_dir, "NMR.csv")
    with open(csv_output_path, "w", newline='') as f:
        writer = csv.writer(f, delimiter=';')
        
        header = ["Molecule", "Energy (hartree)", "Energy (kcal/mol)"] + labels_finals_ordenats
        writer.writerow(header)
        
        for fila_data in dades_per_csv:
            fila_csv = [
                fila_data["Molecule"],
                fila_data["Energy (hartree)"],
                fila_data["Energy (kcal/mol)"]
            ]
            for label in labels_finals_ordenats:
                fila_csv.append(fila_data.get(label, ""))
            writer.writerow(fila_csv)

    # Cleanup: remove .xyz and .mol files from the output folder
    for f in os.listdir(output_dir):
        if f.endswith(".xyz") or f.endswith(".mol"):
            try:
                os.remove(os.path.join(output_dir, f))
            except Exception as e:
                print(f"Unable to delete {f} from {output_dir}: {e}")

    tms_fitxers_out_in_output_dir = sorted([f for f in os.listdir(output_dir) if f.startswith("TMS") and f.endswith(".out")])
    for f in tms_fitxers_out_in_output_dir:
        full_path_f = os.path.join(output_dir, f)
        if tms_fitxer and os.path.abspath(full_path_f) == os.path.abspath(tms_fitxer):
            continue
        try:
            os.remove(full_path_f)
        except Exception as e:
            print(f"Unable to delete {f} from {output_dir}: {e}")


if __name__ == "__main__":
    elements = demanar_elements()
    
    fitxers_input = []
    output_directory = ""

    if len(sys.argv) > 1 and all(f.endswith(".out") for f in sys.argv[1:]):
        fitxers_input = [os.path.abspath(f) for f in sys.argv[1:]]
        if fitxers_input:
            output_directory = os.path.dirname(fitxers_input[0])
        else:
            output_directory = os.getcwd()
    else:
        carpeta_input = input("Introduce the path to the folder containing the .out files (leave blank for current folder): ").strip()
        carpeta_input = carpeta_input if carpeta_input else "."
        
        if not os.path.isdir(carpeta_input):
            print(f"Error: The folder '{carpeta_input}' does not exist or is not a valid directory.")
            sys.exit(1)
        
        fitxers_input = [os.path.join(carpeta_input, f) for f in os.listdir(carpeta_input) if f.endswith(".out")]
        output_directory = carpeta_input

    if not fitxers_input:
        print("No .out files were found in the specified folder. Make sure there are .out files in the directory and that the path is correct.")
        sys.exit(1)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    processar_fitxers(fitxers_input, elements, output_directory)
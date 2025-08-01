####################################################################
# Script: qp.py
# Author: Cristina Berga, https://github.com/CristinaBerga/CompChem-Tools/
# Creation Date: June 2025
####################################################################
#
# This Python script is a tool for the automated extraction and
# visualization of Mulliken charges and spin densities from Gaussian
# output files. It's designed to streamline the analysis of
# electronic properties by generating both structured data and
# clear molecular diagrams.
#
####################################################################
#
# Key Features:
#
# 1. Mulliken Data Extraction: Automatically parses the last block
#    of "Mulliken charges and spin densities with hydrogens summed
#    into heavy atoms" from .out or .log files.
#
# 2. Flexible Geometry Input: Can use molecular coordinates from
#    Gaussian output files or from custom .xyz files, which allows
#    for manual definition of bonds.
#
# 3. Enhanced Molecular Visualization: Generates high-quality SVG
#    images of the molecule with the following variations:
#    - "Original": A standard view with no labels.
#    - "Charges": The molecule with Mulliken charges labeled on
#      each heavy atom.
#    - "SpinDensity": The molecule with spin densities labeled on
#      each heavy atom (if available).
#    The script automatically removes C-H bonds for clarity and
#    improves the layout of linear systems (A=B=C).
#
# 4. Data Export: Creates a single Excel file (qp.xlsx) containing
#    a separate sheet for each molecule, with all extracted
#    charge and spin data.
#
####################################################################
#
# Usage:
#
# - Run the script and follow the prompts to specify the source of
#   molecular coordinates (.out/.log or .xyz files).
#
# - Ensure that all corresponding files are in the specified folder(s).
#
# - Requires the installation of Open Babel for format conversion when
#   processing Gaussian output files.
#
####################################################################

import os
import sys
import re
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D
import pandas as pd
from rdkit.Chem import GetPeriodicTable, BondType

# Attempt to import get_valid_sheet_name. This might vary based on openpyxl version.
try:
    from openpyxl.utils import get_valid_sheet_name
except ImportError:
    try:
        from openpyxl.utils.cell import get_valid_sheet_name # Sometimes it's here
    except ImportError:
        # Fallback function if get_valid_sheet_name cannot be imported at all
        
        def get_valid_sheet_name(sheetname):
            """
            Simple fallback for valid Excel sheet names.
            Removes invalid characters and truncates to 31 characters.
            """
            invalid_chars = r'[\[\]\*\?\/\\:\<\>]' # Invalid Excel sheet name characters
            cleaned_name = re.sub(invalid_chars, '', sheetname)
            return cleaned_name[:31]


# --- Constants for Drawing ---
BOND_LINE_WIDTH = 1.0  # Adjust line width for bonds
CUSTOM_BOND_LENGTH = 1.5 # Adjust bond length for 2D layout (higher value = larger molecule)

# --- Functions for processing custom .xyz files ---
def parse_input_file(filename):
    """
    Parses a custom input file (.xyz with bond section)
    Returns lists of atom symbols, coordinates, and bonds.
    """
    ptable = GetPeriodicTable()
    atoms_symbols = []
    coords = []
    bonds = []
    in_bonds = False

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.lower() == "bonds":
                in_bonds = True
                continue
            if not in_bonds:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        atomic_num = int(parts[0])
                        symbol = ptable.GetElementSymbol(atomic_num)
                        x, y, z = map(float, parts[1:4])
                        atoms_symbols.append(symbol)
                        coords.append((x, y, z))
                    except ValueError:
                        print(f"Warning: Wrong atom line in the file '{filename}': {line}. Ignoring.")
                        continue
                else:
                    print(f"Warning: Incomplete atom line in the file '{filename}': {line}. Ignoring.")
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        i, j = int(parts[0]), int(parts[1])
                        bond_type = parts[2] if len(parts) > 2 else "S"
                        bonds.append((i, j, bond_type))
                    except ValueError:
                        print(f"Warning: Wrong bond line in the file '{filename}': {line}. Ignoring.")
                        continue
                else:
                    print(f"Warning: Incomplete bond line in the file '{filename}': {line}. Ignoring.")

    return atoms_symbols, coords, bonds

def filter_bonds(atoms_symbols, bonds):
    """
    Filter links to ensure that indices are within the range of existing atoms.
    """
    max_index = len(atoms_symbols)
    filtered_bonds = []
    for i, j, btype in bonds:
        if 1 <= i <= max_index and 1 <= j <= max_index:
            filtered_bonds.append((i, j, btype))
        else:
            print(f"Warning: Bond ({i}, {j}, {btype}) ignored because it is out of range (there are {max_index} atoms).")
    return filtered_bonds

def build_rdkit_mol_from_parsed_data(atoms_symbols, coords, bonds):
    """
    Constructs an RDKit Mol object from atom symbols, coordinates, and bonds.
    This function also removes hydrogens attached to carbon and sets implicit Hs.
    """
    mol = Chem.RWMol()

    for symbol in atoms_symbols:
        mol.AddAtom(Chem.Atom(symbol))

    bond_type_map = {
        "S": BondType.SINGLE, "1": BondType.SINGLE, "SINGLE": BondType.SINGLE,
        "D": BondType.DOUBLE, "2": BondType.DOUBLE, "DOUBLE": BondType.DOUBLE,
        "T": BondType.TRIPLE, "3": BondType.TRIPLE, "TRIPLE": BondType.TRIPLE,
        "A": BondType.AROMATIC, "AROMATIC": BondType.AROMATIC
    }

    for i, j, btype in bonds:
        i0, j0 = i - 1, j - 1
        
        if btype.upper() == "H":
            print(f"Skipping hydrogen bond creation between {atoms_symbols[i0]} ({i}) and {atoms_symbols[j0]} ({j}).")
            continue

        if 0 <= i0 < mol.GetNumAtoms() and 0 <= j0 < mol.GetNumAtoms():
            if not mol.GetBondBetweenAtoms(i0, j0): # Avoid duplicate bonds
                bt = bond_type_map.get(btype.upper(), BondType.SINGLE)
                mol.AddBond(i0, j0, bt)
        else:
            print(f"Warning: Bond index out of range when building molecule: {i0}-{j0}. Ignoring.")

    # Add 3D conformation
    conf = Chem.Conformer(len(atoms_symbols))
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)

    # Convert to Mol object for further manipulation
    mol = mol.GetMol()

    # Ensure that implicit hydrogens are correctly calculated
    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True) # Disable implicit hydrogen inference to avoid issues

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True) # Attempt to assign stereochemistry

    return mol

def parse_last_mulliken_block(content):
    """
    Take the Mulliken charges and, if present, the spin densities from the last block
    of the Gaussian output file (with hydrogens summed into heavy atoms).
    Returns dictionaries for plotting and lists of tuples for Excel.
    If spin densities are not found, their respective outputs will be None.
    """
    charges_dict = {}
    spins_dict = {}
    q_labels = []
    q_values = []
    p_labels = []
    p_values = []

    # 1. Attempt to find the block with charges and spin densities
    pattern_with_spin = r"Mulliken charges and spin densities with hydrogens summed into heavy atoms:\s*\n\s*1\s+2\n((?:\s*\d+\s+\w+\s+[-\d\.]+\s+[-\d\.]+\n)+)"
    matches_with_spin = re.findall(pattern_with_spin, content)

    if matches_with_spin:
        last_block = matches_with_spin[-1].strip().splitlines()
        for line in last_block:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    atom_num = int(parts[0])
                    atom_type = parts[1]
                    charge = float(parts[2])
                    spin = float(parts[3])
                    label = f"{atom_type}{atom_num}"

                    charges_dict[label] = charge
                    spins_dict[label] = spin
                    q_labels.append(f"q({label})")
                    q_values.append(charge)
                    p_labels.append(f"p({label})")
                    p_values.append(spin)
                except ValueError as e:
                    print(f"Warning: Could not parse Mulliken line with spin: {line.strip()} - Error: {e}")
                    continue
        return charges_dict, spins_dict, (q_labels, q_values), (p_labels, p_values)
    
    else:
        # 2. If spin densities are not found, attempt to find only charges
        pattern_only_charge = r"Mulliken charges with hydrogens summed into heavy atoms:\s*\n\s*1\s*\n((?:\s*\d+\s+\w+\s+[-\d\.]+\n)+)"
        matches_only_charge = re.findall(pattern_only_charge, content)
        
        if matches_only_charge:
            print("DEBUG: No 'Mulliken charges and spin densities' found, searching for 'Mulliken charges'.")
            last_block = matches_only_charge[-1].strip().splitlines()
            for line in last_block:
                parts = line.split()
                if len(parts) >= 3: # 3 parts: atom_num, atom_type, charge
                    try:
                        atom_num = int(parts[0])
                        atom_type = parts[1]
                        charge = float(parts[2])
                        label = f"{atom_type}{atom_num}"

                        charges_dict[label] = charge
                        q_labels.append(f"q({label})")
                        q_values.append(charge)
                        # p_labels and p_values remain empty/None
                    except ValueError as e:
                        print(f"Warning: Could not parse Mulliken line without spin: {line.strip()} - Error: {e}")
                        continue
            print("DEBUG: Mulliken charges found (without spin density).")
            # Returns None for spin data if not found
            return charges_dict, None, (q_labels, q_values), None
        else:
            print("Warning: No Mulliken charge block found (neither with spin nor only charges).")
            return None, None, None, None

def extract_xyz(file_name, xyz_file):
    """
    Extract the 3D coordinates from the last orientation (standard or input)
    preceding the FIRST "Normal termination" found in the file.
    Saves them to a .xyz file.
    """
    with open(file_name, 'r') as f:
        lines = f.readlines()

    # 1. Find the first "Normal termination"
    normal_termination_idx = -1
    for i, line in enumerate(lines):
        if "Normal termination" in line:
            normal_termination_idx = i
            break
    
    if normal_termination_idx == -1:
        print(f"Warning: No 'Normal termination' found in file {file_name}. Geometry will be based on the last orientation found.")
        # If no normal termination, consider the whole file for the last orientation
        search_end_idx = len(lines)
    else:
        search_end_idx = normal_termination_idx

    # 2. From this line backwards, find the last "Standard orientation" or "Input orientation"
    bloc_inici = -1
    for i in range(search_end_idx - 1, -1, -1):
        if "Standard orientation" in lines[i] or "Input orientation" in lines[i]:
            bloc_inici = i
            break
    
    if bloc_inici == -1:
        print(f"Warning: No 'Standard orientation' or 'Input orientation' found in {file_name}.")
        return False

    bloc_dades_inici = bloc_inici + 5
    simbols = Chem.GetPeriodicTable()
    
    geometria_xyz = []
    count = 0
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
            except (ValueError, IndexError):
                continue
        
    if count == 0:
        print(f"Warning: No valid XYZ coordinates found in the selected block of {file_name}.")
        return False

    with open(xyz_file, 'w') as out_xyz:
        out_xyz.write(f"{count}\n\n")
        for sym, x, y, z in geometria_xyz:
            out_xyz.write(f"{sym:2} {x:>12.6f} {y:>12.6f} {z:>12.6f}\n")
    print(f"XYZ coordinates extracted successfully from {file_name}.")
    return True

def convert_xyz_to_mol(base_name, temp_dir):
    """
    Convert a .xyz (standard) file to .mol using Open Babel.
    Temporary .xyz and .mol files are created in temp_dir.
    """
    xyz_file = os.path.join(temp_dir, f"{base_name}.xyz")
    mol_file = os.path.join(temp_dir, f"{base_name}.mol")
    
    if os.path.exists(xyz_file):
        try:
            # -h to add hydrogens (which we will control later), -k to keep atomic order, -f to overwrite
            subprocess.run(["obabel", "-ixyz", xyz_file, "-O", mol_file, "-h", "-k", "-f"], check=True, capture_output=True, text=True)
            print(f"Molecule {mol_file} created successfully using Open Babel.")
            return mol_file
        except FileNotFoundError:
            print(f"ERROR: Open Babel (obabel o obabel.exe) was not found. Be sure it is installed and in your system PATH.")
            return None
        except subprocess.CalledProcessError as e:
            print(f"ERROR: Open Babel failed for {xyz_file}. Exit code: {e.returncode}")
            print(f"Open Babel output (stdout): {e.stdout.strip()}")
            print(f"Open Babel output (stderr): {e.stderr.strip()}")
            print("Make sure the Open Babel command is correct for your version.")
            return None
    return None

def remove_h_except_on_carbon(mol):
    """
    Eliminate hydrogens that are attached to a carbon atom.
    Keep hydrogens attached to other elements (O, N, S, P, etc.).
    """
    rwmol = Chem.RWMol(mol)
    
    h_to_remove_indices = []
    
    for i in range(rwmol.GetNumAtoms()):
        atom = rwmol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 1: # Is hydrogen?
            for neighbor_atom in atom.GetNeighbors():
                if neighbor_atom.GetAtomicNum() == 6: # The neighbor is Carbon?
                    h_to_remove_indices.append(atom.GetIdx())
                    break
                    
    for h_idx in sorted(h_to_remove_indices, reverse=True):
        rwmol.RemoveAtom(h_idx)
    
    return rwmol.GetMol()

def draw_molecule_with_labels(mol_original, data_labels, title, file_name, label_format="{:.2f}", output_dir="."):
    """
    Generate the SVG image of the molecule with the provided labels (charges or spins),
    removing the hydrogens attached to carbon and forcing linearity in A=B=C.
    If title is "Original" and data_labels is empty, do not show the carbon symbol.
    """
    Chem.RemoveStereochemistry(mol_original)

    mol_to_draw = remove_h_except_on_carbon(mol_original)
    
    try:
        Chem.SanitizeMol(mol_to_draw)
    except Exception as e:
        print(f"WARNING: Could not sanitize the molecule (with H not-C) before 2D generation: {e}")
        print("This may affect the quality of the 2D layout.")

    AllChem.Compute2DCoords(mol_to_draw, clearConfs=True)

    # --- Force linearity for chains A=B=C ---
    mol_conf = mol_to_draw.GetConformer()
    
    for atom_B in mol_to_draw.GetAtoms():
        double_bonds_of_B = [bond for bond in atom_B.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE]
        
        if len(double_bonds_of_B) == 2:
            atom_A = double_bonds_of_B[0].GetOtherAtom(atom_B)
            atom_C = double_bonds_of_B[1].GetOtherAtom(atom_B)
            
            if atom_A.GetIdx() != atom_C.GetIdx():
                idx_A = atom_A.GetIdx()
                idx_B = atom_B.GetIdx()
                idx_C = atom_C.GetIdx()

                pos_A = mol_conf.GetAtomPosition(idx_A)
                pos_B = mol_conf.GetAtomPosition(idx_B)
                pos_C = mol_conf.GetAtomPosition(idx_C)

                vec_AC_x = pos_C.x - pos_A.x
                vec_AC_y = pos_C.y - pos_A.y
                
                len_AC_sq = vec_AC_x**2 + vec_AC_y**2
                
                if len_AC_sq > 1e-6:
                    # Project B onto the line AC
                    vec_AB_x = pos_B.x - pos_A.x
                    vec_AB_y = pos_B.y - pos_A.y
                    
                    dot_product = vec_AB_x * vec_AC_x + vec_AB_y * vec_AC_y
                    t = dot_product / len_AC_sq

                    new_B_x = pos_A.x + t * vec_AC_x
                    new_B_y = pos_A.y + t * vec_AC_y

                    mol_conf.SetAtomPosition(idx_B, Point3D(new_B_x, new_B_y, 0))

    mol_to_draw.SetProp("CustomBondLength", str(CUSTOM_BOND_LENGTH))

    width, height = 2000, 1800 # Adjust to give more space and avoid overlaps
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()

    opts.addAtomIndices = False
    opts.bondLineWidth = BOND_LINE_WIDTH
    opts.addStereoAnnotation = False
    opts.clearBackground = True

    # Save original state of radical electrons
    original_radical_electrons_state = {}
    for atom_idx in range(mol_to_draw.GetNumAtoms()):
        atom = mol_to_draw.GetAtomWithIdx(atom_idx)
        original_radical_electrons_state[atom_idx] = atom.GetNumRadicalElectrons()
        if original_radical_electrons_state[atom_idx] > 0:
            atom.SetNumRadicalElectrons(0) # Temporarily hide radical points

    # Build the mapping original_index -> mol_to_draw_index
    original_idx_to_draw_idx = {}
    draw_atom_counter = 0

    # Iterate over the atoms of the original molecule to create the mapping
    # This is crucial because mol_to_draw has fewer atoms (H on carbon are removed)
    for i, original_atom in enumerate(mol_original.GetAtoms()):
        is_h_on_carbon = False
        if original_atom.GetAtomicNum() == 1: # Is hydrogen?
            for neighbor in original_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6: # The neighbor is Carbon?
                    is_h_on_carbon = True
                    break

        if not is_h_on_carbon: # If the atom is not an H attached to a Carbon (it stays)
            # Find its corresponding atom in mol_to_draw
            # This is based on RDKit maintaining the order of non-removed atoms
            if draw_atom_counter < mol_to_draw.GetNumAtoms():
                original_idx_to_draw_idx[i] = draw_atom_counter
                draw_atom_counter += 1
            else:
                print(f"LOGICAL ERROR: Drawing atom {original_atom.GetSymbol()}{i+1} out of range of mol_to_draw. This may indicate a problem with remove_h_except_on_carbon or with the indexing.")

    # Iterate over the atoms of the ORIGINAL molecule to associate the Mulliken data
    for i, atom_original in enumerate(mol_original.GetAtoms()):
        label_key = f"{atom_original.GetSymbol()}{i+1}"
        
        if i in original_idx_to_draw_idx:
            new_draw_idx = original_idx_to_draw_idx[i]

            # Add charge/spin label
            if data_labels and label_key in data_labels: # Ensure data_labels is not None and the key exists
                value = data_labels[label_key]
                formatted_label = label_format.format(value)

                # Check if a label already exists for this atom index
                if new_draw_idx in opts.atomLabels:
                    current_label = opts.atomLabels[new_draw_idx]
                else:
                    current_label = atom_original.GetSymbol()

                # If the current label is just the atom symbol, replace it. Otherwise, append it.
                if current_label == atom_original.GetSymbol():
                    opts.atomLabels[new_draw_idx] = formatted_label
                else:
                    opts.atomLabels[new_draw_idx] = f"{current_label}\n{formatted_label}"
            else:
                # Handle cases where no specific label is provided, or for the "Original" title
                if title == "Original" and atom_original.GetAtomicNum() == 6:
                    opts.atomLabels[new_draw_idx] = "" # Hide the carbon symbol for the "Original" drawing
                elif new_draw_idx not in opts.atomLabels: # Ensure we don't overwrite if already set
                    opts.atomLabels[new_draw_idx] = atom_original.GetSymbol()


    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol_to_draw)
    drawer.FinishDrawing()

    output_path = os.path.join(output_dir, f"{file_name}_{title}.svg")
    with open(output_path, "w") as f:
        f.write(drawer.GetDrawingText())
    print(f"Image of the molecule with {title} saved to: {output_path}")

    # Restore the radical electrons to their original state
    for atom_idx in range(mol_to_draw.GetNumAtoms()):
        atom = mol_to_draw.GetAtomWithIdx(atom_idx)
        if atom_idx in original_radical_electrons_state:
            atom.SetNumRadicalElectrons(original_radical_electrons_state[atom_idx])

def write_molecule_data_to_excel(charges_lists, spins_lists, filename_base, excel_writer):
    """
    Write the charge and spin density data of a molecule to a new sheet in the Excel file.
    """
    if charges_lists is None:
        print(f"WARNING: No charge data found for Excel for '{filename_base}'.")
        return

    q_labels, q_values = charges_lists

    data = {'q(atom)': q_labels, 'Charges': q_values}
    
    if spins_lists is not None:
        p_labels, p_values = spins_lists
        # Add an empty column for separation if both charges and spins are present
        data[''] = [''] * len(q_labels) 
        data['p(atom)'] = p_labels
        data['Spin density'] = p_values
    else:
        print(f"WARNING: No spin data found for Excel for '{filename_base}'. Only charges will be written.")

    df = pd.DataFrame(data)

    # Use the imported or fallback get_valid_sheet_name
    sheetname = get_valid_sheet_name(filename_base[:31]) 

    try:
        df.to_excel(excel_writer, sheet_name=sheetname, index=False)
        print(f"Data for '{filename_base}' written to Excel sheet '{sheetname}'.")
    except ValueError as e:
        print(f"ERROR: Could not write Excel data for '{filename_base}'. The sheet name '{sheetname}' may be invalid or too long. Error: {e}")
    except Exception as e:
        print(f"ERROR: Could not write Excel data for '{filename_base}'. General error: {e}")


def processar_fitxers_mulliken(llista_fitxers, source_type, xyz_folder_path=None):
    """
    Process output files from Gaussian to extract charges and (optionally) spins,
    generates an Excel file with the data and separate SVG images for each.
    All output files are saved within the input folder.
    Allows choosing the source of the coordinates (from the .out/.log file or from an external .xyz).
    """
    if not llista_fitxers:
        print("ERROR: No files were provided to process.")
        return

    fitxers = sorted(llista_fitxers)

    base_output_folder = os.path.dirname(fitxers[0])

    output_images_dir = os.path.join(base_output_folder, "qp_images")
    os.makedirs(output_images_dir, exist_ok=True)
    print(f"Output images folder: {os.path.abspath(output_images_dir)}")

    output_excel_file = os.path.join(base_output_folder, "qp.xlsx")
    print(f"Creating Excel file: {os.path.abspath(output_excel_file)}")

    temp_files_dir = os.path.join(base_output_folder, "temp_files")
    os.makedirs(temp_files_dir, exist_ok=True)
    print(f"Creating temporary files folder: {os.path.abspath(temp_files_dir)}")

    try:
        with pd.ExcelWriter(output_excel_file, engine='openpyxl') as writer:
            for fitxer in fitxers:
                base = os.path.splitext(os.path.basename(fitxer))[0]
                
                print(f"\nProcessing file: {fitxer}")

                mol = None # Initialize mol to None for each iteration

                if source_type == 'out':
                    xyz_file = os.path.join(temp_files_dir, f"{base}.xyz")
                    mol_file = os.path.join(temp_files_dir, f"{base}.mol")
                    
                    extret = extract_xyz(fitxer, xyz_file)
                    if not extret:
                        print(f"Warning: Could not extract .xyz file from {fitxer}. Skipping processing.")
                        if os.path.exists(xyz_file): os.remove(xyz_file)
                        continue

                    mol_path_from_obabel = convert_xyz_to_mol(base, temp_files_dir)
                    if not mol_path_from_obabel:
                        print(f"Warning: Could not create .mol file for {base}. Skipping processing.")
                        if os.path.exists(xyz_file): os.remove(xyz_file)
                        if os.path.exists(mol_file): os.remove(mol_file)
                        continue

                    mol = Chem.MolFromMolFile(mol_path_from_obabel, removeHs=False, sanitize=False)
                    if mol is None:
                        print(f"Warning: RDKit could not load the .mol file ({mol_path_from_obabel}). Skipping image generation.")
                        if os.path.exists(xyz_file): os.remove(xyz_file)
                        if os.path.exists(mol_file): os.remove(mol_file)
                        continue

                    # Clean up temporary files for the current molecule (if generated here)
                    if os.path.exists(xyz_file): os.remove(xyz_file)
                    if os.path.exists(mol_file): os.remove(mol_file)
                    print(f"Temporary .xyz and .mol files for {base} deleted (if they were used).")


                elif source_type == 'xyz':
                    if xyz_folder_path is None:
                        print("ERROR: The path to the .xyz folder has not been provided.")
                        continue
                    
                    xyz_input_file = os.path.join(xyz_folder_path, f"{base}.xyz")
                    if not os.path.exists(xyz_input_file):
                        print(f"Warning: The .xyz file '{xyz_input_file}' was not found for '{base}'. Skipping processing.")
                        continue
                    
                    try:
                        atoms_symbols, coords, bonds = parse_input_file(xyz_input_file)
                        if not atoms_symbols:
                            print(f"Warning: No atoms could be parsed from the file '{xyz_input_file}'. Skipping processing.")
                            continue
                        
                        filtered_bonds = filter_bonds(atoms_symbols, bonds)
                        mol = build_rdkit_mol_from_parsed_data(atoms_symbols, coords, filtered_bonds)
                        print(f"RDKit molecule successfully constructed from {xyz_input_file}.")

                    except Exception as e:
                        print(f"ERROR: Could not process the .xyz file '{xyz_input_file}': {e}. Skipping processing.")
                        continue
                
                if mol is None:
                    print(f"Warning: Could not generate RDKit molecule for '{base}'. Skipping drawing and Excel export.")
                    continue

                # Charges and spins are always extracted from the .out/.log file
                with open(fitxer, 'r') as f:
                    content = f.read()
                
                charges_dict, spins_dict, charges_lists, spins_lists = parse_last_mulliken_block(content)

                if charges_dict:
                    write_molecule_data_to_excel(charges_lists, spins_lists, base, writer)
                    draw_molecule_with_labels(mol, charges_dict, "Charges", base, output_dir=output_images_dir)
                    
                    if spins_dict:
                        draw_molecule_with_labels(mol, spins_dict, "SpinDensity", base, output_dir=output_images_dir)
                    else:
                        print(f"Warning: No spin density found for {base}. Spin density image will not be generated.")
                else:
                    print(f"Warning: No valid charge data found in {fitxer}. No image or Excel entry will be generated for this molecule.")

                draw_molecule_with_labels(mol, {}, "Original", base, output_dir=output_images_dir) # Always draw the original

        # Clean up the temporary directory if it's empty
        if os.path.exists(temp_files_dir) and not os.listdir(temp_files_dir):
            os.rmdir(temp_files_dir)
            print(f"Temporary files directory '{temp_files_dir}' has been removed.")

        print(f"Excel file '{output_excel_file}' completed and saved to: {os.path.abspath(output_excel_file)}")
    except Exception as e:
        print(f"FATAL ERROR: Could not create or write to Excel file '{output_excel_file}'. Make sure the file is not open and you have write permission. Error: {e}")

    print(f"Processing complete. SVG images saved to folder '{output_images_dir}'.")

if __name__ == "__main__":

    print("Where do you want to extract the coordinates from?")
    print("1. .out/.log files")
    print("2. .xyz files")

    choice = input("1 or 2: ").strip()

    input_files = []
    folder_path = ""
    xyz_folder_path = None
    source_type = None

    if choice == '1':
        source_type = 'out'
        if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
            folder_path = sys.argv[1]
            input_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".out") or f.endswith(".log")]
        else:
            folder_path = input("Enter the path to the folder with .out or .log files: ").strip()
            if not os.path.isdir(folder_path):
                print(f"ERROR: The specified path '{folder_path}' does not exist or is not a valid folder.")
                sys.exit(1)
            input_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".out") or f.endswith(".log")]
        
        if not input_files:
            print("No .out or .log files found to process.")
            sys.exit(0)

    elif choice == '2':
        source_type = 'xyz'
        # For the .xyz option, we still need a list of .out/.log files
        # to extract the charges, and to know which .xyz files to look for.
        if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
            folder_path = sys.argv[1]
            input_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".out") or f.endswith(".log")]
        else:
            folder_path = input("Enter the path to the folder with .out or .log files: ").strip()
            if not os.path.isdir(folder_path):
                print(f"ERROR: The specified path '{folder_path}' does not exist or is not a valid folder.")
                sys.exit(1)
            input_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".out") or f.endswith(".log")]
        
        if not input_files:
            print("No .out or .log files found to process.")
            sys.exit(0)

        xyz_folder_path = input("Enter the path to the folder with the corresponding .xyz files: ").strip()
        if not os.path.isdir(xyz_folder_path):
            print(f"ERROR: The specified path '{xyz_folder_path}' does not exist or is not a valid folder.")
            sys.exit(1)

    else:
        print("Invalid option. Please enter 1 or 2.")
        sys.exit(1)

    # Call the main function with the selected source type
    processar_fitxers_mulliken(input_files, source_type, xyz_folder_path)
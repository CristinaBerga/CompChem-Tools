####################################################################
# Script: xyz.py
# Author: Cristina Berga, https://github.com/CristinaBerga/CompChem-Tools/
# Creation Date: June 2025
####################################################################
#
# This Python script is a molecular visualization tool that
# processes a custom .xyz file format. The script uses RDKit to
# build the molecule and generate a clear SVG image.
#
####################################################################
#
# Key Features:
#
# 1. Custom File Parsing: Reads a specialized .xyz file that
#    contains a list of atoms with their 3D coordinates, followed
#    by a "bonds" section to explicitly define atom connections
#    and their types.
#    The input.xyz is the one that you get when you select "Save 
#    bonding information" in Chemcraft.
#
# 2. Smart Hydrogen Handling: For visualization, the script
#    automatically removes hydrogen atoms that are bonded to
#    carbon, but keeps hydrogens attached to other atoms (e.g.,
#    oxygen, nitrogen) to provide a clearer, more informative
#    structural image.
#
# 3. SVG Visualization: Generates a high-quality, scalable
#    SVG image of the molecule, which is saved in the same
#    directory as the input file.
#
####################################################################
#
# Usage:
#
# - Run the script and provide the name of your custom .xyz
#   input file at the prompt.
#
# - The output image (.svg) will be created with the same name
#   as the input file in the same directory.
#
####################################################################

import os
from rdkit.Chem import GetPeriodicTable
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem

def parse_input_file(filename):
    ptable = GetPeriodicTable()
    atoms = []
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
                    atomic_num = int(parts[0])
                    symbol = ptable.GetElementSymbol(atomic_num)
                    x, y, z = map(float, parts[1:4])
                    atoms.append(symbol)
                    coords.append((x, y, z))
                else:
                    print(f"Atom line wrong format: {line}")
            else:
                parts = line.split()
                if len(parts) >= 2:
                    i, j = int(parts[0]), int(parts[1])
                    bond_type = parts[2] if len(parts) > 2 else "S"
                    bonds.append((i, j, bond_type))
                else:
                    print(f"Bond line wrong format: {line}")

    return atoms, coords, bonds

def filter_bonds(atoms, bonds):
    max_index = len(atoms)
    filtered_bonds = []
    for i, j, btype in bonds:
        if 1 <= i <= max_index and 1 <= j <= max_index:
            filtered_bonds.append((i, j, btype))
        else:
            print(f"Bond ({i}, {j}, {btype}) ignored because it is out of range (there are {max_index} atoms).")
    return filtered_bonds

def build_rdkit_mol(atoms, coords, bonds):
    mol = Chem.RWMol()

    for symbol in atoms:
        mol.AddAtom(Chem.Atom(symbol))

    bond_type_map = {
        "S": BondType.SINGLE,
        "1": BondType.SINGLE,
        "SINGLE": BondType.SINGLE,
        "D": BondType.DOUBLE,
        "2": BondType.DOUBLE,
        "DOUBLE": BondType.DOUBLE,
        "T": BondType.TRIPLE,
        "3": BondType.TRIPLE,
        "TRIPLE": BondType.TRIPLE,
    }

    for i, j, btype in bonds:
        i0, j0 = i - 1, j - 1
        
        if btype.upper() == "H":
            print(f"Skipping hydrogen bond creation between {atoms[i0]} ({i}) and {atoms[j0]} ({j}).")
            continue

        if 0 <= i0 < mol.GetNumAtoms() and 0 <= j0 < mol.GetNumAtoms():
            bt = bond_type_map.get(btype.upper(), BondType.SINGLE)
            mol.AddBond(i0, j0, bt)

    conf = Chem.Conformer(len(atoms))
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf)

    mol = mol.GetMol()

    atoms_to_keep = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "H":
            atoms_to_keep.append(atom.GetIdx())
        else:
            is_bonded_to_carbon = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "C":
                    is_bonded_to_carbon = True
                    break
            if not is_bonded_to_carbon:
                atoms_to_keep.append(atom.GetIdx())

    emol = Chem.EditableMol(Chem.Mol())

    old_to_new = {}
    for new_idx, old_idx in enumerate(atoms_to_keep):
        atom = mol.GetAtomWithIdx(old_idx)
        emol.AddAtom(Chem.Atom(atom.GetSymbol()))
        old_to_new[old_idx] = new_idx

    for bond in mol.GetBonds():
        i_old = bond.GetBeginAtomIdx()
        j_old = bond.GetEndAtomIdx()
        if i_old in old_to_new and j_old in old_to_new:
            i_new = old_to_new[i_old]
            j_new = old_to_new[j_old]
            emol.AddBond(i_new, j_new, bond.GetBondType())

    new_mol = emol.GetMol()

    new_conf = Chem.Conformer(len(atoms_to_keep))
    for new_idx, old_idx in enumerate(atoms_to_keep):
        pos = mol.GetConformer().GetAtomPosition(old_idx)
        new_conf.SetAtomPosition(new_idx, pos)
    new_mol.AddConformer(new_conf)

    for atom in new_mol.GetAtoms():
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(atom.GetTotalNumHs())

    new_mol.UpdatePropertyCache(strict=False)

    return new_mol

def draw_molecule(mol, filename):
    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(1500, 1500)
    opts = drawer.drawOptions()

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    img_content = drawer.GetDrawingText()
    with open(filename, "w") as img_file:
        img_file.write(img_content)

    print(f"Image of the molecule saved to: {filename}")

if __name__ == "__main__":
    filename = input("Enter the name of the file with coordinates and bonds: ").strip()
    atoms, coords, bonds = parse_input_file(filename)
    print(f"Atoms read: {len(atoms)}")
    print(f"Bonds read: {len(bonds)}")

    bonds = filter_bonds(atoms, bonds)
    mol = build_rdkit_mol(atoms, coords, bonds)
    
    img_filename = os.path.splitext(filename)[0] + ".svg" 
    draw_molecule(mol, img_filename)

# -*- coding: utf-8 -*-

####################################################################
# Script: multigap.py
# Author: Cristina Berga, https://github.com/CristinaBerga
# Creation Date: May 2025
####################################################################
#
# This script for Python 2 and 3 is a tool designed to automate
# the extraction and analysis of key data from Gaussian generated
# outputs.
#
####################################################################
#
# Key Features:
#
# 1. Extraction of HOMO and LUMO energies (in Hartree) and their indices.
#
# 2. HOMO-LUMO Gap in kcal/mol.
#
# 3. Extraction of the Total Energy (in Hartree).
#
# 4. All extracted data is exported to a CSV file
#    (multigap.csv), using the semicolon as a delimiter.
#
####################################################################
#
# Usage:
#
# - To process .out files in the current directory:
#   python multigap.py
#
# - To process .out files in a specific folder:
#   python multigap.py folder
#
# The output file `multigap.csv` will be generated in the directory
# where the script is executed.
#
####################################################################

from __future__ import print_function
import os
import sys
import csv

PY2 = sys.version_info[0] == 2

def imprimeix(text):
    """Function to print text compatible with Python 2 and 3."""
    if PY2:
        sys.stdout.write(text.encode('utf-8') + b"\n")
    else:
        print(text)

def extreure_homo_lumo_i_scf(nom_fitxer):
    """
    Extracts HOMO, LUMO values, their indices, and SCF energy from a Gaussian output file.
    Returns: homo, homo_index, lumo, lumo_index, energia_total
    """
    homo = None
    lumo = None
    homo_index = None
    lumo_index = None
    energia_total = None

    with open(nom_fitxer, 'r') as f:
        lines = f.readlines()

    for line in reversed(lines):
        if "SCF Done" in line:
            if "=" in line:
                energia_str = line.split("=")[1].strip().split()[0]
                try:
                    energia_total = float(energia_str)
                except ValueError:
                    energia_total = None
                break

    occ_indices = [i for i, l in enumerate(lines) if "Alpha  occ. eigenvalues" in l or "Beta  occ. eigenvalues" in l]
    virt_indices = [i for i, l in enumerate(lines) if "Alpha virt. eigenvalues" in l or "Beta virt. eigenvalues" in l]

    if not occ_indices or not virt_indices:
        return None, None, None, None, energia_total

    def trobar_bloc(indices):
        """Finds the continuous block of indices (for multiline eigenvalues)."""
        blocs = []
        if not indices:
            return []
        bloc_actual = [indices[0]]
        for idx in indices[1:]:
            if idx == bloc_actual[-1] + 1:
                bloc_actual.append(idx)
            else:
                blocs.append(bloc_actual)
                bloc_actual = [idx]
        blocs.append(bloc_actual)
        return blocs[-1]

    ult_bloc_occ = trobar_bloc(occ_indices)
    ult_bloc_virt = trobar_bloc(virt_indices)

    occ_vals = []
    for i in ult_bloc_occ:
        parts = lines[i].split('--')
        if len(parts) < 2:
            vals = parts[0].split()
        else:
            vals = parts[-1].split()
        occ_vals.extend([float(v) for v in vals])

    virt_vals = []
    for i in ult_bloc_virt:
        parts = lines[i].split('--')
        if len(parts) < 2:
            vals = parts[0].split()
        else:
            vals = parts[-1].split()
        virt_vals.extend([float(v) for v in vals])

    if occ_vals:
        homo = occ_vals[-1]
        homo_index = len(occ_vals)

    if virt_vals:
        lumo = virt_vals[0]
        lumo_index = homo_index + 1 if homo_index is not None else None

    return homo, homo_index, lumo, lumo_index, energia_total


def processar_fitxers_carpeta(carpeta, output_csv="multigap.csv"):
    """
    Processes all .out files in a folder and writes the data to a CSV.
    """
    fitxers = [f for f in os.listdir(carpeta) if f.endswith('.out')]
    fitxers.sort()

    if PY2:
        csvfile = open(output_csv, 'wb')
        writer = csv.writer(csvfile, delimiter=b';')
    else:
        csvfile = open(output_csv, 'w', encoding='utf-8', newline='')
        writer = csv.writer(csvfile, delimiter=';')

    camp_encap = [u'System', u'HOMO (Hartree)', u'Orbital', u'LUMO (Hartree)', u'Orbital', u'Gap (kcal/mol)', u'Electronic Energy (Hartree)']

    if PY2:
        encoded_header = [s.encode('utf-8') for s in camp_encap]
        writer.writerow(encoded_header)
    else:
        writer.writerow(camp_encap)

    for fitxer in fitxers:
        ruta_completa = os.path.join(carpeta, fitxer)
        homo, homo_index, lumo, lumo_index, energia_total = extreure_homo_lumo_i_scf(ruta_completa)

        gap = None
        if homo is not None and lumo is not None:
            gap = (lumo - homo) * 627.5095

        nom_sense_extensio = fitxer[:-4] if fitxer.endswith('.out') else fitxer

        row_data = [
            nom_sense_extensio,
            "%.6f" % homo if homo is not None else "",
            str(homo_index) if homo_index is not None else "",
            "%.6f" % lumo if lumo is not None else "",
            str(lumo_index) if lumo_index is not None else "",
            "%.2f" % gap if gap is not None else "",
            "%.6f" % energia_total if energia_total is not None else ""
        ]

        if PY2:
            encoded_row_data = [unicode(item).encode('utf-8') for item in row_data]
            writer.writerow(encoded_row_data)
        else:
            writer.writerow(row_data)

    csvfile.close()

    if PY2:
        imprimeix(u"Success in '" + unicode(output_csv) + u"'")
    else:
        imprimeix("Success in '{}'".format(output_csv))


if __name__ == "__main__":
    carpeta = "."
    if len(sys.argv) == 2:
        carpeta = sys.argv[1]

    processar_fitxers_carpeta(carpeta)
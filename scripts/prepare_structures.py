#!/usr/bin/env python3

from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]

PDB_DIR = PROJECT / "data" / "pdb"
OUT_DIR = PROJECT / "data" / "processed"

SELECTED_CHAINS = {"A", "B", "C", "D"}

# Convention commune aux deux structures :
#
# chaînes alpha A/C : 141 résidus protéiques -> HEM 142
# chaînes beta  B/D : 146 résidus protéiques -> HEM 147
#
# Cela évite que les quatre HEM de 2HBS soient interprétés
# comme un seul résidu après fusion des chaînes.
HEME_RESSEQ = {
    "A": 142,
    "B": 147,
    "C": 142,
    "D": 147,
}


def normalize_hetatm_heme(line):
    """Normalise le numéro de résidu des HEM dans une ligne HETATM."""

    residue = line[17:20].strip()
    chain = line[21]

    if residue == "HEM" and chain in HEME_RESSEQ:
        new_resseq = HEME_RESSEQ[chain]

        # Colonnes PDB 23-26 = numéro du résidu
        line = line[:22] + f"{new_resseq:4d}" + line[26:]

    return line


def normalize_link_heme(line):
    """Normalise les numéros HEM éventuellement présents dans une ligne LINK."""

    # Premier résidu du LINK
    resname1 = line[17:20].strip()
    chain1 = line[21] if len(line) > 21 else ""

    if resname1 == "HEM" and chain1 in HEME_RESSEQ:
        new_resseq = HEME_RESSEQ[chain1]
        line = line[:22] + f"{new_resseq:4d}" + line[26:]

    # Deuxième résidu du LINK
    resname2 = line[47:50].strip()
    chain2 = line[51] if len(line) > 51 else ""

    if resname2 == "HEM" and chain2 in HEME_RESSEQ:
        new_resseq = HEME_RESSEQ[chain2]

        # Colonnes PDB 53-56 = numéro du second résidu
        line = line[:52] + f"{new_resseq:4d}" + line[56:]

    return line


def prepare_structure(source_name, output_name):
    source = PDB_DIR / source_name
    output = OUT_DIR / output_name

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    link_lines = []
    coordinate_lines = []

    n_atom = 0
    n_heme = 0
    n_ter = 0

    with source.open() as handle:
        for line in handle:
            record = line[:6].strip()

            if record == "LINK":
                chain1 = line[21] if len(line) > 21 else ""
                chain2 = line[51] if len(line) > 51 else ""

                if chain1 in SELECTED_CHAINS and chain2 in SELECTED_CHAINS:
                    line = normalize_link_heme(line)
                    link_lines.append(line)

            elif record == "ATOM":
                chain = line[21]

                if chain in SELECTED_CHAINS:
                    coordinate_lines.append(line)
                    n_atom += 1

            elif record == "HETATM":
                residue = line[17:20].strip()
                chain = line[21]

                if chain in SELECTED_CHAINS and residue == "HEM":
                    line = normalize_hetatm_heme(line)
                    coordinate_lines.append(line)
                    n_heme += 1

            elif record == "TER":
                chain = line[21] if len(line) > 21 else ""

                if chain in SELECTED_CHAINS:
                    coordinate_lines.append(line)
                    n_ter += 1

    with output.open("w") as handle:
        handle.write(
            f"REMARK   Generated from {source_name} "
            f"for HBB molecular dynamics project\n"
        )
        handle.write(
            "REMARK   Selected biological tetramer: chains A B C D\n"
        )
        handle.write(
            "REMARK   Crystallographic waters removed; HEM retained\n"
        )
        handle.write(
            "REMARK   HEM residue numbering normalized per globin chain\n"
        )

        for line in link_lines:
            handle.write(line)

        for line in coordinate_lines:
            handle.write(line)

        handle.write("END\n")

    print(f"{source_name} -> {output}")
    print(f"  Protein ATOM records : {n_atom}")
    print(f"  HEM HETATM records   : {n_heme}")
    print(f"  TER records          : {n_ter}")
    print(f"  LINK records         : {len(link_lines)}")


prepare_structure("2DN2.pdb", "WT_clean.pdb")
prepare_structure("2HBS.pdb", "HBS_clean.pdb")

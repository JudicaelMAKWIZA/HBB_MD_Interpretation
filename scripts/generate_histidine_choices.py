#!/usr/bin/env python3

import re
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]

WT_LOG = (
    PROJECT
    / "tests"
    / "pdb2gmx_wt_gromos54a7_specbond"
    / "pdb2gmx_wt_specbond.log"
)

OUTPUT = PROJECT / "tests" / "histidine_choices_wt_reference.txt"

EXPECTED = {
    "A": [20, 45, 50, 58, 72, 87, 89, 103, 112, 122],
    "B": [2, 63, 77, 92, 97, 116, 117, 143, 146],
    "C": [20, 45, 50, 58, 72, 87, 89, 103, 112, 122],
    "D": [2, 63, 77, 92, 97, 116, 117, 143, 146],
}

# Histidines proximales directement coordonnées au fer de l'hème
HEME_HIS = {
    ("A", 87),
    ("B", 92),
    ("C", 87),
    ("D", 92),
}

# Correspondance avec le menu GROMOS54A7 observé :
# 0 = HISA = proton ND1 = HISD
# 1 = HISB = proton NE2 = HISE
# 2 = HISH = doublement protonée
# 3 = HIS1 = couplée à l'hème
STATE_TO_OPTION = {
    "HISD": 0,
    "HISE": 1,
    "HISH": 2,
}

pattern = re.compile(
    r"Will use (HIS[DHE]) for residue\s+(\d+)"
)

values = []

for line in WT_LOG.read_text().splitlines():
    match = pattern.search(line)

    if match:
        values.append(
            (int(match.group(2)), match.group(1))
        )

expected_total = sum(len(v) for v in EXPECTED.values())

if len(values) != expected_total:
    raise RuntimeError(
        f"Nombre d'histidines inattendu : "
        f"{len(values)} au lieu de {expected_total}"
    )

choices = []
position = 0

print("=== CHOIX HISTIDINES : PROFIL WT DE RÉFÉRENCE ===")

for chain, expected_residues in EXPECTED.items():

    block = values[position:position + len(expected_residues)]

    observed = [resid for resid, state in block]

    if observed != expected_residues:
        raise RuntimeError(
            f"Ordre inattendu pour chaîne {chain}: {observed}"
        )

    for resid, state in block:

        if (chain, resid) in HEME_HIS:
            option = 3
            final_state = "HIS1"
        else:
            option = STATE_TO_OPTION[state]
            final_state = state

        choices.append(str(option))

        print(
            f"Chaîne {chain:1s}  HIS {resid:3d}  "
            f"{state:4s} -> option {option} ({final_state})"
        )

    position += len(expected_residues)

OUTPUT.write_text("\n".join(choices) + "\n")

print()
print(f"Nombre de choix : {len(choices)}")
print(f"Fichier créé    : {OUTPUT}")

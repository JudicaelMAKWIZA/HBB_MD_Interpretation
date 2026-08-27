#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WT="$ROOT/workflow/wt_localff"
HBS="$ROOT/workflow/hbs_localff"
MDP="$ROOT/mdp"
RESULTS="$ROOT/results"

THREADS="${THREADS:-8}"

mkdir -p "$RESULTS"

STATUS="$RESULTS/FINAL_MD_STATUS.txt"
: > "$STATUS"

log()
{
    echo "$*" | tee -a "$STATUS"
}

fail()
{
    log ""
    log "PIPELINE ARRÊTÉ : $*"
    exit 1
}


# ============================================================
# 1. VALIDATION DES INDEX
# ============================================================

validate_index()
{
    local DIR="$1"
    local NAME="$2"

    python3 - "$DIR/em.gro" "$DIR/index_nvt.ndx" "$NAME" <<'PY'
import sys

gro_path, ndx_path, name = sys.argv[1:]

lines = open(gro_path).read().splitlines()
natoms = int(lines[1])
atom_lines = lines[2:2+natoms]

expected_protein = set()
expected_solvent = set()

for i, line in enumerate(atom_lines, 1):
    resname = line[5:10].strip()

    if resname in {"SOL", "NA", "CL"}:
        expected_solvent.add(i)
    else:
        expected_protein.add(i)


groups = {}
current = None

for line in open(ndx_path):
    s = line.strip()

    if not s:
        continue

    if s.startswith("[") and s.endswith("]"):
        current = s[1:-1].strip()
        groups[current] = []
    elif current:
        groups[current].extend(map(int, s.split()))


protein = set(groups.get("Protein_HEM", []))
solvent = set(groups.get("Water_and_ions", []))

if protein != expected_protein:
    raise SystemExit(
        f"{name}: groupe Protein_HEM incorrect"
    )

if solvent != expected_solvent:
    raise SystemExit(
        f"{name}: groupe Water_and_ions incorrect"
    )

if protein & solvent:
    raise SystemExit(
        f"{name}: chevauchement entre groupes"
    )

if len(protein | solvent) != natoms:
    raise SystemExit(
        f"{name}: tous les atomes ne sont pas couverts"
    )

print(
    f"{name}: INDEX OK | total={natoms} | "
    f"Protein_HEM={len(protein)} | "
    f"Water_and_ions={len(solvent)}"
)
PY
}


log "============================================================"
log "VALIDATION DES INDEX"
log "============================================================"

validate_index "$WT" "WT" | tee -a "$STATUS"
validate_index "$HBS" "HbS" | tee -a "$STATUS"


# ============================================================
# 2. FICHIERS MDP
# ============================================================

cat > "$MDP/nvt.mdp" <<'MDPFILE'
; ============================================================
; NVT - 100 ps
; HBB WT / HbS
; ============================================================

define                  = -DPOSRES

integrator              = md
dt                      = 0.002
nsteps                  = 50000

continuation            = no

constraints             = all-bonds
constraint-algorithm    = lincs
lincs-iter              = 1
lincs-order             = 4

cutoff-scheme           = Verlet
nstlist                 = 10

coulombtype             = Reaction-Field
rcoulomb                = 1.4
epsilon-r               = 1
epsilon-rf              = 61

vdwtype                 = Cut-off
rvdw                    = 1.4
DispCorr                = no

tcoupl                  = V-rescale
tc-grps                 = Protein_HEM Water_and_ions
tau-t                   = 0.1 0.1
ref-t                   = 300 300

pcoupl                  = no

pbc                     = xyz

gen-vel                 = yes
gen-temp                = 300
gen-seed                = 20260825
ld-seed                 = 20260825

nstcomm                 = 100
comm-mode               = Linear

nstlog                  = 500
nstenergy               = 500
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0

nstxout-compressed      = 1000
compressed-x-precision  = 1000
compressed-x-grps       = Protein_HEM
MDPFILE


cat > "$MDP/npt.mdp" <<'MDPFILE'
; ============================================================
; NPT avec position restraints - 300 ps
; ============================================================

define                  = -DPOSRES

integrator              = md
dt                      = 0.002
nsteps                  = 150000

continuation            = yes

constraints             = all-bonds
constraint-algorithm    = lincs
lincs-iter              = 1
lincs-order             = 4

cutoff-scheme           = Verlet
nstlist                 = 10

coulombtype             = Reaction-Field
rcoulomb                = 1.4
epsilon-r               = 1
epsilon-rf              = 61

vdwtype                 = Cut-off
rvdw                    = 1.4
DispCorr                = no

tcoupl                  = V-rescale
tc-grps                 = Protein_HEM Water_and_ions
tau-t                   = 0.1 0.1
ref-t                   = 300 300

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau-p                   = 5.0
ref-p                   = 1.0
compressibility         = 4.5e-5

; Important avec position restraints + pression
refcoord-scaling        = all

pbc                     = xyz

gen-vel                 = no
ld-seed                 = 20260825

nstcomm                 = 100
comm-mode               = Linear

nstlog                  = 500
nstenergy               = 500
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0

nstxout-compressed      = 1000
compressed-x-precision  = 1000
compressed-x-grps       = Protein_HEM
MDPFILE


cat > "$MDP/npt_free.mdp" <<'MDPFILE'
; ============================================================
; NPT sans position restraints - 200 ps
; ============================================================

integrator              = md
dt                      = 0.002
nsteps                  = 100000

continuation            = yes

constraints             = all-bonds
constraint-algorithm    = lincs
lincs-iter              = 1
lincs-order             = 4

cutoff-scheme           = Verlet
nstlist                 = 10

coulombtype             = Reaction-Field
rcoulomb                = 1.4
epsilon-r               = 1
epsilon-rf              = 61

vdwtype                 = Cut-off
rvdw                    = 1.4
DispCorr                = no

tcoupl                  = V-rescale
tc-grps                 = Protein_HEM Water_and_ions
tau-t                   = 0.1 0.1
ref-t                   = 300 300

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau-p                   = 5.0
ref-p                   = 1.0
compressibility         = 4.5e-5

pbc                     = xyz

gen-vel                 = no
ld-seed                 = 20260825

nstcomm                 = 100
comm-mode               = Linear

nstlog                  = 500
nstenergy               = 500
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0

nstxout-compressed      = 1000
compressed-x-precision  = 1000
compressed-x-grps       = Protein_HEM
MDPFILE


cat > "$MDP/md_1ns.mdp" <<'MDPFILE'
; ============================================================
; Production pilote - 1 ns
; ============================================================

integrator              = md
dt                      = 0.002

; 500000 x 0.002 ps = 1000 ps = 1 ns
nsteps                  = 500000

continuation            = yes

constraints             = all-bonds
constraint-algorithm    = lincs
lincs-iter              = 1
lincs-order             = 4

cutoff-scheme           = Verlet
nstlist                 = 10

coulombtype             = Reaction-Field
rcoulomb                = 1.4
epsilon-r               = 1
epsilon-rf              = 61

vdwtype                 = Cut-off
rvdw                    = 1.4
DispCorr                = no

tcoupl                  = V-rescale
tc-grps                 = Protein_HEM Water_and_ions
tau-t                   = 0.1 0.1
ref-t                   = 300 300

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau-p                   = 5.0
ref-p                   = 1.0
compressibility         = 4.5e-5

pbc                     = xyz

gen-vel                 = no
ld-seed                 = 20260825

nstcomm                 = 100
comm-mode               = Linear

nstlog                  = 500
nstenergy               = 500

nstxout                  = 0
nstvout                  = 0
nstfout                  = 0

; une image toutes les 2 ps
nstxout-compressed      = 1000
compressed-x-precision  = 1000
compressed-x-grps       = Protein_HEM
MDPFILE


# ============================================================
# 3. VALIDATION GROMPP
# ============================================================

validate_grompp()
{
    local LOG="$1"
    local TPR="$2"

    if grep -Eq \
        '^ERROR|Fatal error|No default G96Angle|No default G96Bond' \
        "$LOG"
    then
        fail "Erreur grompp détectée dans $LOG"
    fi

    local NWARN
    NWARN=$(grep -c '^WARNING [0-9]' "$LOG" || true)

    if [[ "$NWARN" -ne 1 ]]
    then
        fail \
        "$LOG contient $NWARN warning(s), alors qu'un seul warning GROMOS est autorisé."
    fi

    if ! grep -q \
        'The GROMOS force fields have been parametrized' \
        "$LOG"
    then
        fail \
        "Le warning de $LOG n'est pas le warning GROMOS attendu."
    fi

    if [[ ! -s "$TPR" ]]
    then
        fail "TPR absent : $TPR"
    fi
}


# ============================================================
# 4. VALIDATION MDRUN
# ============================================================

validate_mdrun()
{
    local DIR="$1"
    local STAGE="$2"

    if grep -Eqi \
        'LINCS WARNING|Fatal error|Segmentation fault|constraint failure' \
        "$DIR/${STAGE}.log" \
        "$DIR/mdrun_${STAGE}.console.log"
    then
        fail "$STAGE : erreur de dynamique détectée dans $DIR"
    fi

    if ! grep -q 'Finished mdrun' "$DIR/${STAGE}.log"
    then
        fail "$STAGE : mdrun ne s'est pas terminé normalement dans $DIR"
    fi

    [[ -s "$DIR/${STAGE}.gro" ]] \
        || fail "$STAGE : fichier GRO absent"

    [[ -s "$DIR/${STAGE}.edr" ]] \
        || fail "$STAGE : fichier EDR absent"
}


# ============================================================
# 5. GROMPP + MDRUN D'UNE ÉTAPE
# ============================================================

run_stage()
{
    local DIR="$1"
    local LABEL="$2"
    local STAGE="$3"
    local MDPFILE="$4"
    local COORD="$5"
    local REF="$6"
    local CPT="${7:-}"

    log ""
    log "------------------------------------------------------------"
    log "$LABEL : $STAGE"
    log "------------------------------------------------------------"

    cd "$DIR"

    CMD=(
        gmx grompp
        -f "$MDPFILE"
        -c "$COORD"
        -p topol.top
        -n index_nvt.ndx
        -o "${STAGE}.tpr"
        -maxwarn 1
    )

    if [[ -n "$REF" ]]
    then
        CMD+=(-r "$REF")
    fi

    if [[ -n "$CPT" ]]
    then
        CMD+=(-t "$CPT")
    fi

    "${CMD[@]}" \
        2>&1 | tee "grompp_${STAGE}.log"

    validate_grompp \
        "grompp_${STAGE}.log" \
        "${STAGE}.tpr"

    # Permet de reprendre une étape interrompue
    EXTRA=()

    if [[ -f "${STAGE}.cpt" ]] \
       && ! grep -q 'Finished mdrun' "${STAGE}.log" 2>/dev/null
    then
        EXTRA=(
            -cpi "${STAGE}.cpt"
            -append
        )
    fi

    if [[ -f "${STAGE}.log" ]] \
       && grep -q 'Finished mdrun' "${STAGE}.log"
    then
        log "$LABEL $STAGE : déjà terminé, conservation des fichiers."
    else
        gmx mdrun \
            -deffnm "$STAGE" \
            -ntmpi 1 \
            -ntomp "$THREADS" \
            -v \
            "${EXTRA[@]}" \
            2>&1 | tee "mdrun_${STAGE}.console.log"
    fi

    validate_mdrun "$DIR" "$STAGE"

    log "$LABEL $STAGE : OK"
}


# ============================================================
# 6. EXTRACTION THERMODYNAMIQUE
# ============================================================

check_temperature()
{
    local DIR="$1"
    local LABEL="$2"
    local STAGE="$3"

    cd "$DIR"

    printf "Temperature\n0\n" |
        gmx energy \
        -f "${STAGE}.edr" \
        -o "${STAGE}_temperature.xvg" \
        -xvg none \
        > "energy_${STAGE}_temperature.log" 2>&1

    python3 - \
        "${STAGE}_temperature.xvg" \
        "$LABEL" \
        "$STAGE" \
        "$STATUS" <<'PY'
import sys

path, label, stage, status = sys.argv[1:]

data = []

for line in open(path):
    s = line.strip()

    if not s or s.startswith(("#", "@")):
        continue

    f = s.split()

    if len(f) >= 2:
        data.append(float(f[1]))

if len(data) < 5:
    raise SystemExit(
        f"{label} {stage}: données de température insuffisantes"
    )

second_half = data[len(data)//2:]
mean = sum(second_half) / len(second_half)

msg = (
    f"{label} {stage}: "
    f"T moyenne seconde moitié = {mean:.2f} K"
)

print(msg)

with open(status, "a") as f:
    f.write(msg + "\n")

if not (295.0 <= mean <= 305.0):
    raise SystemExit(
        f"{label} {stage}: température hors plage 295-305 K"
    )
PY
}


check_npt()
{
    local DIR="$1"
    local LABEL="$2"
    local STAGE="$3"

    cd "$DIR"

    printf "Temperature\nPressure\nDensity\n0\n" |
        gmx energy \
        -f "${STAGE}.edr" \
        -o "${STAGE}_thermo.xvg" \
        -xvg none \
        > "energy_${STAGE}_thermo.log" 2>&1

    python3 - \
        "${STAGE}_thermo.xvg" \
        "$LABEL" \
        "$STAGE" \
        "$STATUS" <<'PY'
import sys

path, label, stage, status = sys.argv[1:]

rows = []

for line in open(path):
    s = line.strip()

    if not s or s.startswith(("#", "@")):
        continue

    f = s.split()

    if len(f) >= 4:
        rows.append(
            (
                float(f[1]),  # Temperature
                float(f[2]),  # Pressure
                float(f[3]),  # Density
            )
        )

if len(rows) < 10:
    raise SystemExit(
        f"{label} {stage}: données thermodynamiques insuffisantes"
    )

half = rows[len(rows)//2:]

temp = sum(x[0] for x in half) / len(half)
pressure = sum(x[1] for x in half) / len(half)
density = sum(x[2] for x in half) / len(half)

msg = (
    f"{label} {stage}: "
    f"T={temp:.2f} K | "
    f"P={pressure:.2f} bar | "
    f"densité={density:.2f} kg/m3"
)

print(msg)

with open(status, "a") as f:
    f.write(msg + "\n")

if not (295.0 <= temp <= 305.0):
    raise SystemExit(
        f"{label} {stage}: température hors plage"
    )

# Plage volontairement prudente à cause de GROMOS54A7 moderne
if not (950.0 <= density <= 1100.0):
    raise SystemExit(
        f"{label} {stage}: densité physiquement suspecte"
    )
PY
}


# ============================================================
# 7. NVT
# ============================================================

log ""
log "============================================================"
log "NVT 100 ps"
log "============================================================"

run_stage \
    "$WT" "WT" "nvt" \
    "$MDP/nvt.mdp" \
    "em.gro" "em.gro"

run_stage \
    "$HBS" "HbS" "nvt" \
    "$MDP/nvt.mdp" \
    "em.gro" "em.gro"

check_temperature "$WT" "WT" "nvt"
check_temperature "$HBS" "HbS" "nvt"


# ============================================================
# 8. NPT RESTREINT
# ============================================================

log ""
log "============================================================"
log "NPT restreint 300 ps"
log "============================================================"

run_stage \
    "$WT" "WT" "npt" \
    "$MDP/npt.mdp" \
    "nvt.gro" "nvt.gro" "nvt.cpt"

run_stage \
    "$HBS" "HbS" "npt" \
    "$MDP/npt.mdp" \
    "nvt.gro" "nvt.gro" "nvt.cpt"

check_npt "$WT" "WT" "npt"
check_npt "$HBS" "HbS" "npt"


# ============================================================
# 9. NPT SANS CONTRAINTES DE POSITION
# ============================================================

log ""
log "============================================================"
log "NPT libre 200 ps"
log "============================================================"

run_stage \
    "$WT" "WT" "npt_free" \
    "$MDP/npt_free.mdp" \
    "npt.gro" "" "npt.cpt"

run_stage \
    "$HBS" "HbS" "npt_free" \
    "$MDP/npt_free.mdp" \
    "npt.gro" "" "npt.cpt"

check_npt "$WT" "WT" "npt_free"
check_npt "$HBS" "HbS" "npt_free"


# ============================================================
# 10. PRODUCTION PILOTE 1 ns
# ============================================================

log ""
log "============================================================"
log "PRODUCTION PILOTE 1 ns"
log "============================================================"

run_stage \
    "$WT" "WT" "md_1ns" \
    "$MDP/md_1ns.mdp" \
    "npt_free.gro" "" "npt_free.cpt"

run_stage \
    "$HBS" "HbS" "md_1ns" \
    "$MDP/md_1ns.mdp" \
    "npt_free.gro" "" "npt_free.cpt"

check_npt "$WT" "WT" "md_1ns"
check_npt "$HBS" "HbS" "md_1ns"


# ============================================================
# 11. VALIDATION FINALE
# ============================================================

log ""
log "============================================================"
log "PIPELINE FINAL VALIDÉ"
log "============================================================"
log "WT  : NVT + NPT + NPT libre + production 1 ns terminés"
log "HbS : NVT + NPT + NPT libre + production 1 ns terminés"
log ""
log "Aucun LINCS WARNING détecté."
log "Aucune erreur bonded détectée."
log "Aucun nouveau warning grompp accepté."
log ""
log "Étape suivante : RMSD / RMSF / Rg / figures."

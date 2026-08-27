#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

WT="$ROOT/workflow/wt_localff"
HBS="$ROOT/workflow/hbs_localff"

OUT="$ROOT/results/analysis"
FIG="$ROOT/figures"

mkdir -p "$OUT" "$FIG"

STATUS="$OUT/FINAL_ANALYSIS_STATUS.txt"
: > "$STATUS"

log() {
    echo "$*" | tee -a "$STATUS"
}

fail() {
    log ""
    log "ERREUR ANALYSE : $*"
    exit 1
}

# ============================================================
# Analyse d'un système
# ============================================================

analyze_system() {

    local DIR="$1"
    local LABEL="$2"
    local PREFIX="$3"
    local EXPECTED_ATOMS="$4"

    log ""
    log "============================================================"
    log "ANALYSE : $LABEL"
    log "============================================================"

    cd "$DIR"

    # --------------------------------------------------------
    # 1. Vérifier fichiers indispensables
    # --------------------------------------------------------

    for file in \
        md_1ns.tpr \
        md_1ns.xtc \
        md_1ns.edr \
        md_1ns.gro \
        topol.top \
        index_nvt.ndx
    do
        [[ -s "$file" ]] || fail "$LABEL : fichier absent : $file"
    done

    log "$LABEL : fichiers de production présents"

    # --------------------------------------------------------
    # 2. Contrôle d'intégrité de la trajectoire
    # --------------------------------------------------------

    gmx check \
        -f md_1ns.xtc \
        > "$OUT/${PREFIX}_check_xtc.log" 2>&1 \
        || fail "$LABEL : gmx check a échoué"

    log "$LABEL : XTC lisible"

    # --------------------------------------------------------
    # 3. TPR réduit exactement aux atomes Protein_HEM
    #
    # index_nvt.ndx :
    # groupe 0 = Protein_HEM
    # --------------------------------------------------------

    printf "0\n" |
        gmx convert-tpr \
        -s md_1ns.tpr \
        -n index_nvt.ndx \
        -o "$OUT/${PREFIX}_protein_heme.tpr" \
        > "$OUT/${PREFIX}_convert_tpr.log" 2>&1 \
        || fail "$LABEL : convert-tpr a échoué"

    # --------------------------------------------------------
    # 4. Créer indices d'analyse depuis topol.top
    # --------------------------------------------------------

    python3 - \
        topol.top \
        "$OUT/${PREFIX}" \
        "$LABEL" \
        "$EXPECTED_ATOMS" <<'PY'

import sys
from pathlib import Path

top_path = Path(sys.argv[1])
prefix = Path(sys.argv[2])
label = sys.argv[3]
expected_atoms = int(sys.argv[4])

atoms = []

active = False

for line in top_path.read_text().splitlines():

    s = line.strip()

    if s.startswith("["):
        active = (s.lower() == "[ atoms ]")
        continue

    if not active:
        continue

    if not s or s.startswith(";"):
        continue

    f = s.split()

    if len(f) >= 8 and f[0].isdigit():

        atoms.append({
            "nr": int(f[0]),
            "resnr": int(f[2]),
            "resname": f[3],
            "atom": f[4]
        })


if len(atoms) != expected_atoms:
    raise SystemExit(
        f"{label}: {len(atoms)} atomes topology, "
        f"attendu {expected_atoms}"
    )


protein_heme = [
    a["nr"]
    for a in atoms
]

backbone = [
    a["nr"]
    for a in atoms
    if a["resname"] != "HEM"
    and a["atom"] in {"N", "CA", "C"}
]

calpha_records = [
    a
    for a in atoms
    if a["resname"] != "HEM"
    and a["atom"] == "CA"
]

calpha = [
    a["nr"]
    for a in calpha_records
]

heme_atoms = [
    a
    for a in atoms
    if a["resname"] == "HEM"
]


# ------------------------------------------------------------
# Identifier les quatre HEM comme quatre blocs topologiques
# successifs.
#
# On ne compare plus le nombre d'atomes HEM au nombre
# d'atomes lourds du PDB/posre.
# ------------------------------------------------------------

heme_blocks = []
current = []

for a in heme_atoms:

    if (
        current
        and a["resnr"] != current[-1]["resnr"]
    ):
        heme_blocks.append(current)
        current = []

    current.append(a)

if current:
    heme_blocks.append(current)

heme_sizes = [
    len(block)
    for block in heme_blocks
]


# ------------------------------------------------------------
# Contrôles structurels
# ------------------------------------------------------------

if len(calpha) != 574:
    raise SystemExit(
        f"{label}: {len(calpha)} C-alpha "
        f"au lieu de 574"
    )

if len(backbone) != 574 * 3:
    raise SystemExit(
        f"{label}: backbone={len(backbone)} "
        f"au lieu de 1722"
    )

if len(heme_blocks) != 4:
    raise SystemExit(
        f"{label}: {len(heme_blocks)} blocs HEM "
        f"au lieu de 4 ; tailles={heme_sizes}"
    )

if len(set(heme_sizes)) != 1:
    raise SystemExit(
        f"{label}: les quatre HEM n'ont pas "
        f"la même taille topologique : {heme_sizes}"
    )

if sum(heme_sizes) != len(heme_atoms):
    raise SystemExit(
        f"{label}: incohérence interne dans "
        f"le comptage des HEM"
    )


# ------------------------------------------------------------
# Vérification explicite des deux beta6
#
# Ordre des chaînes :
# A = 141 aa
# B = 146 aa
# C = 141 aa
# D = 146 aa
#
# beta6(B) = 141 + 6 = 147
# beta6(D) = 141 + 146 + 141 + 6 = 434
# ------------------------------------------------------------

beta6_positions = (147, 434)

variant_records = [
    calpha_records[pos - 1]
    for pos in beta6_positions
]

expected_variant = (
    "GLU"
    if label == "WT"
    else "VAL"
)

variant_names = [
    a["resname"]
    for a in variant_records
]

if variant_names != [
    expected_variant,
    expected_variant
]:
    raise SystemExit(
        f"{label}: résidus beta6 inattendus : "
        f"{variant_names}, attendu "
        f"{expected_variant}/{expected_variant}"
    )

variant_ca = [
    a["nr"]
    for a in variant_records
]


print(
    f"{label}: "
    f"Protein_HEM={len(protein_heme)}, "
    f"Backbone={len(backbone)}, "
    f"C-alpha={len(calpha)}, "
    f"HEM_sites={len(heme_atoms)}, "
    f"HEM_blocks={heme_sizes}, "
    f"beta6={variant_names}"
)


def write_ndx(path, name, ids):

    with open(path, "w") as f:

        f.write(f"[ {name} ]\n")

        for i in range(0, len(ids), 15):

            f.write(
                " ".join(
                    str(x)
                    for x in ids[i:i+15]
                )
                + "\n"
            )


write_ndx(
    str(prefix) + "_protein_heme.ndx",
    "Protein_HEM",
    protein_heme
)

write_ndx(
    str(prefix) + "_backbone.ndx",
    "Backbone",
    backbone
)

write_ndx(
    str(prefix) + "_calpha.ndx",
    "C_alpha",
    calpha
)

write_ndx(
    str(prefix) + "_variant_ca.ndx",
    "Beta6_CA",
    variant_ca
)


print(
    f"{label}: "
    f"Protein_HEM={len(protein_heme)}, "
    f"Backbone={len(backbone)}, "
    f"C-alpha={len(calpha)}, "
    f"HEM={len(heme_atoms)}, "
    f"Beta6_CA={len(variant_ca)}"
)

PY

    # --------------------------------------------------------
    # 5. Corriger représentation PBC
    # --------------------------------------------------------

    printf "0\n" |
        gmx trjconv \
        -s "$OUT/${PREFIX}_protein_heme.tpr" \
        -f md_1ns.xtc \
        -n "$OUT/${PREFIX}_protein_heme.ndx" \
        -o "$OUT/${PREFIX}_whole.xtc" \
        -pbc nojump \
        > "$OUT/${PREFIX}_trjconv.log" 2>&1 \
        || fail "$LABEL : correction PBC impossible"

    gmx check \
        -f "$OUT/${PREFIX}_whole.xtc" \
        > "$OUT/${PREFIX}_check_whole.log" 2>&1 \
        || fail "$LABEL : trajectoire PBC corrigée invalide"

    log "$LABEL : correction PBC OK"

    # --------------------------------------------------------
    # 6. RMSD backbone
    # --------------------------------------------------------

    printf "0\n0\n" |
        gmx rms \
        -s "$OUT/${PREFIX}_protein_heme.tpr" \
        -f "$OUT/${PREFIX}_whole.xtc" \
        -n "$OUT/${PREFIX}_backbone.ndx" \
        -o "$OUT/${PREFIX}_rmsd_backbone.xvg" \
        -tu ns \
        -xvg none \
        > "$OUT/${PREFIX}_rmsd.log" 2>&1 \
        || fail "$LABEL : RMSD impossible"

    # --------------------------------------------------------
    # 7. RMSF C-alpha
    #
    # On analyse 200-1000 ps afin d'écarter le début
    # de la production pour l'estimation des fluctuations.
    # --------------------------------------------------------

    printf "0\n" |
        gmx rmsf \
        -s "$OUT/${PREFIX}_protein_heme.tpr" \
        -f "$OUT/${PREFIX}_whole.xtc" \
        -n "$OUT/${PREFIX}_calpha.ndx" \
        -o "$OUT/${PREFIX}_rmsf_ca.xvg" \
        -res \
        -b 200 \
        -fit yes \
        -xvg none \
        > "$OUT/${PREFIX}_rmsf.log" 2>&1 \
        || fail "$LABEL : RMSF impossible"

    # --------------------------------------------------------
    # 8. Rayon de giration Protein + HEM
    # --------------------------------------------------------

    printf "0\n" |
        gmx gyrate \
        -s "$OUT/${PREFIX}_protein_heme.tpr" \
        -f "$OUT/${PREFIX}_whole.xtc" \
        -n "$OUT/${PREFIX}_protein_heme.ndx" \
        -o "$OUT/${PREFIX}_rg.xvg" \
        -xvg none \
        > "$OUT/${PREFIX}_rg.log" 2>&1 \
        || fail "$LABEL : rayon de giration impossible"

    # Contrôles des sorties
    for file in \
        "$OUT/${PREFIX}_rmsd_backbone.xvg" \
        "$OUT/${PREFIX}_rmsf_ca.xvg" \
        "$OUT/${PREFIX}_rg.xvg"
    do

        [[ -s "$file" ]] \
            || fail "$LABEL : résultat absent : $file"

    done

    log "$LABEL : RMSD OK"
    log "$LABEL : RMSF OK"
    log "$LABEL : Rg OK"
}


# ============================================================
# WT / HbS
# ============================================================

analyze_system \
    "$WT" \
    "WT" \
    "wt" \
    5804

analyze_system \
    "$HBS" \
    "HbS" \
    "hbs" \
    5800


# ============================================================
# Statistiques + figures
# ============================================================

cd "$ROOT"

python3 - <<'PY'

from pathlib import Path
from statistics import fmean, pstdev
import csv

root = Path.cwd()

out = root / "results" / "analysis"
figdir = root / "figures"


def read_numeric(path, min_cols=2):

    data = []

    for line in Path(path).read_text().splitlines():

        s = line.strip()

        if not s or s.startswith(("#", "@")):
            continue

        f = s.split()

        if len(f) >= min_cols:
            data.append(
                [float(x) for x in f]
            )

    return data


# ------------------------------------------------------------
# Charger RMSD
# x = ns
# ------------------------------------------------------------

wt_rmsd = read_numeric(
    out / "wt_rmsd_backbone.xvg"
)

hbs_rmsd = read_numeric(
    out / "hbs_rmsd_backbone.xvg"
)


# ------------------------------------------------------------
# Charger Rg
# x = ps dans gmx gyrate
# ------------------------------------------------------------

wt_rg = read_numeric(
    out / "wt_rg.xvg"
)

hbs_rg = read_numeric(
    out / "hbs_rg.xvg"
)


# ------------------------------------------------------------
# Charger RMSF
# ------------------------------------------------------------

wt_rmsf = read_numeric(
    out / "wt_rmsf_ca.xvg"
)

hbs_rmsf = read_numeric(
    out / "hbs_rmsf_ca.xvg"
)


if len(wt_rmsf) != 574:
    raise SystemExit(
        f"WT RMSF contient {len(wt_rmsf)} résidus, attendu 574"
    )

if len(hbs_rmsf) != 574:
    raise SystemExit(
        f"HbS RMSF contient {len(hbs_rmsf)} résidus, attendu 574"
    )


# ------------------------------------------------------------
# Thermodynamique produite par finalize_md.sh
# colonnes :
# temps, T, P, densité
# ------------------------------------------------------------

wt_thermo = read_numeric(
    root / "workflow" / "wt_localff"
    / "md_1ns_thermo.xvg",
    4
)

hbs_thermo = read_numeric(
    root / "workflow" / "hbs_localff"
    / "md_1ns_thermo.xvg",
    4
)


def mean_sd(values):

    return (
        fmean(values),
        pstdev(values)
    )


# ------------------------------------------------------------
# Statistiques RMSD :
# considérer t >= 0.2 ns
# ------------------------------------------------------------

def rmsd_stats(data):

    values = [
        row[1]
        for row in data
        if row[0] >= 0.2
    ]

    return {
        "mean": fmean(values),
        "sd": pstdev(values),
        "min": min(values),
        "max": max(values),
        "final": data[-1][1]
    }


# ------------------------------------------------------------
# Rg :
# sortie en ps, analyse >= 200 ps
# ------------------------------------------------------------

def rg_stats(data):

    values = [
        row[1]
        for row in data
        if row[0] >= 200
    ]

    return {
        "mean": fmean(values),
        "sd": pstdev(values),
        "min": min(values),
        "max": max(values),
        "final": data[-1][1]
    }


# ------------------------------------------------------------
# RMSF C-alpha
#
# Topologie :
# A = 141 résidus
# B = 146
# C = 141
# D = 146
#
# β6 chaîne B -> indice séquentiel 147
# β6 chaîne D -> indice séquentiel 434
# ------------------------------------------------------------

BETA6_POSITIONS = [147, 434]


def rmsf_stats(data):

    values = [row[1] for row in data]

    beta6 = [
        values[pos - 1]
        for pos in BETA6_POSITIONS
    ]

    local = []

    for pos in BETA6_POSITIONS:

        start = max(1, pos - 5)
        end = min(len(values), pos + 5)

        local.extend(
            values[start - 1:end]
        )

    return {
        "mean": fmean(values),
        "sd": pstdev(values),
        "max": max(values),
        "beta6_B": beta6[0],
        "beta6_D": beta6[1],
        "beta6_mean": fmean(beta6),
        "local_pm5_mean": fmean(local)
    }


def thermo_stats(data):

    # Derniers 80 % de la production
    start = len(data) // 5

    subset = data[start:]

    temp = [r[1] for r in subset]
    pressure = [r[2] for r in subset]
    density = [r[3] for r in subset]

    tm, ts = mean_sd(temp)
    pm, ps = mean_sd(pressure)
    dm, ds = mean_sd(density)

    return {
        "T_mean": tm,
        "T_sd": ts,
        "P_mean": pm,
        "P_sd": ps,
        "D_mean": dm,
        "D_sd": ds
    }


wr = rmsd_stats(wt_rmsd)
hr = rmsd_stats(hbs_rmsd)

wg = rg_stats(wt_rg)
hg = rg_stats(hbs_rg)

wf = rmsf_stats(wt_rmsf)
hf = rmsf_stats(hbs_rmsf)

wt = thermo_stats(wt_thermo)
ht = thermo_stats(hbs_thermo)


# ------------------------------------------------------------
# CSV final
# ------------------------------------------------------------

rows = [

    (
        "RMSD backbone moyen (0.2-1 ns)",
        wr["mean"],
        hr["mean"],
        "nm"
    ),

    (
        "RMSD backbone SD",
        wr["sd"],
        hr["sd"],
        "nm"
    ),

    (
        "RMSD final",
        wr["final"],
        hr["final"],
        "nm"
    ),

    (
        "RMSF C-alpha moyen",
        wf["mean"],
        hf["mean"],
        "nm"
    ),

    (
        "RMSF beta6 chaine B",
        wf["beta6_B"],
        hf["beta6_B"],
        "nm"
    ),

    (
        "RMSF beta6 chaine D",
        wf["beta6_D"],
        hf["beta6_D"],
        "nm"
    ),

    (
        "RMSF beta6 moyen",
        wf["beta6_mean"],
        hf["beta6_mean"],
        "nm"
    ),

    (
        "RMSF voisinage beta6 +/-5",
        wf["local_pm5_mean"],
        hf["local_pm5_mean"],
        "nm"
    ),

    (
        "Rayon de giration moyen",
        wg["mean"],
        hg["mean"],
        "nm"
    ),

    (
        "Rayon de giration SD",
        wg["sd"],
        hg["sd"],
        "nm"
    ),

    (
        "Temperature moyenne",
        wt["T_mean"],
        ht["T_mean"],
        "K"
    ),

    (
        "Pression moyenne",
        wt["P_mean"],
        ht["P_mean"],
        "bar"
    ),

    (
        "Pression SD",
        wt["P_sd"],
        ht["P_sd"],
        "bar"
    ),

    (
        "Densite moyenne",
        wt["D_mean"],
        ht["D_mean"],
        "kg/m3"
    ),
]


csv_path = out / "final_metrics.csv"

with csv_path.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(f)

    writer.writerow([
        "Metrique",
        "WT",
        "HbS",
        "Delta_HbS_WT",
        "Unite"
    ])

    for name, w, h, unit in rows:

        writer.writerow([
            name,
            f"{w:.6f}",
            f"{h:.6f}",
            f"{h-w:.6f}",
            unit
        ])


# ------------------------------------------------------------
# Résumé lisible
# ------------------------------------------------------------

summary = f"""
============================================================
RESULTATS FINAUX DE DYNAMIQUE MOLECULAIRE
============================================================

Production :
    WT  = 1.0 ns
    HbS = 1.0 ns

------------------------------------------------------------
THERMODYNAMIQUE
------------------------------------------------------------

WT
    Temperature = {wt['T_mean']:.2f} +/- {wt['T_sd']:.2f} K
    Pression    = {wt['P_mean']:.2f} +/- {wt['P_sd']:.2f} bar
    Densite     = {wt['D_mean']:.2f} +/- {wt['D_sd']:.2f} kg/m3

HbS
    Temperature = {ht['T_mean']:.2f} +/- {ht['T_sd']:.2f} K
    Pression    = {ht['P_mean']:.2f} +/- {ht['P_sd']:.2f} bar
    Densite     = {ht['D_mean']:.2f} +/- {ht['D_sd']:.2f} kg/m3

------------------------------------------------------------
RMSD BACKBONE
Moyenne calculee sur 0.2-1.0 ns
------------------------------------------------------------

WT
    RMSD moyen = {wr['mean']:.4f} +/- {wr['sd']:.4f} nm
    RMSD final = {wr['final']:.4f} nm

HbS
    RMSD moyen = {hr['mean']:.4f} +/- {hr['sd']:.4f} nm
    RMSD final = {hr['final']:.4f} nm

Delta moyen HbS-WT = {hr['mean']-wr['mean']:+.4f} nm

------------------------------------------------------------
RMSF C-ALPHA
Calcul sur 200-1000 ps
------------------------------------------------------------

WT
    RMSF moyen              = {wf['mean']:.4f} nm
    Beta6 chaine B          = {wf['beta6_B']:.4f} nm
    Beta6 chaine D          = {wf['beta6_D']:.4f} nm
    Beta6 moyen             = {wf['beta6_mean']:.4f} nm
    Voisinage beta6 +/-5    = {wf['local_pm5_mean']:.4f} nm

HbS
    RMSF moyen              = {hf['mean']:.4f} nm
    Beta6 chaine B          = {hf['beta6_B']:.4f} nm
    Beta6 chaine D          = {hf['beta6_D']:.4f} nm
    Beta6 moyen             = {hf['beta6_mean']:.4f} nm
    Voisinage beta6 +/-5    = {hf['local_pm5_mean']:.4f} nm

Delta beta6 HbS-WT = {hf['beta6_mean']-wf['beta6_mean']:+.4f} nm

------------------------------------------------------------
RAYON DE GIRATION
Moyenne calculee sur 200-1000 ps
------------------------------------------------------------

WT
    Rg = {wg['mean']:.4f} +/- {wg['sd']:.4f} nm

HbS
    Rg = {hg['mean']:.4f} +/- {hg['sd']:.4f} nm

Delta Rg HbS-WT = {hg['mean']-wg['mean']:+.4f} nm

============================================================
IMPORTANT
============================================================

Ces valeurs decrivent une simulation comparative pilote
de 1 ns.

Elles ne constituent pas une estimation de DeltaG et ne
demontrent pas, seules, une destabilisation thermodynamique.

WT et HbS proviennent de deux structures experimentales
differentes ; toute difference doit donc etre interpretee
comme comparative/exploratoire.

============================================================
"""

summary_path = (
    out / "FINAL_ANALYSIS_SUMMARY.txt"
)

summary_path.write_text(
    summary,
    encoding="utf-8"
)

print(summary)


# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

try:

    import matplotlib.pyplot as plt

except ImportError:

    print(
        "\nATTENTION : matplotlib absent."
        "\nLes calculs sont termines, mais les figures "
        "n'ont pas ete generees."
    )

    raise SystemExit(0)


# RMSD
plt.figure(figsize=(8, 5))

plt.plot(
    [r[0] for r in wt_rmsd],
    [r[1] for r in wt_rmsd],
    label="WT"
)

plt.plot(
    [r[0] for r in hbs_rmsd],
    [r[1] for r in hbs_rmsd],
    label="HbS"
)

plt.xlabel("Temps (ns)")
plt.ylabel("RMSD du squelette (nm)")
plt.title("Comparaison du RMSD : WT vs HbS")
plt.legend()
plt.tight_layout()

plt.savefig(
    figdir / "rmsd_WT_vs_HbS.png",
    dpi=300
)

plt.savefig(
    figdir / "rmsd_WT_vs_HbS.pdf"
)

plt.close()


# Rg
plt.figure(figsize=(8, 5))

plt.plot(
    [r[0] / 1000 for r in wt_rg],
    [r[1] for r in wt_rg],
    label="WT"
)

plt.plot(
    [r[0] / 1000 for r in hbs_rg],
    [r[1] for r in hbs_rg],
    label="HbS"
)

plt.xlabel("Temps (ns)")
plt.ylabel("Rayon de giration (nm)")
plt.title("Compacité structurale : WT vs HbS")
plt.legend()
plt.tight_layout()

plt.savefig(
    figdir / "rg_WT_vs_HbS.png",
    dpi=300
)

plt.savefig(
    figdir / "rg_WT_vs_HbS.pdf"
)

plt.close()


# RMSF
positions = list(
    range(1, 575)
)

plt.figure(figsize=(10, 5))

plt.plot(
    positions,
    [r[1] for r in wt_rmsf],
    label="WT"
)

plt.plot(
    positions,
    [r[1] for r in hbs_rmsf],
    label="HbS"
)

# Limites des chaînes A/B/C/D
for boundary in [
    141.5,
    287.5,
    428.5
]:
    plt.axvline(
        boundary,
        linestyle="--",
        linewidth=0.7,
        alpha=0.5
    )

# Deux positions beta6
for position in BETA6_POSITIONS:
    plt.axvline(
        position,
        linestyle=":",
        linewidth=1.0
    )

plt.xlabel(
    "Position séquentielle dans le tétramère"
)

plt.ylabel(
    "RMSF Cα (nm)"
)

plt.title(
    "Flexibilité résiduelle : WT vs HbS"
)

plt.legend()
plt.tight_layout()

plt.savefig(
    figdir / "rmsf_WT_vs_HbS.png",
    dpi=300
)

plt.savefig(
    figdir / "rmsf_WT_vs_HbS.pdf"
)

plt.close()


print("\nFIGURES GENEREES :")

for name in [
    "rmsd_WT_vs_HbS.png",
    "rmsf_WT_vs_HbS.png",
    "rg_WT_vs_HbS.png"
]:

    print(" -", figdir / name)

PY


# ============================================================
# Validation finale
# ============================================================

log ""
log "============================================================"
log "ANALYSE FINALE TERMINEE"
log "============================================================"

log "RMSD : OK"
log "RMSF : OK"
log "Rg   : OK"

log ""
log "Résumé :"
log "$OUT/FINAL_ANALYSIS_SUMMARY.txt"

log "CSV :"
log "$OUT/final_metrics.csv"

log "Figures :"
log "$FIG/rmsd_WT_vs_HbS.png"
log "$FIG/rmsf_WT_vs_HbS.png"
log "$FIG/rg_WT_vs_HbS.png"


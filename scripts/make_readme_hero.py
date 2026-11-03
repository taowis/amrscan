#!/usr/bin/env python3
"""Build the animated README banner for AMRScan.

    python3 scripts/make_readme_hero.py docs/images/amrscan-hero.svg [teal|blue|gram]

Four scenes on one loop. The first shows the workflow: a scan line passes along
the chromosome of an E. coli cell, three determinants found in the study cohort
light up (blaCTX-M-15 in 889 genomes, gyrA S83L in 2,298, tet(A) in 2,509), each
adds a row to a report with the drug the study maps it to, and the run's
provenance is recorded. The other three follow the E. coli domain-shift study,
with numbers from docs/ecoli-domain-shift-error-decomposition.md and
docs/ecoli-representation-sensitivity.md:
  - random split versus holdout of a whole population;
  - sampled k-mers (1 in 200) cover gyrA S83L in 0 of 2,298 calls, 2,295 at 1 in 50;
  - nalidixic acid holdout AUROC with k-mers only and with known mutations added,
    against the random-split AUROC on the same test isolates.
"""
import math
import random
import re
import sys

W, H = 1280, 480
SCENES = [10.0, 7.0, 7.0, 7.0]           # seconds per scene
STARTS = [sum(SCENES[:i]) for i in range(len(SCENES))]
T = sum(SCENES)

KEYS = ["#0c1a22", "#10262e", "#14363b", "#2a9d8f", "#1d4a4f", "#0f2a30", "#e9fffb", "#7fd1c4",
        "#9fe3d6", "#5fb8ab", "#3f7f86", "#ff8a65", "#2a1d1f", "#ffd2c2", "#0b1820", "#2e4f55",
        "#5f8a84", "#d9efe9", "#8fb3ad", "#0b161c", "#24434a", "#eef6f3", "#86d6c8", "#e6f4f1", "#7f9a96"]
THEMES = {
    "teal": KEYS,
    "blue": ["#0a1530", "#0e1d42", "#132a5c", "#3d7bff", "#1e3a7a", "#10224d", "#eaf2ff", "#7fb2ff",
             "#a9cbff", "#5e93f0", "#3f5f9e", "#ff8a65", "#2b1a1d", "#ffd2c2", "#081128", "#2a3f73",
             "#6f86b8", "#e3ecff", "#9fb2d9", "#081128", "#22386a", "#f2f6ff", "#8fbfff", "#e3ecff", "#8394bd"],
    "gram": ["#0b0910", "#120d19", "#1a1222", "#c2478a", "#4a1d3a", "#2a1024", "#fff0f7", "#f08cc0",
             "#ffb3d9", "#d86aa6", "#8a4a72", "#ffc94d", "#2a2312", "#ffe9b0", "#0e0b12", "#3a2a40",
             "#8e6f8a", "#f3e6ef", "#bfa2b7", "#0e0b12", "#33263a", "#f7f0f5", "#f39ac8", "#f3e6ef", "#9a8a96"],
}
ACC, ACC2, ACC3, ACC4, HOT = "#7fd1c4", "#5fb8ab", "#3f7f86", "#9fe3d6", "#ff8a65"
TXT, MUTED, FAINT, PANEL, LINE = "#d9efe9", "#8fb3ad", "#5f8a84", "#0b1820", "#2e4f55"
CHIPBG, CHIPTX, BEAM = "#2a1d1f", "#ffd2c2", "#e9fffb"

css = []


def P(scene, sec):
    """percentage of the loop at `sec` seconds into `scene`"""
    return min(100.0, max(0.0, (STARTS[scene] + sec) / T * 100))


def window(name, scene, a, b, fade=0.3):
    pa, pb, f = P(scene, a), P(scene, b), fade / T * 100
    css.append(f".{name}{{opacity:0;animation:{name} {T}s linear infinite}}"
               f"@keyframes {name}{{0%,{pa:.2f}%{{opacity:0}}{pa + f:.2f}%,{pb - f:.2f}%{{opacity:1}}{pb:.2f}%,100%{{opacity:0}}}}")


def caption(step, text):
    return (f'<text x="724" y="62" class="sans" font-size="13" fill="{HOT}" letter-spacing="2" font-weight="700">{step}</text>'
            f'<text x="724" y="88" class="sans" font-size="21" font-weight="600" fill="{TXT}">{text}</text>')


def helix(x0, x1, hy, amp, per, flow=False):
    def strand(phase):
        pts, x = [], x0
        while x <= x1 + 0.1:
            pts.append(f"{x:.1f},{hy + amp * math.sin(2 * math.pi * (x - x0) / per + phase):.1f}")
            x += 3
        return "M " + " L ".join(pts)
    out, x = [], x0 + 8
    while x < x1 - 4:
        a = hy + amp * math.sin(2 * math.pi * (x - x0) / per)
        b = hy + amp * math.sin(2 * math.pi * (x - x0) / per + math.pi)
        depth = abs(math.cos(2 * math.pi * (x - x0) / per))
        out.append(f'<line x1="{x:.1f}" y1="{a:.1f}" x2="{x:.1f}" y2="{b:.1f}" stroke="{ACC}" stroke-opacity="{0.18 + 0.3 * (1 - depth):.2f}" stroke-width="2.2"/>')
        x += 12
    out.append(f'<path d="{strand(0)}" fill="none" stroke="{ACC4}" stroke-width="3" stroke-linecap="round"/>')
    dash = ' class="flow" stroke-dasharray="20 12"' if flow else ""
    out.append(f'<path{dash} d="{strand(math.pi)}" fill="none" stroke="{ACC2}" stroke-width="3" stroke-linecap="round"/>')
    return "".join(out)


def gene(name):
    for g in ("bla", "gyrA", "tet"):
        if name.startswith(g):
            return f'<tspan font-style="italic">{g}</tspan>' + name[len(g):]
    return name


# ---- scene 1: call the determinants ----------------------------------------
window("s0", 0, 0, SCENES[0], fade=0.4)
CX0, CX1, CY0, CY1 = 732, 1220, 116, 256
HY = (CY0 + CY1) / 2
HX0, HX1 = CX0 + 40, CX1 - 40
SCAN_A, SCAN_B = 0.6, 6.0
LOCI = [(842, "blaCTX-M-15", "ceftriaxone", -1), (976, "gyrA S83L", "nalidixic acid", +1),
        (1110, "tet(A)", "tetracycline", -1)]


def t_at(x):
    return SCAN_A + (x - HX0) / (HX1 - HX0) * (SCAN_B - SCAN_A)


css.append(f".scan{{animation:scan {T}s linear infinite}}"
           f"@keyframes scan{{0%,{P(0, SCAN_A):.2f}%{{transform:translateX(0);opacity:0}}{P(0, SCAN_A + .3):.2f}%{{opacity:1}}"
           f"{P(0, SCAN_B - .3):.2f}%{{opacity:1}}{P(0, SCAN_B):.2f}%{{transform:translateX({HX1 - HX0}px);opacity:0}}100%{{transform:translateX({HX1 - HX0}px);opacity:0}}}}")
loci, rows = [], []
END0 = SCENES[0] - 0.4
RY0 = 318
for i, (x, name, drug, side) in enumerate(LOCI):
    t = t_at(x)
    p0, p1, p2, pe, pz = P(0, t - .05), P(0, t + .15), P(0, t + .5), P(0, END0), P(0, SCENES[0])
    css.append(f".hit{i}{{animation:hit{i} {T}s ease-out infinite;transform-box:fill-box;transform-origin:center}}"
               f"@keyframes hit{i}{{0%,{p0:.2f}%{{opacity:0;transform:scale(.6)}}{p1:.2f}%{{opacity:1;transform:scale(1.15)}}{p2:.2f}%{{opacity:1;transform:scale(1)}}{pe:.2f}%{{opacity:1}}{pz:.2f}%,100%{{opacity:0}}}}"
               f".ring{i}{{animation:ring{i} {T}s ease-out infinite;transform-box:fill-box;transform-origin:center}}"
               f"@keyframes ring{i}{{0%,{p0:.2f}%{{opacity:0;transform:scale(.3)}}{P(0, t + .1):.2f}%{{opacity:.85}}{P(0, t + 1.1):.2f}%,100%{{opacity:0;transform:scale(2.4)}}}}"
               f".row{i}{{animation:row{i} {T}s cubic-bezier(.2,.7,.2,1) infinite}}"
               f"@keyframes row{i}{{0%,{P(0, t + .15):.2f}%{{opacity:0;transform:translateX(-12px)}}{P(0, t + .6):.2f}%{{opacity:1;transform:none}}{pe:.2f}%{{opacity:1}}{pz:.2f}%,100%{{opacity:0}}}}")
    cy = HY + side * 72
    w = 14 + len(name) * 9.6
    loci.append(f"""<circle class="ring{i}" cx="{x}" cy="{HY}" r="20" fill="none" stroke="{HOT}" stroke-width="2.5"/>
    <g class="hit{i}">
      <rect x="{x - 22}" y="{HY - 30}" width="44" height="60" rx="10" fill="{HOT}" fill-opacity=".18" stroke="{HOT}" stroke-width="2"/>
      <rect x="{x - w / 2:.1f}" y="{cy - 14}" width="{w:.1f}" height="28" rx="14" fill="{CHIPBG}" stroke="{HOT}" stroke-width="1.5"/>
      <text x="{x}" y="{cy + 5}" text-anchor="middle" class="mono" font-size="15" fill="{CHIPTX}">{gene(name)}</text>
    </g>""")
    y = RY0 + i * 27
    rows.append(f"""<g class="row{i}">
      <text x="760" y="{y}" class="mono" font-size="15" fill="{TXT}">{gene(name)}</text>
      <text x="920" y="{y}" class="mono" font-size="15" fill="{MUTED}">{drug}</text>
      <rect x="1166" y="{y - 16}" width="30" height="22" rx="6" fill="{HOT}"/>
      <text x="1181" y="{y}" text-anchor="middle" class="mono" font-size="15" font-weight="700" fill="{CHIPBG}">R</text>
    </g>""")
window("stamp", 0, SCAN_B + .6, SCENES[0])
s0 = f"""<g class="s0">
  {caption("1 · CALL RESISTANCE DETERMINANTS", "Every call recorded with how it was made")}
  <path class="flag" d="M {CX1 - 2} {HY - 14} C {CX1 + 8} {HY - 30}, {CX1 + 16} {HY + 2}, {CX1 + 24} {HY - 12}" fill="none" stroke="{ACC3}" stroke-width="2" stroke-linecap="round"/>
  <path class="flag" d="M {CX0 + 2} {HY + 10} C {CX0 - 8} {HY - 6}, {CX0 - 14} {HY + 26}, {CX0 - 22} {HY + 12}" fill="none" stroke="{ACC3}" stroke-width="2" stroke-linecap="round" style="animation-delay:-.6s"/>
  <rect x="{CX0}" y="{CY0}" width="{CX1 - CX0}" height="{CY1 - CY0}" rx="{(CY1 - CY0) / 2}" fill="url(#cellfill)" stroke="{ACC}" stroke-opacity=".5" stroke-width="2"/>
  <g clip-path="url(#cellclip)">{helix(HX0, HX1, HY, 24, 92, flow=True)}</g>
  {''.join(loci)}
  <g class="scan"><rect x="{HX0 - 28}" y="{CY0 + 6}" width="28" height="{CY1 - CY0 - 12}" fill="url(#beam)"/>
    <line x1="{HX0}" y1="{CY0 + 6}" x2="{HX0}" y2="{CY1 - 6}" stroke="{BEAM}" stroke-width="2" stroke-opacity=".9"/></g>
  <text x="760" y="{RY0 - 24}" class="sans" font-size="12" fill="{FAINT}" letter-spacing="1.5">DETERMINANT</text>
  <text x="920" y="{RY0 - 24}" class="sans" font-size="12" fill="{FAINT}" letter-spacing="1.5">ANTIMICROBIAL</text>
  <text x="1160" y="{RY0 - 24}" class="sans" font-size="12" fill="{FAINT}" letter-spacing="1.5">CALL</text>
  {''.join(rows)}
  <g class="stamp"><line x1="760" y1="{RY0 + 72}" x2="1196" y2="{RY0 + 72}" stroke="{LINE}"/>
    <path d="M 762 {RY0 + 92} l 6 6 l 11 -13" fill="none" stroke="{ACC}" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>
    <text x="788" y="{RY0 + 97}" class="mono" font-size="13" fill="{ACC}">tool · database · input sha256 recorded</text></g>
</g>"""

# ---- scene 2: random split versus holdout ----------------------------------
window("s1", 1, 0, SCENES[1], fade=0.4)
COLS, ROWS = 14, 6
GX, GY, DX, DY = 746, 140, 33, 33
groups = [ACC, ACC2, ACC4, ACC3]
rng = random.Random(7)
cells = []
for r in range(ROWS):
    for c in range(COLS):
        g = 0 if c < 4 else 1 if c < 7 else 2 if c < 11 else 3
        if rng.random() < 0.12:
            g = min(3, max(0, g + rng.choice([-1, 1])))
        cells.append((r, c, g))
rand_test = set(rng.sample(range(len(cells)), 18))
window("rr", 1, 0.5, 3.4)
window("rg", 1, 3.7, 6.7)
dots = []
for k, (r, c, g) in enumerate(cells):
    x, y = GX + c * DX, GY + r * DY
    dots.append(f'<circle cx="{x}" cy="{y}" r="8" fill="{groups[g]}" fill-opacity=".85"/>')
    if k in rand_test:
        dots.append(f'<circle class="rr" cx="{x}" cy="{y}" r="13" fill="none" stroke="{HOT}" stroke-width="2.2"/>')
    if g == 2:
        dots.append(f'<circle class="rg" cx="{x}" cy="{y}" r="13" fill="none" stroke="{HOT}" stroke-width="2.2"/>')
s1 = f"""<g class="s1">
  {caption("2 · HOLD OUT A POPULATION", "Test on genomes unlike the training data")}
  {''.join(dots)}
  <g class="rr"><circle cx="740" cy="352" r="7" fill="none" stroke="{HOT}" stroke-width="2.2"/>
    <text x="756" y="357" class="sans" font-size="16" fill="{MUTED}">Random split: test isolates from every population</text></g>
  <g class="rg"><circle cx="740" cy="352" r="7" fill="none" stroke="{HOT}" stroke-width="2.2"/>
    <text x="756" y="357" class="sans" font-size="16" fill="{MUTED}">Holdout: a whole lineage, country, project or period</text></g>
  <line x1="724" y1="378" x2="1226" y2="378" stroke="{LINE}"/>
  <text x="724" y="404" class="sans" font-size="15" fill="{TXT}">For most of the 24 drugs a random split overstated AUROC</text>
  <text x="724" y="426" class="sans" font-size="15" fill="{TXT}">only slightly (median gap 0.005 to 0.025); a few drugs lost more.</text>
</g>"""

# ---- scene 3: sampled k-mers miss gyrA S83L -------------------------------
window("s2", 2, 0, SCENES[2], fade=0.4)
KX0, KX1, KY = 744, 1204, 205
MUT_X = 985
sparse = [758, 838, 912, 1062, 1140]
dense = [x for x in range(752, 1196, 21) if x not in sparse]
window("ks", 2, 0.3, 7)
window("c200", 2, 0.3, 3.4)
window("miss", 2, 0.7, 3.4)
window("kd", 2, 3.6, 7)
window("hitm", 2, 3.9, 7)
ticks = "".join(f'<rect x="{x - 9}" y="{KY + 42}" width="18" height="8" rx="3" fill="{ACC}"/>' for x in sparse)
ticks_d = "".join(f'<rect x="{x - 9}" y="{KY + 42}" width="18" height="8" rx="3" fill="{ACC}" fill-opacity=".75"/>' for x in dense)
s2 = f"""<g class="s2">
  {caption("3 · LOOK AT THE FEATURES", "Sampled k-mers can miss a point mutation")}
  {helix(KX0, KX1, KY, 24, 92)}
  <line x1="{MUT_X}" y1="{KY - 40}" x2="{MUT_X}" y2="{KY + 34}" stroke="{HOT}" stroke-width="3"/>
  <rect x="{MUT_X - 54}" y="{KY - 74}" width="108" height="28" rx="14" fill="{CHIPBG}" stroke="{HOT}" stroke-width="1.5"/>
  <text x="{MUT_X}" y="{KY - 55}" text-anchor="middle" class="mono" font-size="15" fill="{CHIPTX}">{gene("gyrA S83L")}</text>
  <g class="ks">{ticks}</g>
  <g class="kd">{ticks_d}<rect x="{MUT_X - 10}" y="{KY + 40}" width="20" height="12" rx="4" fill="{HOT}"/></g>
  <g class="c200"><text x="744" y="{KY + 76}" class="sans" font-size="14" fill="{FAINT}">k-mers kept: 1 in 200</text></g>
  <g class="kd"><text x="744" y="{KY + 76}" class="sans" font-size="14" fill="{FAINT}">k-mers kept: 1 in 50</text></g>
  <g class="miss"><text x="744" y="{KY + 130}" class="sans" font-size="30" font-weight="700" fill="{HOT}">0 of 2,298</text>
    <text x="916" y="{KY + 130}" class="sans" font-size="16" fill="{MUTED}">S83L calls covered by a sampled k-mer</text></g>
  <g class="hitm"><text x="744" y="{KY + 130}" class="sans" font-size="30" font-weight="700" fill="{ACC}">2,295 of 2,298</text>
    <text x="972" y="{KY + 130}" class="sans" font-size="16" fill="{MUTED}">covered at four times the density</text></g>
</g>"""

# ---- scene 4: nalidixic acid holdout AUROC ---------------------------------
window("s3", 3, 0, SCENES[3], fade=0.4)
# holdout, k-mers only (R0), + 62 known mutations (R2), random split on the same test isolates (R0 + R0 gap)
ROWS3 = [("lineage", 0.776, 0.968, 0.880), ("country", 0.765, 0.920, 0.875),
         ("project", 0.698, 0.915, 0.875), ("geography", 0.693, 0.911, 0.868),
         ("time", 0.714, 0.946, 0.845)]
AX0, AX1, V0, V1 = 850, 1180, 0.65, 1.0


def ax(v):
    return AX0 + (v - V0) / (V1 - V0) * (AX1 - AX0)


bars = []
for i, (h, r0, r2, rnd) in enumerate(ROWS3):
    y = 150 + i * 44
    d = 0.8 + i * 0.3
    css.append(f".g{i}{{animation:g{i} {T}s cubic-bezier(.3,.7,.2,1) infinite;transform-box:fill-box;transform-origin:0% 50%}}"
               f"@keyframes g{i}{{0%,{P(3, d):.2f}%{{transform:scaleX(0)}}{P(3, d + 1.1):.2f}%,100%{{transform:scaleX(1)}}}}")
    window(f"p{i}", 3, d + .9, SCENES[3], fade=.25)
    bars.append(f"""<text x="{AX0 - 16}" y="{y + 5}" text-anchor="end" class="sans" font-size="15" fill="{MUTED}">{h}</text>
    <line x1="{AX0}" y1="{y}" x2="{AX1}" y2="{y}" stroke="{LINE}" stroke-width="1"/>
    <line x1="{ax(rnd):.1f}" y1="{y - 12}" x2="{ax(rnd):.1f}" y2="{y + 12}" stroke="{TXT}" stroke-width="2" stroke-dasharray="3 3"/>
    <rect class="g{i}" x="{ax(r0):.1f}" y="{y - 3}" width="{ax(r2) - ax(r0):.1f}" height="6" rx="3" fill="{HOT}" fill-opacity=".55"/>
    <circle cx="{ax(r0):.1f}" cy="{y}" r="7" fill="{PANEL}" stroke="{ACC}" stroke-width="2.4"/>
    <g class="p{i}"><circle cx="{ax(r2):.1f}" cy="{y}" r="7.5" fill="{HOT}"/>
      <text x="{ax(r2) + 14:.1f}" y="{y + 5}" class="mono" font-size="13" fill="{CHIPTX}">+{r2 - r0:.3f}</text></g>""")
axis = "".join(f'<text x="{ax(v):.1f}" y="382" text-anchor="middle" class="mono" font-size="12" fill="{FAINT}">{v:.1f}</text>'
               for v in (0.7, 0.8, 0.9, 1.0))
s3 = f"""<g class="s3">
  {caption("4 · ADD THE KNOWN MUTATIONS", "Nalidixic acid AUROC under each holdout")}
  {''.join(bars)}
  {axis}
  <circle cx="740" cy="418" r="6" fill="{PANEL}" stroke="{ACC}" stroke-width="2.2"/><text x="752" y="423" class="sans" font-size="14" fill="{MUTED}">k-mers only</text>
  <circle cx="860" cy="418" r="6.5" fill="{HOT}"/><text x="872" y="423" class="sans" font-size="14" fill="{MUTED}">+ known mutations</text>
  <line x1="1022" y1="410" x2="1022" y2="426" stroke="{TXT}" stroke-width="2" stroke-dasharray="3 3"/><text x="1032" y="423" class="sans" font-size="14" fill="{MUTED}">random split, k-mers only</text>
</g>"""

# ---- left side, frame, terminal ------------------------------------------
COMMANDS = ["nextflow run main.nf -profile test,docker",
            "scripts/fetch_amrfinder_db.sh --engine docker",
            "nextflow run main.nf -profile gadi",
            "nextflow run main.nf --help"]
CW, TX, TY, CT = 10.8, 132, 336, 4.5
n = len(COMMANDS)
cyc, share = n * CT, 100 / n
term = []
for i, c in enumerate(COMMANDS):
    k, tl = len(c), len(c) * CW
    tp = k * 0.055 / cyc * 100
    css.append(f".c{i}{{opacity:0;animation:show{i} {cyc}s linear infinite}}"
               f"@keyframes show{i}{{0%,{i * share:.2f}%{{opacity:0}}{i * share + 0.01:.2f}%,{(i + 1) * share - 0.5:.2f}%{{opacity:1}}{(i + 1) * share:.2f}%,100%{{opacity:0}}}}"
               f".v{i}{{animation:type{i} {cyc}s infinite}}"
               f"@keyframes type{i}{{0%,{i * share:.2f}%{{transform:translateX(0);animation-timing-function:steps({k},end)}}{i * share + tp:.2f}%,100%{{transform:translateX({tl:.1f}px)}}}}")
    term.append(f'<g class="c{i}"><text x="{TX}" y="{TY}" class="mono" font-size="18" fill="#e6f4f1" textLength="{tl:.1f}" lengthAdjust="spacingAndGlyphs">{c}</text>'
                f'<g class="v{i}"><rect x="{TX - 1}" y="{TY - 20}" width="{tl + 24:.1f}" height="28" fill="#0b161c"/>'
                f'<rect class="caret" x="{TX}" y="{TY - 17}" width="10" height="22" rx="2" fill="{ACC}"/></g></g>')

css.append(f"""
.flow{{animation:flow 2.4s linear infinite}}@keyframes flow{{to{{stroke-dashoffset:-32}}}}
.flag{{animation:wave 1.4s ease-in-out infinite alternate;transform-box:fill-box;transform-origin:0% 50%}}@keyframes wave{{from{{transform:rotate(-5deg)}}to{{transform:rotate(5deg)}}}}
.caret{{animation:blink 1s steps(1) infinite}}@keyframes blink{{50%{{opacity:0}}}}
.rise{{animation:rise 1.2s cubic-bezier(.2,.7,.2,1) both}}.d1{{animation-delay:.15s}}.d2{{animation-delay:.3s}}
@keyframes rise{{from{{opacity:0;transform:translateY(14px)}}to{{opacity:1;transform:none}}}}
.prog{{animation:prog {T}s linear infinite;transform-box:fill-box;transform-origin:0% 50%}}@keyframes prog{{from{{transform:scaleX(0)}}to{{transform:scaleX(1)}}}}
@media (prefers-reduced-motion: reduce){{*{{animation:none!important}}.s1,.s2,.s3,.scan,.prog{{opacity:0}}
  .s0,.hit0,.hit1,.hit2,.row0,.row1,.row2,.stamp,.c0{{opacity:1}}.v0{{transform:translateX(442.8px)}}}}""")

svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" role="img" aria-labelledby="t d">
  <title id="t">AMRScan</title>
  <desc id="d">Four scenes. 1: a scan line passes along the chromosome of an E. coli cell; blaCTX-M-15, gyrA S83L and tet(A) light up and are added to a report as resistant calls for ceftriaxone, nalidixic acid and tetracycline, with the tool, database and input hashes recorded. 2: a random split draws test isolates from every population, a holdout keeps a whole lineage, country, project or period out; for most of the 24 drugs a random split overstated AUROC only slightly. 3: k-mers sampled at 1 in 200 cover gyrA S83L in 0 of 2,298 calls, and 2,295 of 2,298 at 1 in 50. 4: adding known mutations raises nalidixic acid holdout AUROC by 0.155 to 0.232.</desc>
  <style>
.sans{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}}
.mono{{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,"Liberation Mono",monospace}}
{chr(10).join(css)}
  </style>
  <defs>
    <linearGradient id="bg" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="#0c1a22"/><stop offset=".55" stop-color="#10262e"/><stop offset="1" stop-color="#14363b"/></linearGradient>
    <radialGradient id="glow" cx="76%" cy="45%" r="45%"><stop offset="0" stop-color="#2a9d8f" stop-opacity=".2"/><stop offset="1" stop-color="#2a9d8f" stop-opacity="0"/></radialGradient>
    <linearGradient id="cellfill" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#1d4a4f" stop-opacity=".55"/><stop offset="1" stop-color="#0f2a30" stop-opacity=".55"/></linearGradient>
    <linearGradient id="beam" x1="0" y1="0" x2="1" y2="0"><stop offset="0" stop-color="{BEAM}" stop-opacity="0"/><stop offset="1" stop-color="{BEAM}" stop-opacity=".28"/></linearGradient>
    <clipPath id="cellclip"><rect x="{CX0}" y="{CY0}" width="{CX1 - CX0}" height="{CY1 - CY0}" rx="{(CY1 - CY0) / 2}"/></clipPath>
    <clipPath id="card"><rect width="{W}" height="{H}" rx="24"/></clipPath>
    <clipPath id="term"><rect x="92" y="300" width="560" height="54" rx="12"/></clipPath>
    <pattern id="dots" width="22" height="22" patternUnits="userSpaceOnUse"><circle cx="2" cy="2" r="1.1" fill="{ACC}" fill-opacity=".07"/></pattern>
  </defs>
  <g clip-path="url(#card)">
    <rect width="{W}" height="{H}" fill="url(#bg)"/><rect width="{W}" height="{H}" fill="url(#dots)"/><rect width="{W}" height="{H}" fill="url(#glow)"/>
    <rect x="0" y="{H - 4}" width="{W}" height="4" fill="{LINE}"/>
    <rect class="prog" x="0" y="{H - 4}" width="{W}" height="4" fill="{HOT}" fill-opacity=".7"/>
  </g>
  <rect x=".5" y=".5" width="{W - 1}" height="{H - 1}" rx="24" fill="none" stroke="#24434a"/>
  <rect x="700" y="30" width="552" height="420" rx="18" fill="{PANEL}" fill-opacity=".55" stroke="#24434a"/>

  <g class="rise"><text x="92" y="190" class="sans" font-size="104" font-weight="800" fill="#eef6f3" letter-spacing="-3">AMRScan</text></g>
  <g class="rise d1"><text x="94" y="236" class="sans" font-size="25" font-weight="600" fill="#86d6c8">Resistance calls you can trace back to the run,</text>
    <text x="94" y="268" class="sans" font-size="25" font-weight="600" fill="#86d6c8">tested on populations the model has not seen.</text></g>
  <g class="rise d2">
    <rect x="92" y="300" width="560" height="54" rx="12" fill="#0b161c" stroke="#24434a"/>
    <text x="110" y="{TY}" class="mono" font-size="18" fill="{HOT}">›</text>
    <g clip-path="url(#term)">{''.join(term)}</g>
  </g>
  <text x="94" y="404" class="sans" font-size="16" fill="#7f9a96">Nextflow · NCBI AMRFinderPlus · 10,521 <tspan font-style="italic">E. coli</tspan> genomes · 26 countries</text>

  {s0}
  {s1}
  {s2}
  {s3}
</svg>
"""
theme = sys.argv[2] if len(sys.argv) > 2 else "teal"
m = dict(zip(KEYS, THEMES[theme]))
svg = re.sub(r"#[0-9a-f]{6}", lambda x: m.get(x.group(0), x.group(0)), svg)
open(sys.argv[1], "w").write(svg)
print("wrote", sys.argv[1], len(svg), "bytes")

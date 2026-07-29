#!/usr/bin/env python3
"""Convert $$...$$ display math in TEP-BH site components to numbered
\\begin{equation}\\label{eq:...}\\end{equation} blocks wrapped in
<div class="equation">, matching the TEP-COS / TEP-KIN convention.

- Each component file gets its own per-file sequential counter.
- Equations inside <p>...</p> split the paragraph into <p>text</p>,
  <div class="equation">...</div>, <p>text</p> chains (block-level div
  cannot live inside <p>).
- Equations outside <p> are wrapped directly.
- Empty <p></p> fragments produced by the split are dropped.
"""

import os
import re

COMPONENTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "site", "components",
)

TAG_PREFIX = {
    "0_abstract.html": "abstract",
    "01_introduction.html": "intro",
    "02_tep_theory.html": "tep",
    "03_fixed_background_theorem.html": "fbt",
    "04_coupled_solution.html": "coupled",
    "05_global_geometry.html": "global",
    "06_causal_structure.html": "causal",
    "07_stability.html": "stab",
    "08_predictions.html": "pred",
    "09_interpretation.html": "interp",
    "10_conclusion.html": "concl",
    "17_references.html": "ref",
    "18_acknowledgements.html": "ack",
    "A_appendix_conventions.html": "appA",
    "B_appendix_field_equations.html": "appB",
    "C_appendix_asymptotics.html": "appC",
    "D_appendix_curvature.html": "appD",
    "E_appendix_geodesics.html": "appE",
    "F_appendix_volume.html": "appF",
    "G_appendix_perturbations.html": "appG",
    "H_appendix_raytracing.html": "appH",
    "I_appendix_comparison.html": "appI",
    "J_appendix_reproducibility.html": "appJ",
    "K_appendix_hayward_benchmark.html": "appK",
    "L_appendix_corrected_observables.html": "appL",
}

EQN_RE = re.compile(r"\$\$(.*?)\$\$", re.DOTALL)
P_RE = re.compile(r"<p>(.*?)</p>", re.DOTALL)


def make_block(content, label):
    content = content.strip()
    return (
        '<div class="equation">\n'
        f'    \\begin{{equation}} \\label{{{label}}}\n'
        f'    {content}\n'
        '\\end{equation}\n'
        '</div>'
    )


def convert_text(text, prefix, counter):
    """Replace all $$...$$ in `text` with numbered div blocks."""
    def repl(m):
        counter[0] += 1
        label = f"eq:{prefix}_{counter[0]}"
        return make_block(m.group(1), label)
    return EQN_RE.sub(repl, text)


def convert_file(path, prefix):
    with open(path, "r", encoding="utf-8") as f:
        html = f.read()

    counter = [0]

    def replace_in_p(m):
        p_content = m.group(1)
        if "$$" not in p_content:
            return m.group(0)
        parts = EQN_RE.split(p_content)  # text, eqn, text, eqn, ...
        if len(parts) == 1:
            return m.group(0)
        pieces = []
        for i, part in enumerate(parts):
            if i % 2 == 0:
                txt = part.strip()
                if txt:
                    pieces.append(f"<p>\n{txt}\n</p>")
            else:
                counter[0] += 1
                label = f"eq:{prefix}_{counter[0]}"
                pieces.append(make_block(part, label))
        return "\n\n".join(pieces)

    new_html = P_RE.sub(replace_in_p, html)

    # Catch any $$...$$ that lived outside <p> (e.g. bare in a <div>).
    new_html = convert_text(new_html, prefix, counter)

    with open(path, "w", encoding="utf-8") as f:
        f.write(new_html)
    return counter[0]


def main():
    total = 0
    for fname in sorted(os.listdir(COMPONENTS_DIR)):
        if not fname.endswith(".html"):
            continue
        prefix = TAG_PREFIX.get(fname)
        if prefix is None:
            print(f"  skip (no prefix mapping): {fname}")
            continue
        path = os.path.join(COMPONENTS_DIR, fname)
        n = convert_file(path, prefix)
        total += n
        print(f"  {fname:42s} -> {n:3d} equations  (labels eq:{prefix}_1..{prefix}_{n})")
    print(f"TOTAL: {total} equations converted")


if __name__ == "__main__":
    main()

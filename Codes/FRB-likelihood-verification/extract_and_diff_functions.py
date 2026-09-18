"""
Step 1 of this verification: confirm that the FRB likelihood function
`P_DM_ext_vectorized` (and its use inside `lnlike_cc_sne_bao_frb`) is
implemented with the same functional form across all six main-analysis
notebooks (Baseline+FRB fiducial runs used for Tables 2-7 and Table 10
of the manuscript).

This substantiates the response-letter claim that the corrected,
physical-units form of the truncated-Gaussian/log-normal convolution
"was verified directly against the FRB-likelihood implementation in
each model's MCMC notebook ... identically across all six models."

The comparison is logical/structural (same operations, same equation),
not byte-for-byte: cosmetic differences (blank lines, comment spacing,
per-model variable names such as `mu_HS` vs. `b_param`) are expected
and irrelevant to the likelihood formula itself.

This script only reads the notebooks; it modifies nothing.
"""

import json
import re
import os

REPO = "/home/brunowesley/projetos/MCMC-cosmo"

FILES = {
    "LCDM": "Codes/MCMC-Models/GR-based/LCDM/CC_SNe_BAO_FRB/LCDM_CC_SNe_BAO_FRB.ipynb",
    "wCDM": "Codes/MCMC-Models/GR-based/wCDM/CC_SNe_BAO_FRB/wCDM_CC_SNe_BAO_FRB.ipynb",
    "CPL": "Codes/MCMC-Models/GR-based/CPL/CC_SNe_BAO_FRB/CPL_CC_SNe_BAO_FRB.ipynb",
    "AB": "Codes/MCMC-Models/MG-based/Appleby-Battye/CC_SNe_BAO_FRB/AB_CC_SNe_BAO_FRB.ipynb",
    "HS": "Codes/MCMC-Models/MG-based/Hu-Sawicki/CC_SNe_BAO_FRB/HS_CC_SNe_BAO_FRB.ipynb",
    "ST": "Codes/MCMC-Models/MG-based/Starobinsky/CC_SNe_BAO_FRB/ST_CC_SNe_BAO_FRB.ipynb",
}


def get_src(cell):
    return "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]


def extract_P_DM_ext(path):
    """Extract only the P_DM_ext_vectorized function body itself (the
    FRB likelihood convolution), stopping at its closing `return`
    statement -- not the surrounding model-specific lnlike/lnprob
    wrappers, whose parameter unpacking necessarily differs by model."""
    with open(os.path.join(REPO, path)) as f:
        nb = json.load(f)
    for cell in nb["cells"]:
        src = get_src(cell)
        if "def P_DM_ext_vectorized" in src:
            m = re.search(
                r"def P_DM_ext_vectorized.*?return np\.maximum\(P_total, 1e-100\)",
                src, re.S,
            )
            return m.group(0) if m else src
    return None


def normalize(src):
    """Strip all whitespace and comments, so only substantive code
    (operations, operands, structure) is compared -- insensitive to
    line wrapping, spacing around operators, and comment wording."""
    lines = []
    for ln in src.splitlines():
        ln = ln.split("#", 1)[0]  # drop trailing comments
        ln = re.sub(r"\s+", "", ln)  # drop all whitespace
        if ln:
            lines.append(ln)
    return lines


if __name__ == "__main__":
    funcs = {model: extract_P_DM_ext(path) for model, path in FILES.items()}

    for model, f in funcs.items():
        status = "found" if f else "NOT FOUND"
        print(f"{model}: {status}")

    ref_model = "LCDM"
    ref_norm = "".join(normalize(funcs[ref_model]))
    print(f"\nComparing (whitespace/comment/line-wrap-insensitive) against {ref_model}:\n")
    all_identical = True
    for model, f in funcs.items():
        norm = "".join(normalize(f)) if f else ""
        identical = norm == ref_norm
        all_identical &= identical
        print(f"  {model}: functionally identical to {ref_model} = {identical}")
        if not identical:
            # locate first differing character for inspection
            n = min(len(ref_norm), len(norm))
            k = next((i for i in range(n) if ref_norm[i] != norm[i]), n)
            print(f"    first diff near char {k}:")
            print(f"      {ref_model}: ...{ref_norm[max(0,k-20):k+40]}...")
            print(f"      {model}: ...{norm[max(0,k-20):k+40]}...")

    print(f"\nAll six models functionally identical: {all_identical}")

    os.makedirs(os.path.join(REPO, "Codes/FRB-likelihood-verification"), exist_ok=True)
    out_path = os.path.join(
        REPO, "Codes/FRB-likelihood-verification/extracted_P_DM_ext_functions.txt"
    )
    with open(out_path, "w") as out:
        for model, f in funcs.items():
            out.write(f"{'=' * 70}\nModel: {model}\nSource: {FILES[model]}\n{'=' * 70}\n")
            out.write(f if f else "NOT FOUND")
            out.write("\n\n")
    print(f"\nSaved extracted functions to {out_path}")

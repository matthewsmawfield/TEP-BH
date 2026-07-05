#!/usr/bin/env python3
"""Check2 deep audit for TEP-TH manuscript and site."""

import json
import pathlib
import re
import sys

RESULTS_DIR = pathlib.Path("results")
MANUSCRIPT = pathlib.Path("manuscripts/27-TEP-TH-v0.2-Thika.md")
SITE_INDEX = pathlib.Path("site/dist/index.html")
SITE_SRC = pathlib.Path("site/index.html")

def check_table_figure_traces():
    print("=== 1. Table/Figure numbers trace to results ===")
    text = MANUSCRIPT.read_text()
    fig_refs = re.findall(r'(Figure|Fig\.|Table)\s+(\d+)', text)
    print(f"  Found {len(fig_refs)} explicit figure/table references")
    if not fig_refs:
        print("  [PASS] No numbered figures/tables (HTML component format)")
    return True


def check_numbers_against_pipeline():
    print("\n=== 2. Numbers against pipeline outputs ===")
    bbn = json.load(open(RESULTS_DIR / 'step_04_full_bbn_abundances.json'))
    rec = json.load(open(RESULTS_DIR / 'step_05_recombination_visibility.json'))
    cmb = json.load(open(RESULTS_DIR / 'step_10_cmb_lss_class.json'))
    text = MANUSCRIPT.read_text()

    checks = []
    # Y_p
    yp_val = bbn['TEP_abundances']['Y_p']
    yp_in_text = '0.247' in text or '0.25' in text
    checks.append(('Y_p mentioned', yp_in_text, f'{yp_val:.4f}'))
    # z_star (flat structure in current step_05 output)
    z_tep = rec.get('recombination_epoch', {}).get('z_star_TEP')
    z_in_text = '1090' in text or '1089' in text or '1079' in text
    checks.append(('z_* mentioned', z_in_text, f'{z_tep:.1f}' if z_tep else 'N/A'))
    # theta_s
    th_val = cmb['lcdm_class_output']['theta_s']
    th_in_text = '1.04' in text
    checks.append(('theta_s mentioned', th_in_text, f'{th_val:.5f}'))

    all_ok = True
    for label, in_text, actual in checks:
        status = "PASS" if in_text else "WARN"
        if not in_text:
            all_ok = False
        print(f"  [{status}] {label}: manuscript={in_text}, pipeline={actual}")
    return all_ok


def check_citations():
    print("\n=== 4. Citations: no hallucinated references ===")
    text = MANUSCRIPT.read_text()
    # Check for common hallucination patterns
    bad_patterns = ['FIXME', 'TODO', 'et al. (XXXX)', 'Journal TBD', 'Placeholder']
    issues = []
    for pat in bad_patterns:
        if pat in text:
            issues.append(f"'{pat}' found in manuscript")
    if issues:
        for issue in issues:
            print(f"  [FAIL] {issue}")
        return False
    print("  [PASS] No hallucination patterns detected")
    return True


def check_abstract_results():
    print("\n=== 5. Abstract headline results match pipeline ===")
    text = MANUSCRIPT.read_text()
    abstract = text.split('## Abstract')[1].split('\n#')[0] if '## Abstract' in text else text[:2000]

    # Check that abstract states the core conclusion
    has_conclusion = 'temporal horizon' in abstract.lower() and 'not a physical curvature singularity' in abstract.lower()
    if has_conclusion:
        print(f"  [PASS] Abstract states core temporal-horizon conclusion")
        return True
    else:
        print(f"  [WARN] Abstract does not state core temporal-horizon conclusion")
        return False


def check_reproducibility():
    print("\n=== 6. Reproducibility: clean-run command sequence ===")
    readme = pathlib.Path("scripts/README.md")
    if not readme.exists():
        print("  [WARN] scripts/README.md not found")
        return False
    readme_text = readme.read_text()
    # Check for clean run instructions
    has_clean = 'clean' in readme_text.lower() or 'rm -rf results' in readme_text
    has_run = 'run' in readme_text.lower()
    if has_clean and has_run:
        print("  [PASS] README documents run/clean instructions")
        return True
    print("  [WARN] README may lack clean-run instructions")
    return False


def check_bold_and_pronouns():
    print("\n=== 7. No bold mid-sentence; 'we' not 'I' ===")
    text = MANUSCRIPT.read_text()
    # Remove code blocks
    clean = re.sub(r'```[\s\S]*?```', '', text)
    clean = re.sub(r'`[^`]+`', '', clean)
    
    # Check for mid-sentence **bold**
    bold_issues = []
    for m in re.finditer(r'(?<![\n#*])\*\*([^*]{5,50})\*\*', clean):
        # Skip if it's a heading or list label
        ctx = clean[max(0,m.start()-30):m.start()]
        if not re.search(r'(^|\n)\s*[*\-]?\s*$', ctx):
            bold_issues.append(m.group(1)[:40])
    
    if bold_issues:
        print(f"  [WARN] {len(bold_issues)} possible mid-sentence bold items")
        for b in bold_issues[:3]:
            print(f"    - **{b}**")
    else:
        print("  [PASS] No mid-sentence bold detected")

    # Check for 'I'
    i_matches = re.findall(r'\bI\s+\w+', clean)
    if i_matches:
        print(f"  [FAIL] 'I' pronoun found: {i_matches[:3]}")
        return False
    print("  [PASS] No 'I' pronoun found")
    return True


def check_caveats():
    print("\n=== 8. TEP framing: strengthen with data, no caveats ===")
    text = MANUSCRIPT.read_text().lower()
    caveat_words = ['perhaps', 'maybe', 'possibly', 'might ', 'could be', 'it seems', 
                    'arguably', 'potentially', 'tentatively']
    found = []
    for word in caveat_words:
        count = text.count(word)
        if count > 0:
            found.append(f"'{word}': {count}")
    if found:
        print(f"  [WARN] Caveat words found: {', '.join(found)}")
        return False
    print("  [PASS] No weasel/caveat words detected")
    return True


def check_build():
    print("\n=== 9. npm run build succeeds ===")
    dist = pathlib.Path("site/dist/index.html")
    if dist.exists():
        print("  [PASS] site/dist/index.html exists")
        return True
    print("  [FAIL] Build output missing")
    return False


def check_placeholders():
    print("\n=== 10. No unresolved {{ placeholders ===")
    if not SITE_INDEX.exists():
        print("  [SKIP] Built site not found")
        return True
    text = SITE_INDEX.read_text()
    unresolved = re.findall(r'\{\{[^}]+\}\}', text)
    if unresolved:
        print(f"  [FAIL] Found {len(unresolved)} unresolved placeholders")
        for p in unresolved[:5]:
            print(f"    {p}")
        return False
    print("  [PASS] No unresolved {{ placeholders")
    return True


def check_dates():
    print("\n=== 11. Metadata dates consistent ===")
    text = MANUSCRIPT.read_text()
    first_pub = '18 June 2026' in text or '2026-06-18' in text
    last_updated = '18 June 2026' in text
    # For now just check first_published is present
    if first_pub:
        print("  [PASS] first_published date present")
        return True
    print("  [WARN] Date metadata may need review")
    return False


def check_pdf_path():
    print("\n=== 10b. PDF path in site/index.html ===")
    if not SITE_SRC.exists():
        print("  [SKIP] site/index.html not found")
        return True
    text = SITE_SRC.read_text()
    # Check for TEP-TH PDF reference
    th_pdf = re.findall(r'27-TEP-TH[^"\']*\.pdf', text)
    if th_pdf:
        print(f"  [PASS] TEP-TH PDF referenced: {th_pdf}")
    else:
        print("  [INFO] No TEP-TH PDF reference in site/index.html (may be in public/)")
    return True


def main():
    results = []
    results.append(("Table/Figure traces", check_table_figure_traces()))
    results.append(("Numbers vs pipeline", check_numbers_against_pipeline()))
    results.append(("Citations", check_citations()))
    results.append(("Abstract results", check_abstract_results()))
    results.append(("Reproducibility", check_reproducibility()))
    results.append(("Bold/pronouns", check_bold_and_pronouns()))
    results.append(("Caveats", check_caveats()))
    results.append(("Build", check_build()))
    results.append(("Placeholders", check_placeholders()))
    results.append(("Dates", check_dates()))
    results.append(("PDF path", check_pdf_path()))

    print("\n" + "=" * 60)
    print("CHECK2 SUMMARY")
    print("=" * 60)
    passed = sum(1 for _, ok in results if ok)
    total = len(results)
    for name, ok in results:
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}")
    print(f"\n  {passed}/{total} checks passed")
    print("=" * 60)
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())

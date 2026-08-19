---
name: regen-report
description: Regenerate the step-2 HTML report by running phylohpc report
---

Regenerate the step-2 HTML report.

Steps:
1. Run `phylohpc report --results_dir results --species_tree data/species_tree.full.newick --family_info genefam.csv --output report2.html` from the project root
2. If it succeeds, run `stat -c '%s bytes' report2.html` and report the size to the user
3. If it fails, show the error output so the user can diagnose it

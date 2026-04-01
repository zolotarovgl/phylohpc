---
name: regen-report
description: Regenerate the step-2 HTML report by running make_report.sh
---

Regenerate the step-2 HTML report.

Steps:
1. Run `bash make_report.sh` from the project root
2. If it succeeds, run `stat -c '%s bytes' report2.html` and report the size to the user
3. If it fails, show the error output so the user can diagnose it

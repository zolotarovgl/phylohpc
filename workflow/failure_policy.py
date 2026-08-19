#!/usr/bin/env python3
"""Failure-accounting policy for step2.nf (family failures from ANY of its processes --
ALN/PHY/PVM/PVM_PREV/GR_watcher, not GeneRax only, since 2026-08-19).

Two pure functions, both unit-testable without Nextflow:

  - tolerance_limit(n_families, max_failed_frac, max_failed_count):
        the number of failures the run tolerates before it is FAILED.
  - decide(n_failed, n_families, max_failed_frac, max_failed_count, strict):
        'ok' | 'fail' -- whether the run passes the tolerance.
  - write_failures_tsv(path, rows): the exact 5-column
        family\texit_code\tfirst_error\tworkdir\tprocess format. The `process` column
        (added 2026-08-19) records WHICH stage failed -- ALN/PHY/PVM/PVM_PREV/GR_watcher
        -- because before this every failure looked like a GeneRax failure by omission.

Also usable as a CLI so step2.nf can shell out to the SAME decision rather than
re-implementing the arithmetic in Groovy (two copies of a threshold is how it drifts):

    python3 failure_policy.py decide N_FAILED N_FAMILIES MAX_FRAC MAX_COUNT STRICT
        -> prints "ok" or "fail" and exits 0 either way; the caller reads stdout.

    python3 failure_policy.py write-tsv OUTPATH < rows.tsv
        -> rows.tsv is family\texit_code\tfirst_error\tworkdir\tprocess, no header;
           writes OUTPATH with the header prepended.
"""
import math
import sys


def tolerance_limit(n_families, max_failed_frac, max_failed_count):
    """The max number of failures a run of n_families tolerates."""
    if n_families <= 0:
        return 0
    frac_limit = math.ceil(n_families * max_failed_frac)
    return min(int(max_failed_count), int(frac_limit))


def decide(n_failed, n_families, max_failed_frac, max_failed_count, strict=False):
    """Return 'ok' or 'fail'."""
    if n_failed <= 0:
        return 'ok'
    if strict:
        return 'fail'
    limit = tolerance_limit(n_families, max_failed_frac, max_failed_count)
    return 'ok' if n_failed <= limit else 'fail'


HEADER = "family\texit_code\tfirst_error\tworkdir\tprocess"


def write_failures_tsv(path, rows):
    """rows: iterable of (family, exit_code, first_error, workdir, process) tuples/lists."""
    lines = [HEADER]
    for family, exit_code, first_error, workdir, process in rows:
        lines.append(f"{family}\t{exit_code}\t{first_error}\t{workdir}\t{process}")
    text = "\n".join(lines) + "\n"
    with open(path, "w") as fh:
        fh.write(text)
    return text


def _main(argv):
    if len(argv) < 2:
        print("usage: failure_policy.py decide|write-tsv ...", file=sys.stderr)
        return 2

    cmd = argv[1]
    if cmd == 'decide':
        n_failed, n_families = int(argv[2]), int(argv[3])
        max_frac, max_count = float(argv[4]), int(argv[5])
        strict = argv[6].lower() in ('1', 'true', 'yes') if len(argv) > 6 else False
        print(decide(n_failed, n_families, max_frac, max_count, strict))
        return 0
    elif cmd == 'write-tsv':
        outpath = argv[2]
        rows = []
        for line in sys.stdin:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 5:
                print(f"malformed row (expected 5 fields, got {len(parts)}): {line!r}", file=sys.stderr)
                return 1
            rows.append(parts)
        write_failures_tsv(outpath, rows)
        return 0
    else:
        print(f"unknown command: {cmd}", file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(_main(sys.argv))

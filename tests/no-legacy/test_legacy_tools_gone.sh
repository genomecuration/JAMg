#!/usr/bin/env bash
# Phase 5 static-policy test. Two-prong:
#   1. Every path enumerated in forbidden_files.txt must be absent from HEAD.
#   2. No production code path (bin/, PerlLib/, workflow/, tests/) contains
#      references to the dropped predictors (SNAP / geneid / GlimmerHMM /
#      HHblits / fathom / .zff). Allowlists noted inline.
#
# Per plan §5.1 (review-#2 I16): cwd-aware, runs from repo root.

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

fail=0
FORBIDDEN=tests/no-legacy/forbidden_files.txt

# 1. orphan paths: must not exist
while IFS= read -r f; do
    [[ -z "$f" || "$f" =~ ^# ]] && continue
    if [[ -e "$f" ]]; then
        echo "STILL PRESENT: $f"
        fail=1
    fi
done < "$FORBIDDEN"

# 2. residual code references. Production paths only.
#    - bin/ + PerlLib/ + tests/ scanned (dropped predictors must not live here)
#    - workflow/rules/rnaseq.smk legitimately references gsnap (optional RNA aligner)
#    - Augustus's --uniqueGeneId CLI flag + renamed $gene_identifier variable allowlisted
#    - fathom + .zff caught here (Critical-1 + remove-fathom directive)
hits=$(
    grep -rIn --include='*.pl' --include='*.pm' --include='*.py' --include='*.smk' --include='*.sh' \
      -E '\b(SNAP|snap|GeneID|geneid|GlimmerHMM|glimmer|HHblits|hhblits|HHsuite|hhsuite|HHsearch|hhsearch|fathom|\.zff\b)\b' \
      bin/ PerlLib/ \
      $(find tests -type f \( -name '*.pl' -o -name '*.pm' -o -name '*.py' -o -name '*.smk' -o -name '*.sh' -o -name '*.t' \) \
         ! -path 'tests/no-legacy/*' ! -path 'tests/e2e/*') \
      workflow/Snakefile workflow/rules/common.smk workflow/rules/repeats.smk \
      workflow/rules/tgg.smk workflow/rules/pasa.smk workflow/rules/golden.smk \
      workflow/rules/genemark.smk workflow/rules/proteins.smk workflow/rules/augustus.smk \
      workflow/rules/evm.smk workflow/rules/ogs.smk \
      2>/dev/null \
      | grep -vE ':\s*#' \
      | grep -vE 'uniqueGeneId|gene_identifier' \
      || true
)
# Separate scan: gsnap token allowed in workflow/rules/rnaseq.smk ONLY.
gsnap_leak=$(
    grep -rIn --include='*.pl' --include='*.pm' --include='*.py' --include='*.smk' --include='*.sh' \
      -E '\b(GSNAP|gsnap)\b' \
      bin/ PerlLib/ workflow/Snakefile workflow/rules/common.smk \
      $(ls workflow/rules/*.smk | grep -v rnaseq.smk) \
      2>/dev/null \
      | grep -vE ':\s*#' \
      || true
)
hits="$hits$gsnap_leak"
if [[ -n "$hits" ]]; then
    echo "RESIDUAL REFS FOUND:"
    echo "$hits"
    fail=1
fi
exit $fail

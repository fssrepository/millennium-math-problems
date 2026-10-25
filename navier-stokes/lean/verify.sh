#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

lake_bin="${LAKE:-lake}"
if ! command -v "$lake_bin" >/dev/null 2>&1; then
  if [[ -z "${LAKE:-}" && -x "$HOME/.elan/bin/lake" ]]; then
    lake_bin="$HOME/.elan/bin/lake"
  else
    echo 'Lake not found. Install elan using the command in README.md.' >&2
    exit 2
  fi
fi
mkdir -p .verification

# Every proof source must be imported by the audited root. An unimported new
# file must not silently escape the build or axiom audit.
python3 - <<'PY'
from pathlib import Path
import re
sources = {str(p.with_suffix('')).replace('/', '.')
           for p in Path('NavierStokesAudit').rglob('*.lean')}
imports = set(re.findall(r'^public import (NavierStokesAudit\.[\w.]+)$',
                        Path('NavierStokesAudit.lean').read_text(), re.M))
if sources != imports:
    raise SystemExit(f'Proof import inventory mismatch: missing={sources-imports}, stale={imports-sources}')
print(f'Proof inventory PASS: {len(sources)} modules')
PY

# The manifest pins every dependency. Do not run `lake update` during verification.
# Fetch only the mathlib modules imported by this project, not all of mathlib.
mapfile -t mathlib_modules < <(python3 - <<'PY'
from pathlib import Path
import re
modules = set()
for p in [Path('NavierStokesAudit.lean'), *Path('NavierStokesAudit').rglob('*.lean')]:
    modules.update(re.findall(r'^(?:public )?(?:meta )?import (Mathlib\.[\w.]+)$', p.read_text(), re.M))
print('\n'.join(sorted(modules)))
PY
)
MATHLIB_NO_CACHE_ON_UPDATE=1 "$lake_bin" exe cache get "${mathlib_modules[@]}"
"$lake_bin" build
"$lake_bin" env lean -DwarningAsError=true tests/lean/Audit.lean | tee .verification/audit.log

expect_rejection() {
  local source="$1" diagnostic="$2" log_file=".verification/$(basename "$1").log"
  if "$lake_bin" env lean -DwarningAsError=true "$source" >"$log_file" 2>&1; then
    echo "FAIL: Lean accepted the negative regression $source" >&2
    exit 1
  fi
  if ! grep -Fq -- "$diagnostic" "$log_file"; then
    cat "$log_file" >&2
    echo "FAIL: $source failed for an unexpected reason" >&2
    exit 1
  fi
  echo "Expected rejection PASS: $source"
}

expect_rejection tests/lean/rejections/FalseBound.lean 'unsolved goals'
expect_rejection tests/lean/rejections/FalseMajorantBound.lean 'unsolved goals'
expect_rejection tests/lean/rejections/Unfinished.lean 'declaration uses'
expect_rejection tests/lean/rejections/ExtraAxiom.lean 'PNT audit: unexpected axiom'
bash scripts/check-cpp.sh
echo 'MATH-19 PASS: unrestricted PNT-2/7/8/12/14 bounds disproved; proof audit and C++ consistency checks passed.'

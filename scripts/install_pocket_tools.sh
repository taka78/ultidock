#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 1 ]]; then
  printf 'Usage: %s [fpocket|p2rank|all]\n' "$0" >&2
  exit 2
fi
method="${1:-all}"
case "$method" in
  fpocket|p2rank|all) ;;
  *) printf 'Unknown pocket tool: %s (expected fpocket, p2rank, or all)\n' "$method" >&2; exit 2 ;;
esac

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
external="$root/external"
mkdir -p "$external/bin" "$external/tools"

if [[ "$method" == fpocket || "$method" == all ]]; then
  if [[ ! -x "$external/fpocket/bin/fpocket" ]]; then
    if [[ ! -d "$external/fpocket" ]]; then
      git clone --depth 1 --branch 4.2.3 https://github.com/Discngine/fpocket.git "$external/fpocket"
    fi
    # fpocket 4.2.3 declares an array of pointers where its parser needs a
    # character buffer. GCC 15 rejects the resulting strcpy argument type.
    python3 - "$external/fpocket/src/fparams.c" <<'PY'
from pathlib import Path
import sys

source = Path(sys.argv[1])
original = source.read_text(encoding="utf-8")
updated = original
for old, new in (
    ("char *residue_string[M_MAX_CUSTOM_POCKET_LEN];", "char residue_string[M_MAX_CUSTOM_POCKET_LEN];"),
    ("strcpy(&residue_string, pt);", "strcpy(residue_string, pt);"),
):
    if old in updated:
        updated = updated.replace(old, new, 1)
    elif new not in updated:
        raise SystemExit(f"Unexpected fpocket 4.2.3 source in {source}: missing {old!r}")
if updated != original:
    source.write_text(updated, encoding="utf-8")
PY
    make -C "$external/fpocket"
  fi
  printf 'fpocket: %s\n' "$external/fpocket/bin/fpocket"
fi

if [[ "$method" == p2rank || "$method" == all ]]; then
  if [[ ! -x "$external/tools/p2rank_2.5/prank" ]]; then
    archive="$external/tools/p2rank_2.5.tar.gz"
    curl -fL https://github.com/rdk/p2rank/releases/download/2.5/p2rank_2.5.tar.gz -o "$archive"
    tar -xzf "$archive" -C "$external/tools"
  fi

  cat > "$external/bin/prank" <<'WRAPPER'
#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../tools/p2rank_2.5" && pwd)"
if [[ -z "${JAVA_HOME:-}" && -x /usr/lib/jvm/java-21-openjdk-amd64/bin/java ]]; then
  export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
fi
exec "$root/prank" "$@"
WRAPPER
  chmod +x "$external/bin/prank"

  printf 'P2Rank: %s\n' "$external/bin/prank"
fi

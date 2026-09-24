#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
external="$root/external"
mkdir -p "$external/bin" "$external/tools"

if [[ ! -x "$external/fpocket/bin/fpocket" ]]; then
  if [[ ! -d "$external/fpocket" ]]; then
    git clone --depth 1 --branch 4.2.3 https://github.com/Discngine/fpocket.git "$external/fpocket"
  fi
  make -C "$external/fpocket"
fi

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

printf 'fpocket: %s\nP2Rank: %s\n' "$external/fpocket/bin/fpocket" "$external/bin/prank"

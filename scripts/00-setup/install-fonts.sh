#!/usr/bin/env bash
# Install the figure font, TeX Gyre Heros 2.004 (Helvetica metrics; GUST e-foundry, GUST Font License), for the
# current user. make-figures.R uses it when installed and falls back to the default sans font otherwise.
#
#   pixi run fonts
#
# The zip is checked against its sha256 (recorded 2026-11-17; the two .otf files it gives are identical to the ones
# the published figures were drawn with). HTTPS from some clusters is reset during the TLS handshake, so curl retries.
set -euo pipefail
URL=https://www.gust.org.pl/projects/e-foundry/tex-gyre/heros/qhv2.004otf.zip
SHA256=755954b7349b803fc1c3d82fe9d9c4f7cf66467af718eaaf4f78ae1a09bf265d
DEST=${XDG_DATA_HOME:-$HOME/.local/share}/fonts

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
for _ in $(seq 1 30); do
  curl -sSfL --connect-timeout 15 --retry 5 --retry-delay 1 --retry-all-errors -o "$tmp/qhv.zip" "$URL" && break
  sleep 1
done
[[ -s $tmp/qhv.zip ]] || { echo "error: could not download $URL" >&2; exit 1; }
echo "$SHA256  $tmp/qhv.zip" | sha256sum -c --quiet || { echo "error: checksum mismatch for $URL" >&2; exit 1; }

mkdir -p "$DEST"
# unzip is not always installed; Python's zipfile is
python3 -I - "$tmp/qhv.zip" "$DEST" <<'EOF'
import sys, zipfile, pathlib
z, dest = zipfile.ZipFile(sys.argv[1]), pathlib.Path(sys.argv[2])
for name in z.namelist():
    base = name.rsplit("/", 1)[-1]
    if base.startswith("texgyreheros-") and base.endswith(".otf"):   # regular, bold, italic, bolditalic
        (dest / base).write_bytes(z.read(name)); print("installed", dest / base)
EOF
command -v fc-cache >/dev/null && fc-cache -f "$DEST" >/dev/null
fc-list : family 2>/dev/null | grep -q "TeX Gyre Heros" && echo "TeX Gyre Heros is available to R and Graphviz."

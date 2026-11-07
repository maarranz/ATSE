#!/usr/bin/env bash
# Runs after every project render: collects the Typst PDFs into docs/handouts/.
# Reports what it did; never stops the render.
out="${QUARTO_PROJECT_OUTPUT_DIR:-docs}"
mkdir -p "$out/handouts"
n=0
for dir in "$out/typst" "$out/chapters"; do
  [ -d "$dir" ] || continue
  while IFS= read -r -d '' f; do
    cp "$f" "$out/handouts/" && n=$((n + 1))
  done < <(find "$dir" -name '*.pdf' -print0)
done
if [ "$n" -eq 0 ]; then
  echo "copy-pdfs: WARNING: no PDFs found under $out/typst or $out/chapters" >&2
else
  echo "copy-pdfs: copied $n PDF(s) to $out/handouts/"
fi
exit 0

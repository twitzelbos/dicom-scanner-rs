#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

DRY_RUN=0

get_tag() {
  local tag="$1"
  local file="$2"
  dcmdump +P "$tag" "$file" 2>/dev/null | awk -F'[][]' '{print $2}' | tr -d '[:space:]'
}

for dir in */; do
  echo "==> Processing: $dir"

  # collect candidate files (skip DICOMDIR)
  files=()
  for f in "$dir"*; do
    [ -f "$f" ] || continue
    [ "$(basename "$f")" = "DICOMDIR" ] && continue
    files+=("$f")
  done
  [ "${#files[@]}" -eq 0 ] && continue

  mapfile="$(mktemp)"
  trap 'rm -f "$mapfile"' EXIT

  # Build mapping: src|final
  # final must be <= 8 chars: IM + 6 digits
  used="$(mktemp)"
  trap 'rm -f "$mapfile" "$used"' EXIT

  i=1
  for f in "${files[@]}"; do
    inst="$(get_tag "0020,0013" "$f")"
    # default to counter if missing or non-numeric
    if [[ "${inst:-}" =~ ^[0-9]+$ ]]; then
      n="$inst"
    else
      n="$i"
    fi

    # try IM%06d; if collision, bump until free
    while :; do
      final=$(printf "IM%06d" "$n")
      if ! grep -qx "$final" "$used" 2>/dev/null; then
        echo "$final" >> "$used"
        break
      fi
      n=$((n+1))
    done

    printf "%s|%s\n" "$f" "$final" >> "$mapfile"
    i=$((i+1))
  done

  # Phase 1: temp names
  while IFS='|' read -r src final; do
    tmp="$dir.tmp_$(basename "$src")"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "  mv -- '$src' '$tmp'"
    else
      mv -- "$src" "$tmp"
    fi
  done < "$mapfile"

  # Phase 2: temp -> final (8-char)
  while IFS='|' read -r src final; do
    tmp="$dir.tmp_$(basename "$src")"
    dest="$dir$final"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "  mv -- '$tmp' '$dest'"
    else
      mv -- "$tmp" "$dest"
    fi
  done < "$mapfile"

  rm -f "$mapfile" "$used"
done

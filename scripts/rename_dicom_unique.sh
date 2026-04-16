#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

DRY_RUN=0

get_tag() {
  local tag="$1"
  local file="$2"
  dcmdump +P "$tag" "$file" 2>/dev/null | awk -F'[][]' '{print $2}' | tr -d '[:space:]'
}

sanitize_uid() {
  # UIDs are already [0-9.] but trim just in case
  echo "$1" | tr -cd '0-9.'
}

for dir in */; do
  echo "==> Processing: $dir"

  mapfile="$(mktemp)"
  trap 'rm -f "$mapfile"' EXIT

  for f in "$dir"*; do
    [ -f "$f" ] || continue
    [ "$(basename "$f")" = "DICOMDIR" ] && continue

    inst="$(get_tag "0020,0013" "$f")"
    uid="$(get_tag "0008,0018" "$f")"   # SOP Instance UID
    uid="$(sanitize_uid "$uid")"

    [ -z "${uid:-}" ] && { echo "  Skipping (no SOPInstanceUID): $(basename "$f")"; continue; }
    [ -z "${inst:-}" ] && inst="0"

    final=$(printf "IMG%04d_%s.dcm" "$inst" "$uid")
    printf "%s|%s\n" "$f" "$final" >> "$mapfile"
  done

  # Detect filename collisions (should be extremely unlikely now)
  if awk -F'|' '{print $2}' "$mapfile" | sort | uniq -d | grep -q .; then
    echo "ERROR: filename collision even with SOPInstanceUID in $dir"
    rm -f "$mapfile"
    exit 1
  fi

  # Phase 1: move to temp names to avoid overwrites
  while IFS='|' read -r src final; do
    tmp="$dir.tmp_$(basename "$src")"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "  mv -- '$src' '$tmp'"
    else
      mv -- "$src" "$tmp"
    fi
  done < "$mapfile"

  # Phase 2: move temp to final names
  while IFS='|' read -r src final; do
    tmp="$dir.tmp_$(basename "$src")"
    dest="$dir$final"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "  mv -- '$tmp' '$dest'"
    else
      mv -- "$tmp" "$dest"
    fi
  done < "$mapfile"

  rm -f "$mapfile"
done

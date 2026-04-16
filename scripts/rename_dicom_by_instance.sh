#!/opt/homebrew/bin/bash
set -euo pipefail
shopt -s nullglob

# Rename DICOM files inside each subdirectory using Instance Number (0020,0013)
# Usage:
#   ./rename_dicom_by_instance.sh /path/to/series_root
# If no arg is given, uses current directory.

root="${1:-.}"

# change to 1 for a dry-run (prints mv commands only)
DRY_RUN=0

# Helper: extract a single tag value from one DICOM file
get_tag() {
  local tag="$1"
  local file="$2"
  dcmdump +P "$tag" "$file" 2>/dev/null | awk -F'[][]' '{print $2}'
}

# Iterate over series directories (or just treat root as a series dir if no subdirs)
dirs=()
while IFS= read -r -d '' d; do dirs+=("$d"); done < <(find "$root" -mindepth 1 -maxdepth 1 -type d -print0)
if [ "${#dirs[@]}" -eq 0 ]; then
  dirs=("$root")
fi

for dir in "${dirs[@]}"; do
  # gather all files in this dir (non-recursive)
  files=()
  while IFS= read -r -d '' f; do files+=("$f"); done < <(find "$dir" -maxdepth 1 -type f -print0)

  [ "${#files[@]}" -eq 0 ] && continue

  echo "==> Processing: $dir"

  # Build rename plan: tmpname -> finalname
  declare -A tmp_for=()
  declare -A final_for=()
  declare -A seen_final=()

  for f in "${files[@]}"; do
    # Skip DICOMDIR or non-DICOM
    base="$(basename "$f")"
    [ "$base" = "DICOMDIR" ] && continue

    inst="$(get_tag "0020,0013" "$f")" || inst=""
    if [ -z "${inst:-}" ]; then
      echo "  Skipping (no Instance Number): $base"
      continue
    fi

    # Instance Number can have spaces; trim
    inst="$(echo "$inst" | tr -d '[:space:]')"

    # Final name
    final="$(printf "IMG%04d.dcm" "$inst")"

    # Detect collisions (two files with same instance number)
    if [ -n "${seen_final[$final]:-}" ]; then
      echo "  ERROR: collision on $final (multiple files share Instance Number $inst) in $dir"
      echo "         Use Option B (UID-based) or include SOPInstanceUID in name."
      exit 1
    fi
    seen_final["$final"]=1

    # Temporary name to avoid overwriting during renames
    tmp="$(printf ".tmp_%s_%s" "$inst" "$base")"

    tmp_for["$f"]="$tmp"
    final_for["$f"]="$final"
  done

  # Phase 1: move to temp names
  for f in "${!tmp_for[@]}"; do
    tmp_path="$dir/${tmp_for[$f]}"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "  mv -- '$f' '$tmp_path'"
    else
      mv -- "$f" "$tmp_path"
    fi
  done

  # Phase 2: move temp names to final names
  for f in "${!final_for[@]}"; do
    tmp_path="$dir/${tmp_for[$f]}"
    final_path="$dir/${final_for[$f]}"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "  mv -- '$tmp_path' '$final_path'"
    else
      mv -- "$tmp_path" "$final_path"
    fi
  done
done

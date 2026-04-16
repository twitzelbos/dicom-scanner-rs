#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

DRY_RUN=0

for d in SERIES*/; do
  # strip trailing slash for basename
  name="${d%/}"

  # extract digits after SERIES
  num="${name#SERIES}"

  # if it’s not digits, skip
  case "$num" in
    ''|*[!0-9]*)
      echo "Skipping (unexpected dir name): $d"
      continue
      ;;
  esac

  new="SER${num}"   # "SER" + digits => 7 or 8 chars for 4-5 digit nums

  if [ "$name" != "$new" ]; then
    echo "Renaming $name -> $new"
    if [ "$DRY_RUN" -eq 1 ]; then
      echo "mv -- '$name' '$new'"
    else
      mv -- "$name" "$new"
    fi
  fi
done

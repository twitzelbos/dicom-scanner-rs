#!/usr/bin/env bash

shopt -s nullglob

for dir in */; do
  # pick the first DICOM file in the folder
  first_dcm=$(find "$dir" -type f | head -n 1)
  [ -z "$first_dcm" ] && continue

  # extract Series Number (0020,0011)
  series_num=$(dcmdump +P "0020,0011" "$first_dcm" \
    | awk -F'[][]' '{print $2}')

  # skip if missing
  [ -z "$series_num" ] && continue

  # zero-pad to 4 digits
  newdir=$(printf "SERIES%04d" "$series_num")

  # avoid collisions
  if [ "$dir" != "$newdir/" ]; then
    echo "Renaming $dir → $newdir/"
    mv "$dir" "$newdir"
  fi
done

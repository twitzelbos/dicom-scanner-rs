#!/bin/bash

# sort_downloaded_zips.sh
# Sorts ZIP files containing DICOM data by extracting MRN and renaming files

# Check if correct number of arguments provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_path> <output_path>"
    echo "  input_path:  Directory to recursively search for ZIP files"
    echo "  output_path: Directory where sorted ZIP files will be copied"
    exit 1
fi

INPUT_PATH="$1"
OUTPUT_PATH="$2"

# Path to the dicom_scanner tool (adjust if needed)
DICOM_SCANNER="$(dirname "$0")/../target/debug/dicom_scanner"

# Check if dicom_scanner exists
if [ ! -f "$DICOM_SCANNER" ]; then
    echo "Error: dicom_scanner not found at $DICOM_SCANNER"
    echo "Please build the project first with 'cargo build'"
    exit 1
fi

# Check if input path exists
if [ ! -d "$INPUT_PATH" ]; then
    echo "Error: Input path does not exist: $INPUT_PATH"
    exit 1
fi

# Create output path if it doesn't exist
mkdir -p "$OUTPUT_PATH"

# Counter for processed files
PROCESSED=0
SKIPPED=0
ERRORS=0

echo "Searching for ZIP files in: $INPUT_PATH"
echo "Output directory: $OUTPUT_PATH"
echo "----------------------------------------"

# Use process substitution instead of pipe to avoid subshell
while IFS= read -r zipfile; do
    echo "Processing: $zipfile"

    # Extract MRN using dicom_scanner
    MRN=$("$DICOM_SCANNER" --mrn -f "$zipfile" 2>/dev/null | head -n 1)

    # Check if MRN was extracted
    if [ -z "$MRN" ]; then
        echo "  WARNING: No MRN found, skipping file"
        ((SKIPPED++))
        continue
    fi

    echo "  Found MRN: $MRN"

    # Clean MRN for filename (replace any problematic characters)
    # Replace forward slashes, backslashes, and other problematic chars with underscores
    CLEAN_MRN=$(echo "$MRN" | tr '/' '_' | tr '\\' '_' | tr ':' '_' | tr '*' '_' | tr '?' '_' | tr '"' '_' | tr '<' '_' | tr '>' '_' | tr '|' '_')

    # Determine output filename
    OUTPUT_FILE="$OUTPUT_PATH/${CLEAN_MRN}.zip"

    # Check if file already exists and find unique name
    if [ -f "$OUTPUT_FILE" ]; then
        COUNTER=1
        while [ -f "$OUTPUT_PATH/${CLEAN_MRN}_${COUNTER}.zip" ]; do
            ((COUNTER++))
        done
        OUTPUT_FILE="$OUTPUT_PATH/${CLEAN_MRN}_${COUNTER}.zip"
        echo "  File exists, using: $(basename "$OUTPUT_FILE")"
    fi

    # Copy the file
    cp "$zipfile" "$OUTPUT_FILE"

    if [ $? -eq 0 ]; then
        echo "  Copied to: $(basename "$OUTPUT_FILE")"
        ((PROCESSED++))
    else
        echo "  ERROR: Failed to copy file"
        ((ERRORS++))
    fi

    echo ""
done < <(find "$INPUT_PATH" -type f -name "*.zip")

echo "----------------------------------------"
echo "Processing complete!"
echo "Files processed: $PROCESSED"
echo "Files skipped: $SKIPPED"
echo "Errors: $ERRORS"
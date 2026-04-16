# DICOM Scanner (Rust)

A high-performance DICOM file scanner built in Rust that processes ZIP archives containing medical imaging files.

## Features

- Fast parallel DICOM file detection in ZIP archives
- Deep metadata extraction from DICOM files
- MR imaging-specific parameter extraction
- GE manufacturer private tag support
- Medical Record Number (MRN) extraction mode
- Batch processing capabilities
- Extract and organize DICOM files by Study/Series hierarchy

## Installation

```bash
# Clone the repository
git clone git@github.com:twitzelbos/dicom-scanner-rs.git
cd dicom-scanner-rs

# Build the project
cargo build --release

# (Optional) Set up pre-commit hooks
./scripts/setup-pre-commit.sh
```

## Usage

### Basic scanning
```bash
dicom_scanner --file archive.zip
```

### Extract only MRN
```bash
dicom_scanner --mrn --file archive.zip
```

### Extract and organize DICOM files by series
```bash
dicom_scanner --file archive.zip --output organized_dicoms
```

This will:
1. Extract all DICOM files from the ZIP archive
2. Create a directory structure organized by Study and Series with descriptive names:
   ```
   organized_dicoms/
   ├── Study_[StudyDescription]_[ShortUID]/
   │   ├── Series_[Number]_[SeriesDescription]_[ShortUID]/
   │   │   ├── IMG00001.dcm
   │   │   ├── IMG00002.dcm
   │   │   └── ...
   │   └── Series_[Number]_[SeriesDescription]_[ShortUID]/
   │       └── ...
   ├── series_metadata.csv
   ```
3. Preserve original filenames where possible
4. Use descriptive folder names from DICOM metadata (StudyDescription, SeriesDescription, SeriesNumber)
5. Include shortened UIDs for uniqueness
6. Generate a CSV file (`series_metadata.csv`) with metadata for all series, including:
   - **SeriesDescription**: Human-readable series description
   - **ProtocolName** (tag 0018,1030): Protocol name used for acquisition
   - **SeriesNumber**: Series number within the study
   - **Modality**: Imaging modality (MR, CT, etc.)
   - **SeriesInstanceUID**: Unique identifier for the series
   - **AcquisitionType**: 2D or 3D acquisition (for MR)
   - **PixelSpacing**: In-plane pixel spacing (mm)
   - **SliceThickness**: Slice thickness (mm)
   - **FOV**: Calculated field of view (mm x mm)
   - **TR**: Repetition Time (ms) for MR
   - **TE**: Echo Time (ms) for MR
   - **TI**: Inversion Time (ms) for MR
   - **FlipAngle**: Flip angle in degrees
   - **NumberOfAverages**: Number of signal averages (NEX)
   - **EchoTrainLength**: Echo train length for fast spin echo sequences
   - **ParallelImagingFactor**: Parallel imaging acceleration factor (GE ASSET, Siemens GRAPPA/iPAT, etc.)
   - **MagneticFieldStrength**: Scanner field strength in Tesla
   - **ImageType**: Image type (ORIGINAL\\PRIMARY, DERIVED\\SECONDARY, etc.)
   - **SpacingBetweenSlices**: Spacing between slices (mm)
   - **AcquisitionTime**: Time of acquisition
   - **AcquisitionDuration**: Total acquisition/scan duration in seconds (when available)
   - **DerivationDescription**: How the series was derived (if applicable)
   - **ReferencedSeriesUID**: Source series UID for derived series
   - **FileCount**: Number of DICOM files in the series

## Scripts

The `scripts/` folder contains utility scripts:

- `sort_downloaded_zips.sh`: Sort DICOM ZIP files by MRN

## License

MIT

## Development

This project uses pre-commit hooks to ensure code quality:

- **ShellCheck** & **shfmt**: For shell script linting and formatting
- **Clippy** & **rustfmt**: For Rust linting and formatting
- **cargo test**: Runs all tests before commit

### Setting up development environment

```bash
# Install pre-commit hooks
./scripts/setup-pre-commit.sh

# Run hooks manually
pre-commit run --all-files

# Skip hooks if needed (not recommended)
git commit --no-verify
```

## Contributing

Contributions welcome! Please feel free to submit a Pull Request.

All PRs must pass:
- Clippy linting
- Code formatting checks
- All tests
- Pre-commit hooks

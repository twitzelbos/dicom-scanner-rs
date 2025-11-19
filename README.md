# DICOM Scanner (Rust)

A high-performance DICOM file scanner built in Rust that processes ZIP archives containing medical imaging files.

## Features

- Fast parallel DICOM file detection in ZIP archives
- Deep metadata extraction from DICOM files
- MR imaging-specific parameter extraction
- GE manufacturer private tag support
- Medical Record Number (MRN) extraction mode
- Batch processing capabilities

## Installation

```bash
# Clone the repository
git clone git@github.com:twitzelbos/dicom-scanner-rs.git
cd dicom-scanner-rs

# Build the project
cargo build --release
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

## Scripts

The `scripts/` folder contains utility scripts:

- `sort_downloaded_zips.sh`: Sort DICOM ZIP files by MRN

## License

MIT

## Contributing

Contributions welcome! Please feel free to submit a Pull Request.
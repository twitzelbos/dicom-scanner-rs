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

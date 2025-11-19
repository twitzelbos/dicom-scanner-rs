# DICOM Scanner Project - Claude's Notes

## Project Overview
A high-performance Rust-based DICOM file scanner that processes ZIP archives containing medical imaging files. The project emphasizes parallel processing and efficient metadata extraction, particularly for MR (Magnetic Resonance) imaging data.

## Technical Stack
- **Language**: Rust (Edition 2024)
- **Build System**: Cargo
- **Key Dependencies**:
  - `dicom` (0.8.1) - DICOM file parsing
  - `zip` (0.6) - ZIP archive handling
  - `rayon` (1.10) - Parallel processing
  - `clap` (4.0) - CLI argument parsing

## Project Structure
```
dicom_scanner/
├── src/
│   └── main.rs              # Single source file (1,766 lines)
├── sample-archive/          # Test DICOM data (~500+ files)
│   ├── SER00001-SER00005/  # Series directories
│   ├── README.txt
│   ├── LOG.txt
│   └── RWC_dicom_list.csv
├── Cargo.toml              # Project configuration
├── Cargo.lock              # Dependency lock
├── .gitignore              # Ignores /target
├── DICOMDIR                # Raw DICOM directory file
└── *.zip                   # Test archives
```

## Core Functionality

### 1. Two-Phase Scanning Approach
- **Phase 1**: Fast DICOM detection using magic bytes ("DICM" at offset 128)
- **Phase 2**: Deep metadata extraction from validated DICOM files

### 2. Key Data Structures

#### DicomCandidate (Fast scan result)
- Basic file information: index, name, compressed/uncompressed size
- Used for quick DICOM file identification

#### DeepDicomCandidate (Deep scan result)
- Comprehensive metadata including:
  - Study/Series/SOP Instance UIDs
  - Manufacturer and modality information
  - MR-specific parameters (if applicable)
  - GE vendor-specific tags (if applicable)

#### DicomAcqMatrix
- Custom parser for DICOM acquisition matrix values
- Extracts non-zero value pairs from 4-value matrix

### 3. Main Components

#### Fast Scanner (`scan_dicom_candidates_parallel`)
- Parallel processing using Rayon
- Reads only first 132 bytes per file
- Validates DICOM magic bytes
- Returns list of potential DICOM files

#### Deep Scanner (`deep_scan_dicom_candidates_parallel`)
- Full DICOM parsing up to pixel data
- Extracts comprehensive metadata
- Special handling for:
  - MR modality images
  - GE manufacturer equipment
  - Enhanced MR images (SOP Class: 1.2.840.10008.5.1.4.1.1.4.1)

#### GE Private Tag Scanner (`scan_gems_parm_01`)
- Extracts 100+ GE-specific private DICOM tags
- Tags range: 0x0043,0x0010 to 0x0043,0x10BF
- Includes SAR, dB/dt, gradient parameters, coil information

### 4. MR Imaging Support
Extensive extraction of MR-specific parameters:
- Echo Time (TE), Repetition Time (TR)
- Flip Angle, Slice Thickness
- SAR (Specific Absorption Rate), dB/dt
- Receive Coil Name, Pixel Bandwidth
- Acquisition Matrix, Phase Encoding Direction
- Scanning Sequence, Sequence Variant
- B1 RMS values

## Command Line Interface
```bash
dicom_scanner --file <ZIP_PATH>
# or
dicom_scanner -f <ZIP_PATH>

# Output only MRN (Medical Record Number)
dicom_scanner --file <ZIP_PATH> --mrn
```

### Options
- `--file` or `-f`: Path to ZIP file containing DICOM files (required)
- `--mrn`: Output only the MRN (Patient ID) from the DICOM files

When using `--mrn`:
- Only unique patient IDs are output (one per line)
- All other output is suppressed
- Patient IDs with "N/A" values are filtered out

## Output Information
1. List of detected DICOM files
2. Performance metrics:
   - File read time
   - Detection time
   - Scan rate (files/second)
3. Compression statistics
4. Study/series organization
5. Detailed metadata for each DICOM file

## Performance Optimizations
- Parallel processing for header scanning
- Memory-efficient ZIP handling (load once, stream access)
- Selective DICOM parsing (stops at pixel data)
- Two-phase approach minimizes unnecessary parsing

## Test Data Available
- `sample-archive/`: 5 series with 500+ DICOM files
- Test archives: archive.zip, dicom.zip, female.zip, jb_prenuvo.zip
- Raw DICOMDIR file

## Current Limitations
1. Only processes ZIP archives (no direct file/directory support)
2. Single-file codebase (could benefit from modularization)
3. No unit tests
4. Limited documentation
5. No configuration file support
6. Output format is text-only (no JSON/CSV export)

## Potential Improvements
1. **Code Organization**: Split main.rs into modules
2. **Testing**: Add unit and integration tests
3. **Features**:
   - Direct DICOM file/directory scanning
   - Multiple output formats (JSON, CSV)
   - Configuration file support
   - Progress indicators for large archives
4. **Documentation**: Add inline documentation and examples
5. **Error Handling**: More descriptive error messages
6. **Performance**: Streaming ZIP processing for very large archives

## Recent Changes (2025-11-19)
1. **Added --mrn option**:
   - Extracts and outputs ONLY the Medical Record Number (Patient ID)
   - Absolutely no other output - no debug messages, no progress information
   - Filters out duplicate MRNs and "N/A" values
   - Perfect for piping to other commands or scripts
2. **Improved error handling**:
   - Fixed crash when DICOM acquisition matrix values are malformed
   - Now returns "N/A" instead of panicking on invalid matrix data

## Notes
- Project appears to be in active development
- Strong focus on GE MR imaging equipment
- Designed for batch processing of medical imaging archives
- Git repository initialized but no commits yet

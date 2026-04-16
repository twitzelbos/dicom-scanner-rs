# Siemens Enhanced MR: Derivation Reference Findings

## Problem

Siemens Enhanced MR DICOM files (SOP Class `1.2.840.10008.5.1.4.1.1.4.1`) from inline
reconstructions (ADC, CALC_BVAL, DIXON Water/Fat/FF, MPR, WB_COMP) do **not** encode
source-to-derived relationships in any DICOM tag — standard or private.

## Evidence (MAGNETOM Sola, syngo MR XA51)

Test dataset: 223 Enhanced MR files, 128 series. Examined DWI (Series 67, ORIGINAL)
and ADC (Series 68, DERIVED) in detail.

### Standard derivation tags

| Tag | Name | Status in ADC file |
|-----|------|--------------------|
| (0008,9124) | DerivationImageSequence | Present in all 48 per-frame functional groups, but **EMPTY (0 items)** |
| (0008,2112) | SourceImageSequence | **Not present** |
| (0008,9215) | DerivationCodeSequence | **Not present** |
| (0008,1115) | ReferencedSeriesSequence | **Not present** at top level |
| (0008,1250) | RelatedSeriesSequence | **Not present** |

### Siemens private tags

| Tag Group | Content | UIDs or cross-references? |
|-----------|---------|---------------------------|
| 0009 | Syngo Index Service (2 tags) | No |
| 0021 | MR SDR 01 — slice timing, diffusion params (3 top-level tags + per-frame sequences) | No |
| 0029 | CSA Image/Series headers | **Not present** in Enhanced MR objects |
| 0051 | Not present | N/A |

### Raw byte search

The DWI SeriesInstanceUID, SOPInstanceUID, and even the unique numeric fragment
`432165235` appear **nowhere** in the ADC file's raw bytes. The reverse is also true.

### ReferencedImageEvidenceSequence

Both DWI and ADC reference the **same two localizer series** (purpose code 121311 =
"Localizer"). These are planning references, not derivation sources.

## What DWI and ADC DO share

### Identical fields (strong linking signals)

| Field | Shared Value |
|-------|-------------|
| FrameOfReferenceUID | `1.3.12.2.1107.5.2.51.186540.2.20260325121421691.0.0.0` |
| ProtocolName | Identical |
| Geometry (Rows/Cols/Frames/SliceThickness/PixelSpacing/Orientation) | 208x256, 48 frames, 6mm, 1.797mm, axial |
| SharedFunctionalGroupsSequence | ALL sub-sequences identical (coils, TR/TE, GRAPPA, FOV, Siemens private block) |
| ReferencedImageEvidenceSequence | Same 2 localizer series |
| StudyInstanceUID | Same study |

### Differing fields (directional signals)

| Field | DWI (Series 67) | ADC (Series 68) |
|-------|-----------------|-----------------|
| ImageType | ORIGINAL/PRIMARY/DIFFUSION/NONE | DERIVED/PRIMARY/DIFFUSION/ADC |
| SeriesDescription | `<name>` | `<name>_ADC` |
| AcquisitionDateTime | 13:02:58 | 13:03:35 (+37s) |
| InstanceCreationTime | 13:03:01 | 13:03:38 (+37s) |
| SeriesTime | 13:03:43.010 | 13:03:43.024 (+14ms) |

### SeriesInstanceUID timestamp encoding

Siemens UIDs follow the structure:

```
1.3.12.2.1107.5.2.51.<DeviceSerial>.<YYYYMMDDHHMMSS><Counter>.0.0.0
```

Segment 10 (0-indexed 9) embeds a timestamp. DWI and ADC share the same
14-character timestamp prefix `20260325130224` with differing counter suffixes.

### Timestamp grouping reliability

| Scenario | Result |
|----------|--------|
| STIR + WB_COMP + WB_COMP_FILT | Correctly grouped (same timestamp) |
| DIXON opp/in/fat/water outputs | **Fragmented** across adjacent seconds (e.g., 13:17:59 vs 13:18:00) |
| Unrelated secondary captures | **Collide** on same timestamp (false positives) |

The timestamp is a useful hint but requires a tolerance window (~2 seconds) and
additional signals to be reliable.

## Recommended heuristic for automated derivation resolution

When DerivationImageSequence is empty and no explicit source references exist,
match DERIVED series to ORIGINAL series using multi-factor scoring:

1. **FrameOfReferenceUID** — must match (same spatial coordinate system)
2. **UID timestamp** — within +/- 2 seconds (same acquisition pipeline)
3. **SeriesDescription** — ORIGINAL description is a prefix/substring of DERIVED description
4. **ImageType** — ORIGINAL on source, DERIVED on target
5. **Geometry** — matching Rows, Columns, NumberOfFrames, SliceThickness, PixelSpacing, ImageOrientation
6. **Temporal order** — DERIVED AcquisitionDateTime is after ORIGINAL

No single factor is sufficient. A combination of FrameOfReferenceUID + description
pattern + ImageType direction provides the strongest signal.

## Scope

These findings are specific to:
- **Scanner**: Siemens MAGNETOM Sola
- **Software**: syngo MR XA51
- **SOP Class**: Enhanced MR Image Storage (1.2.840.10008.5.1.4.1.1.4.1)
- **Reconstruction types observed**: ADC, CALC_BVAL, DIXON (Water/Fat/InPhase/OppPhase/FF), WB_COMP, WB_COMP_FILT, MPR

Other Siemens software versions or scanner models may behave differently.

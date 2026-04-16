use std::{
    io::{Cursor, Read},
    path::{Path, PathBuf},
    str::FromStr,
    sync::Arc,
    time::{Duration, Instant},
};

use dicom::core::VR;
use dicom::core::value::Value;
use dicom::object::mem::InMemDicomObject;
use dicom::{
    core::Tag,
    dictionary_std::tags::{self, PATIENT_ID},
};

use dicom::object::StandardDataDictionary;
use rayon::prelude::*;
use zip::ZipArchive;

use clap::Parser;

use dicom::object::OpenFileOptions;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to ZIP file containing DICOM files
    #[arg(short, long)]
    file: PathBuf,

    /// Output only the MRN (Medical Record Number) of the study
    #[arg(long)]
    mrn: bool,

    /// Output directory to extract and organize DICOM files by series
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output directory to extract XProtocol data from Siemens DICOM files
    #[arg(long)]
    xprot: Option<PathBuf>,

    /// Analyze derivation relationships between DICOM series
    #[arg(long)]
    derivations: bool,
}

#[derive(Debug, Clone)]
pub struct DicomCandidate {
    pub index: usize,
    pub name: String,
    pub compressed_size: u64,
    pub uncompressed_size: u64,
}

pub struct DeepDicomCandidate {
    pub index: usize,
    pub name: String,
    pub compressed_size: u64,
    pub uncompressed_size: u64,
    pub study_instance_uid: String,
    pub series_instance_uid: String,
    pub sop_instance_uid: String,
    pub manufacturer: String,
    pub modality: String,
    pub patient_id: String,
    pub study_description: String,
    pub series_description: String,
    pub series_number: String,
    pub protocol_name: String,
    pub acquisition_type: String, // 2D or 3D
    pub pixel_spacing: String,
    pub slice_thickness: String,
    pub acquisition_time: String,
    pub rows: String,
    pub columns: String,
    pub repetition_time: String, // TR
    pub echo_time: String,       // TE
    pub inversion_time: String,  // TI
    pub derivation_description: String,
    pub referenced_series_uid: String, // For tracking derived series relationships
    pub acquisition_duration: String,  // Protocol duration in seconds
    // MR-specific technical parameters
    pub flip_angle: String,
    pub number_of_averages: String, // NEX
    pub echo_train_length: String,
    pub parallel_imaging_factor: String,
    pub magnetic_field_strength: String,
    // Image set information
    pub spacing_between_slices: String,
    pub image_type: String,
}

#[derive(Debug)]
struct DicomAcqMatrix {
    values: [u32; 4],
}

impl DicomAcqMatrix {
    /// Returns the two non-zero values, based on the valid DICOM acquisition matrix pattern.
    pub fn extract_pair(&self) -> (u32, u32) {
        let [a, b, c, d] = self.values;
        if a == 0 && d == 0 {
            (b, c)
        } else if b == 0 && c == 0 {
            (a, d)
        } else {
            unreachable!(
                "Struct invariant violated: DicomAcqMatrix is only constructed when valid"
            );
        }
    }
}

impl std::str::FromStr for DicomAcqMatrix {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split('\\').collect();

        if parts.len() != 4 {
            return Err(format!("Expected 4 values, got {}", parts.len()));
        }

        let mut values = [0u32; 4];
        for (i, part) in parts.iter().enumerate() {
            values[i] = part
                .trim()
                .parse::<u32>()
                .map_err(|e| format!("Failed to parse '{}' as u32: {}", part, e))?;
        }

        let [a, b, c, d] = values;
        if !((a == 0 && d == 0) || (b == 0 && c == 0)) {
            return Err(format!(
                "Invalid acquisition matrix: must have (1st and 4th == 0) or (2nd and 3rd == 0), got: {:?}",
                values
            ));
        }

        Ok(DicomAcqMatrix { values })
    }
}

// Calculate the acquistion resolution based on DICOM tags as best as we can
pub fn calculate_acq_resolution(
    acq_mtx: String,
    dcm_rows: String,
    dcm_cols: String,
    pixel_spacing: String,
) -> String {
    // println!("acq_mtx: {}, pixel_spacing {}", acq_mtx, pixel_spacing);

    // Handle cases where acquisition matrix is invalid
    let acq_matrix = match DicomAcqMatrix::from_str(&acq_mtx) {
        Ok(matrix) => matrix,
        Err(_) => {
            // If parsing fails, return N/A
            return "N/A".to_string();
        }
    };
    let (rows, cols) = acq_matrix.extract_pair();

    let pixel_spacing = pixel_spacing
        .split('\\')
        .map(|s| s.parse::<f32>().unwrap_or(0.0))
        .collect::<Vec<f32>>();

    let dcm_rows = dcm_rows.parse::<f32>().unwrap_or(0.0);

    let dcm_cols = dcm_cols.parse::<f32>().unwrap_or(0.0);

    let fov_x = pixel_spacing[0] * dcm_rows;
    let fov_y = pixel_spacing[1] * dcm_cols;

    let resolution_x = fov_x / rows as f32;
    let resolution_y: f32 = fov_y / cols as f32;

    format!("{} x {} mm", resolution_x, resolution_y)
}

/*
fn get_tag_string<D: DicomObject>(dcm: &D, tag: Tag) -> String {
    dcm.element(tag)
        .ok()
        .and_then(|e| e.to_str().ok())
        .map(|s| s.to_string())
        .unwrap_or_else(|| "N/A".to_string())
}

fn get_tag_float<D: DicomObject>(dcm: &D, tag: Tag) -> f32 {
    dcm.element(tag)
        .ok()
        .and_then(|e| e.value().to_float32().ok())
        .unwrap_or(f32::NAN)
}


fn get_element_str(obj: &impl dicom::object::DicomObject, tag: Tag) -> Option<String> {
    obj.element(tag)
        .ok()?
        .value()
        .to_str()
        .ok()
        .map(|s| s.to_string())
}
*/

pub fn scan_gems_parm_01(
    dcm_object: &InMemDicomObject<StandardDataDictionary>,
    suppress_output: bool,
) {
    // VR: LO
    let gehc_private_creator_id = dcm_object
        .element(Tag(0x0043, 0x0010))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _bitmap_of_prescan_options = dcm_object
        .element(Tag(0x0043, 0x1001))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _gradient_offset_x = dcm_object
        .element(Tag(0x0043, 0x1002))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _gradient_offset_y = dcm_object
        .element(Tag(0x0043, 0x1003))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _gradient_offset_z = dcm_object
        .element(Tag(0x0043, 0x1004))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS, no longer used in DV26
    let _image_is_original = dcm_object
        .element(Tag(0x0043, 0x1005))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _number_of_epi_shots = dcm_object
        .element(Tag(0x0043, 0x1006))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _views_per_segment = dcm_object
        .element(Tag(0x0043, 0x1007))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _respiratory_rate_bpm = dcm_object
        .element(Tag(0x0043, 0x1008))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _respiratory_trigger_point = dcm_object
        .element(Tag(0x0043, 0x1009))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _type_of_receiver_used = dcm_object
        .element(Tag(0x0043, 0x100A))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let peak_dbdt = dcm_object
        .element(Tag(0x0043, 0x100B))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let dbdt_limits_percent = dcm_object
        .element(Tag(0x0043, 0x100C))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let psd_estimatated_limit = dcm_object
        .element(Tag(0x0043, 0x100D))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let psd_estimated_limit_tps = dcm_object
        .element(Tag(0x0043, 0x100E))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS, no longer used in DV26
    let sar_avg_head = dcm_object
        .element(Tag(0x0043, 0x100F))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US, no longer used in DV26
    let _window_value = dcm_object
        .element(Tag(0x0043, 0x1010))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _ge_image_integrity = dcm_object
        .element(Tag(0x0043, 0x101C))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS, no longer used in DV26
    let _level_value = dcm_object
        .element(Tag(0x0043, 0x101D))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB, no longer used in DV26
    let _unique_image_identifier = dcm_object
        .element(Tag(0x0043, 0x1028))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB
    let _histogram_tables = dcm_object
        .element(Tag(0x0043, 0x1029))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB
    let _user_defined_data = dcm_object
        .element(Tag(0x0043, 0x102A))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS[4], no longer used in DV26
    let _private_scan_options = dcm_object
        .element(Tag(0x0043, 0x102B))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _effective_echo_spacing = dcm_object
        .element(Tag(0x0043, 0x102C))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SH
    // (String slop field 1 in legacy GE MR images)
    let _filter_mode = dcm_object
        .element(Tag(0x0043, 0x102D))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SH
    let _string_slop_field_2 = dcm_object
        .element(Tag(0x0043, 0x102E))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS (image_type)
    let _raw_data_type = dcm_object
        .element(Tag(0x0043, 0x102F))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _vas_collapse_flag = dcm_object
        .element(Tag(0x0043, 0x1030))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[2], not used in DV26
    let _ra_coord_of_target_recon_center = dcm_object
        .element(Tag(0x0043, 0x1031))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _vas_flags = dcm_object
        .element(Tag(0x0043, 0x1032))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: FL
    let _neg_scanspacing = dcm_object
        .element(Tag(0x0043, 0x1033))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: IS
    let _offset_frequency = dcm_object
        .element(Tag(0x0043, 0x1034))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UL
    let _user_usage_tag = dcm_object
        .element(Tag(0x0043, 0x1035))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UL
    let _user_fill_map_msw = dcm_object
        .element(Tag(0x0043, 0x1036))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UL
    let _user_fill_map_lsw = dcm_object
        .element(Tag(0x0043, 0x1037))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: FL[24]
    let _user_data25_48 = dcm_object
        .element(Tag(0x0043, 0x1038))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: IS[4]
    let _slop_int_6_9 = dcm_object
        .element(Tag(0x0043, 0x1039))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: IS[8]
    let _slop_int_10_17 = dcm_object
        .element(Tag(0x0043, 0x1060))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SH, not used in DV26
    let _scanner_study_id = dcm_object
        .element(Tag(0x0043, 0x1062))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VS: DS[3 or 4]
    // 3 on single gradient coil systems, on multiple gradient coil systems the 4th value is the selected gradient coil
    let _scanner_table_entry = dcm_object
        .element(Tag(0x0043, 0x106f))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: ST
    let _paradigm_description = dcm_object
        .element(Tag(0x0043, 0x1071))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UI
    let _paradigm_uid = dcm_object
        .element(Tag(0x0043, 0x1072))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US
    let _experiment_type = dcm_object
        .element(Tag(0x0043, 0x1073))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US
    let _number_of_rest_volumes = dcm_object
        .element(Tag(0x0043, 0x1074))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US
    let _number_of_active_volumes = dcm_object
        .element(Tag(0x0043, 0x1075))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US
    let _number_of_dummy_scans = dcm_object
        .element(Tag(0x0043, 0x1076))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SH
    let _application_name = dcm_object
        .element(Tag(0x0043, 0x1077))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SH
    let _application_version = dcm_object
        .element(Tag(0x0043, 0x1078))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US
    let _slices_per_volume = dcm_object
        .element(Tag(0x0043, 0x1079))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US
    let _expected_time_points = dcm_object
        .element(Tag(0x0043, 0x107A))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: FL[1...n]
    let _regressor_values = dcm_object
        .element(Tag(0x0043, 0x107B))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: FL
    let _delay_after_slice_group = dcm_object
        .element(Tag(0x0043, 0x107C))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: US
    let _recon_mode_flag_word = dcm_object
        .element(Tag(0x0043, 0x107D))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1...n]
    let _pacc_specific_information = dcm_object
        .element(Tag(0x0043, 0x107E))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1...n]
    let _private_data = dcm_object
        .element(Tag(0x0043, 0x107F))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1...n]
    let _coil_id_data = dcm_object
        .element(Tag(0x0043, 0x1080))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _ge_coil_name = dcm_object
        .element(Tag(0x0043, 0x1081))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1...n]
    let _system_configuration_information = dcm_object
        .element(Tag(0x0043, 0x1082))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[2]
    let _asset_r_factors = dcm_object
        .element(Tag(0x0043, 0x1083))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[5]
    let _additional_asset_data = dcm_object
        .element(Tag(0x0043, 0x1084))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UT
    let _debug_data_text = dcm_object
        .element(Tag(0x0043, 0x1085))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB
    let _debug_data_bin = dcm_object
        .element(Tag(0x0043, 0x1086))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UT
    let _software_version_long = dcm_object
        .element(Tag(0x0043, 0x1087))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UI
    let _pure_cal_series_uid = dcm_object
        .element(Tag(0x0043, 0x1088))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[3]
    let _gov_body_dbdt_sar_def = dcm_object
        .element(Tag(0x0043, 0x1089))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: CS
    let _private_inplace_pe_dir = dcm_object
        .element(Tag(0x0043, 0x108A))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB, not used in DV26
    let _fmri_binary_data_block = dcm_object
        .element(Tag(0x0043, 0x108B))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[6]
    let _voxel_location = dcm_object
        .element(Tag(0x0043, 0x108C))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    //VR: DS[7n]
    let _sat_band_locations = dcm_object
        .element(Tag(0x0043, 0x108D))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[3]
    let _spectro_prescan_values = dcm_object
        .element(Tag(0x0043, 0x108E))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[3]
    let _spectro_parameters = dcm_object
        .element(Tag(0x0043, 0x108F))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1..n]
    let _sar_definition = dcm_object
        .element(Tag(0x0043, 0x1090))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _sar_value = dcm_object
        .element(Tag(0x0043, 0x1091))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _image_error_text = dcm_object
        .element(Tag(0x0043, 0x1092))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _spectro_quantitation_values = dcm_object
        .element(Tag(0x0043, 0x1093))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _spectro_ratio_values = dcm_object
        .element(Tag(0x0043, 0x1094))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _prescan_reuse_string = dcm_object
        .element(Tag(0x0043, 0x1095))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: CS
    let _content_qualification = dcm_object
        .element(Tag(0x0043, 0x1096))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[8]
    let _image_filtering_parameters = dcm_object
        .element(Tag(0x0043, 0x1097))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: UI
    let _asset_acquisition_calibration_uid = dcm_object
        .element(Tag(0x0043, 0x1098))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1..n]
    let _extended_options = dcm_object
        .element(Tag(0x0043, 0x1099))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: IS
    let _rx_stack_identification = dcm_object
        .element(Tag(0x0043, 0x109A))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let _npw_factor = dcm_object
        .element(Tag(0x0043, 0x109B))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB
    let _research_tag_1 = dcm_object
        .element(Tag(0x0043, 0x109C))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB
    let _research_tag_2 = dcm_object
        .element(Tag(0x0043, 0x109D))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB
    let _research_tag_3 = dcm_object
        .element(Tag(0x0043, 0x109E))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: OB
    let _research_tag_4 = dcm_object
        .element(Tag(0x0043, 0x109F))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SQ
    let _spectroscopy_pixel_sequence = dcm_object
        .element(Tag(0x0043, 0x10A0))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SQ
    let _spectroscopy_default_display_sequence = dcm_object
        .element(Tag(0x0043, 0x10A1))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VS: DS[1..n]
    let _mef_data = dcm_object
        .element(Tag(0x0043, 0x10A2))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: CS
    let _asl_contrast_technique = dcm_object
        .element(Tag(0x0043, 0x10A3))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _detailed_text_for_asl_labeling = dcm_object
        .element(Tag(0x0043, 0x10A4))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: IS
    let _duration_of_label_or_ctrl_pulse = dcm_object
        .element(Tag(0x0043, 0x10A5))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS, not used in DV26
    let _offset_frequency_fastb1map = dcm_object
        .element(Tag(0x0043, 0x10A6))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let _motion_encoding_factor = dcm_object
        .element(Tag(0x0043, 0x10A7))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[3]
    let _dual_drive_mode_amplitude_attenuation_phase_offset = dcm_object
        .element(Tag(0x0043, 0x10A8))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1..n]
    let _threed_cal_data = dcm_object
        .element(Tag(0x0043, 0x10A9))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1..n]
    let _additional_filtering_parameters = dcm_object
        .element(Tag(0x0043, 0x10AA))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _silenz_data = dcm_object
        .element(Tag(0x0043, 0x10AB))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1..n], reserved for future use
    let _qmap_delay_data = dcm_object
        .element(Tag(0x0043, 0x10AC))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _other_recovery_times_values = dcm_object
        .element(Tag(0x0043, 0x10AD))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[1..n]
    let _other_recovery_times_labels = dcm_object
        .element(Tag(0x0043, 0x10AE))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _additional_echo_times = dcm_object
        .element(Tag(0x0043, 0x10AF))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: FL
    let _rescan_time_in_acquisition = dcm_object
        .element(Tag(0x0043, 0x10B0))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _excitation_mode = dcm_object
        .element(Tag(0x0043, 0x10B1))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _advanced_eddy_correction = dcm_object
        .element(Tag(0x0043, 0x10B3))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: SS
    let _mrf_transmit_gain = dcm_object
        .element(Tag(0x0043, 0x10B4))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _mr_table_position_information = dcm_object
        .element(Tag(0x0043, 0x10B2))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[7]
    let _multiband_parameters = dcm_object
        .element(Tag(0x0043, 0x10B6))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO[4]
    let _compressed_sensing_parameters = dcm_object
        .element(Tag(0x0043, 0x10B7))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let _grad_comp_parameters = dcm_object
        .element(Tag(0x0043, 0x10B8))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _parallel_transmit_information = dcm_object
        .element(Tag(0x0043, 0x10B9))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS
    let _echo_spacing = dcm_object
        .element(Tag(0x0043, 0x10BA))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _pixel_information = dcm_object
        .element(Tag(0x0043, 0x10BB))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: IS
    let _heart_beats_pattern = dcm_object
        .element(Tag(0x0043, 0x10BC))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: LO
    let _hyper_kat_factor = dcm_object
        .element(Tag(0x0043, 0x10BD))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    // VR: DS[1..n]
    let _delta_transmit_gain = dcm_object
        .element(Tag(0x0043, 0x10BF))
        .map_or("N/A".to_string(), |e| {
            e.value().to_str().unwrap().to_string()
        });

    if !suppress_output {
        println!(
            "GEHC Private Creator ID: {} Peak dB/dt: {} dB/dt limits: {}% PSD estimated limit: {} Tps: {} SAR avg head: {}",
            gehc_private_creator_id,
            peak_dbdt,
            dbdt_limits_percent,
            psd_estimatated_limit,
            psd_estimated_limit_tps,
            sar_avg_head
        );
    }
}

// --- Siemens XProtocol extraction and recursive SQ traversal ---

pub fn is_csa_header(data: &[u8]) -> bool {
    data.len() >= 4 && &data[0..4] == b"SV10"
}

pub fn parse_csa_header(data: &[u8]) -> Option<String> {
    if !is_csa_header(data) {
        return None;
    }

    let mut offset = 8; // skip SV10 signature (4 bytes) + unused (4 bytes)

    if data.len() < offset + 4 {
        return None;
    }

    let num_elements = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?) as usize;
    offset += 4;

    // skip unused 4 bytes after num_elements
    offset += 4;

    for _ in 0..num_elements {
        if data.len() < offset + 64 + 4 + 4 + 4 + 4 {
            break;
        }

        let name_bytes = &data[offset..offset + 64];
        let name = String::from_utf8_lossy(name_bytes)
            .trim_end_matches('\0')
            .to_string();
        offset += 64;

        // skip VM (4 bytes), VR (4 bytes), sync pattern (4 bytes)
        offset += 12;

        let num_items = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?) as usize;
        offset += 4;

        let mut element_data = Vec::new();
        for _ in 0..num_items {
            if data.len() < offset + 4 {
                break;
            }
            let item_len = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?) as usize;
            offset += 4;

            // skip 3 unused u32 fields
            if data.len() < offset + 12 {
                break;
            }
            offset += 12;

            if item_len > 0 && offset + item_len <= data.len() {
                element_data.extend_from_slice(&data[offset..offset + item_len]);
            }

            // items are padded to 4-byte boundaries
            let padded_len = (item_len + 3) & !3;
            if offset + padded_len > data.len() {
                break;
            }
            offset += padded_len;
        }

        if name.contains("Protocol") || name == "MrPhoenixProtocol" {
            if let Some(start_pos) = element_data
                .windows(10)
                .position(|w| w.starts_with(b"### ASCCONV") || w.starts_with(b"<XProtocol>"))
            {
                let content = &element_data[start_pos..];
                let trimmed_len = content
                    .iter()
                    .rposition(|&b| b != 0)
                    .map(|i| i + 1)
                    .unwrap_or(0);
                return Some(String::from_utf8_lossy(&content[..trimmed_len]).to_string());
            } else if !element_data.is_empty() {
                let trimmed_len = element_data
                    .iter()
                    .rposition(|&b| b != 0)
                    .map(|i| i + 1)
                    .unwrap_or(0);
                if trimmed_len > 0 {
                    return Some(String::from_utf8_lossy(&element_data[..trimmed_len]).to_string());
                }
            }
        }
    }

    None
}

pub fn parse_dicom_sequence(data: &[u8]) -> Option<String> {
    let mut offset = 0;

    if data.len() < 8 || data[0] != 0xFE || data[1] != 0xFF || data[2] != 0x00 || data[3] != 0xE0 {
        return None;
    }
    offset += 4;

    let item_length = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?) as usize;
    offset += 4;

    let item_end = if item_length == 0xFFFFFFFF {
        data.len()
    } else {
        offset + item_length
    };

    while offset + 8 <= item_end && offset + 8 <= data.len() {
        let group = u16::from_le_bytes(data[offset..offset + 2].try_into().ok()?);
        let element = u16::from_le_bytes(data[offset + 2..offset + 4].try_into().ok()?);
        offset += 4;

        if group == 0xFFFE && (element == 0xE00D || element == 0xE0DD) {
            offset += 4;
            continue;
        }

        let vr_bytes = if offset + 2 <= data.len() {
            [data[offset], data[offset + 1]]
        } else {
            break;
        };
        let vr = String::from_utf8_lossy(&vr_bytes).to_string();
        offset += 2;

        let valid_vr = vr_bytes[0].is_ascii_uppercase() && vr_bytes[1].is_ascii_uppercase();

        let value_length = if !valid_vr {
            offset -= 2;
            if offset + 4 > data.len() {
                break;
            }
            let len = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?) as usize;
            offset += 4;
            len
        } else if vr == "OB" || vr == "OW" || vr == "SQ" || vr == "UN" || vr == "UT" {
            offset += 2;
            if offset + 4 > data.len() {
                break;
            }
            u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?) as usize
        } else {
            if offset + 2 > data.len() {
                break;
            }
            u16::from_le_bytes(data[offset..offset + 2].try_into().ok()?) as usize
        };

        if !valid_vr {
            // offset already advanced above
        } else if vr == "OB" || vr == "OW" || vr == "SQ" || vr == "UN" || vr == "UT" {
            offset += 4;
        } else {
            offset += 2;
        }

        if offset + value_length <= data.len() {
            let tag_data = &data[offset..offset + value_length];

            if let Some(xprot_pos) = tag_data
                .windows(11)
                .position(|window| window.starts_with(b"<XProtocol>"))
            {
                let content = &tag_data[xprot_pos..];
                let trimmed_len = content
                    .iter()
                    .rposition(|&b| b != 0)
                    .map(|i| i + 1)
                    .unwrap_or(0);
                return Some(String::from_utf8_lossy(&content[..trimmed_len]).to_string());
            }
        }

        offset += value_length;
    }

    None
}

pub fn extract_tag_recursive(
    obj: &InMemDicomObject<StandardDataDictionary>,
    tag_to_find: Tag,
) -> Option<String> {
    if let Ok(tag_element) = obj.element(tag_to_find)
        && let Ok(bytes) = tag_element.to_bytes()
    {
        if tag_to_find == Tag(0x0021, 0x10fe) {
            if bytes.len() >= 4
                && bytes[0] == 0xFE
                && bytes[1] == 0xFF
                && bytes[2] == 0x00
                && bytes[3] == 0xE0
            {
                if let Some(parsed) = parse_dicom_sequence(&bytes) {
                    return Some(parsed);
                }
            } else if is_csa_header(&bytes)
                && let Some(parsed) = parse_csa_header(&bytes)
            {
                return Some(parsed);
            }
        }
        let trimmed_len = bytes
            .iter()
            .rposition(|&b| b != 0)
            .map(|i| i + 1)
            .unwrap_or(0);
        return Some(String::from_utf8_lossy(&bytes[..trimmed_len]).to_string());
    }

    for element in obj.iter() {
        if element.vr() == VR::SQ
            && let Value::Sequence(sequence) = element.value()
        {
            for item in sequence.items() {
                if let Some(value) = extract_tag_recursive(item, tag_to_find) {
                    return Some(value);
                }
            }
        }
    }

    None
}

pub fn get_tag_string(obj: &InMemDicomObject<StandardDataDictionary>, tag: Tag) -> String {
    if let Ok(element) = obj.element(tag)
        && let Ok(s) = element.value().to_str()
    {
        return s.to_string();
    }

    for element in obj.iter() {
        if element.vr() == VR::SQ
            && let Value::Sequence(sequence) = element.value()
        {
            for item in sequence.items() {
                let result = get_tag_string(item, tag);
                if result != "N/A" {
                    return result;
                }
            }
        }
    }

    "N/A".to_string()
}

#[derive(Debug, Clone)]
struct DerivationCode {
    code_value: String,
    coding_scheme: String,
    code_meaning: String,
}

#[derive(Debug, Clone)]
struct SeriesDerivationInfo {
    series_instance_uid: String,
    #[allow(dead_code)]
    series_description: String,
    series_number: String,
    modality: String,
    image_type: String,
    is_derived: bool,
    file_count: usize,
    derivation_description: String,
    derivation_codes: Vec<DerivationCode>,
    frame_of_reference_uid: String,
    referenced_series_uids: Vec<String>,
    source_sop_uids: Vec<String>,
    referenced_image_sop_uids: Vec<String>,
}

#[derive(Debug, Clone)]
struct DerivationEdge {
    source_series_uid: String,
    evidence_tags: Vec<&'static str>,
}

fn extract_all_from_sq(
    obj: &InMemDicomObject<StandardDataDictionary>,
    sq_tag: Tag,
    inner_tag: Tag,
) -> Vec<String> {
    let mut results = Vec::new();
    if let Ok(sq_element) = obj.element(sq_tag)
        && let Value::Sequence(sequence) = sq_element.value()
    {
        for item in sequence.items() {
            if let Ok(inner_element) = item.element(inner_tag)
                && let Ok(s) = inner_element.value().to_str()
            {
                let trimmed = s.trim().to_string();
                if !trimmed.is_empty() {
                    results.push(trimmed);
                }
            }
        }
    }
    results
}

fn extract_derivation_codes(obj: &InMemDicomObject<StandardDataDictionary>) -> Vec<DerivationCode> {
    let mut codes = Vec::new();
    if let Ok(sq_element) = obj.element(tags::DERIVATION_CODE_SEQUENCE)
        && let Value::Sequence(sequence) = sq_element.value()
    {
        for item in sequence.items() {
            let code_value = item
                .element(tags::CODE_VALUE)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap_or_default().trim().to_string()
                });
            let coding_scheme = item
                .element(tags::CODING_SCHEME_DESIGNATOR)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap_or_default().trim().to_string()
                });
            let code_meaning = item
                .element(tags::CODE_MEANING)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap_or_default().trim().to_string()
                });
            codes.push(DerivationCode {
                code_value,
                coding_scheme,
                code_meaning,
            });
        }
    }
    codes
}

fn format_series_label(
    uid: &str,
    descriptions: &std::collections::HashMap<String, String>,
    numbers: &std::collections::HashMap<String, String>,
) -> String {
    let desc = descriptions.get(uid).map(|s| s.as_str()).unwrap_or("N/A");
    let num = numbers.get(uid).map(|s| s.as_str()).unwrap_or("?");
    let short_uid = uid.split('.').next_back().unwrap_or(uid);
    format!(
        "Series {} \"{}\" [...{}]",
        num,
        desc,
        &short_uid[..8.min(short_uid.len())]
    )
}

/// Parse DICOM objects from either a ZIP archive or a directory of DICOM files.
/// Returns a Vec of parsed DICOM objects.
fn load_dicom_objects(
    input_path: &Path,
) -> Result<Vec<InMemDicomObject<StandardDataDictionary>>, Box<dyn std::error::Error>> {
    let mut objects = Vec::new();

    if input_path.is_dir() {
        fn visit_dir(
            dir: &Path,
            objects: &mut Vec<InMemDicomObject<StandardDataDictionary>>,
        ) -> Result<(), Box<dyn std::error::Error>> {
            for entry in std::fs::read_dir(dir)? {
                let entry = entry?;
                let path = entry.path();
                if path.is_dir() {
                    visit_dir(&path, objects)?;
                } else if let Ok(obj) = OpenFileOptions::new()
                    .read_until(tags::PIXEL_DATA)
                    .open_file(&path)
                {
                    objects.push(obj.into_inner());
                }
            }
            Ok(())
        }
        visit_dir(input_path, &mut objects)?;
    } else {
        // Try as ZIP first, fall back to single DICOM file
        let mut file_data = Vec::new();
        std::fs::File::open(input_path)?.read_to_end(&mut file_data)?;

        if let Ok(mut archive) = ZipArchive::new(Cursor::new(&file_data)) {
            for i in 0..archive.len() {
                let file = match archive.by_index(i) {
                    Ok(f) => f,
                    Err(_) => continue,
                };
                if file.size() < 132 || file.is_dir() {
                    continue;
                }
                if let Ok(obj) = OpenFileOptions::new()
                    .read_until(tags::PIXEL_DATA)
                    .from_reader(file)
                {
                    objects.push(obj.into_inner());
                }
            }
        } else if let Ok(obj) = OpenFileOptions::new()
            .read_until(tags::PIXEL_DATA)
            .open_file(input_path)
        {
            objects.push(obj.into_inner());
        }
    }

    Ok(objects)
}

pub fn analyze_derivations(input_path: &Path) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::{HashMap, HashSet};

    println!("Loading DICOM files from {}...", input_path.display());
    let objects = load_dicom_objects(input_path)?;
    println!("Parsed {} DICOM files", objects.len());

    // Pass 1: Build SOP→Series map and collect representative per series
    let mut sop_to_series: HashMap<String, String> = HashMap::new();
    let mut series_file_counts: HashMap<String, usize> = HashMap::new();
    let mut series_representatives: HashMap<String, usize> = HashMap::new(); // series_uid → index in objects
    let mut series_descriptions: HashMap<String, String> = HashMap::new();
    let mut series_numbers: HashMap<String, String> = HashMap::new();

    for (idx, dcm_object) in objects.iter().enumerate() {
        let sop_uid = get_tag_string(dcm_object, tags::SOP_INSTANCE_UID);
        let series_uid = get_tag_string(dcm_object, tags::SERIES_INSTANCE_UID);

        if sop_uid != "N/A" && series_uid != "N/A" {
            sop_to_series.insert(sop_uid, series_uid.clone());
        }

        *series_file_counts.entry(series_uid.clone()).or_insert(0) += 1;

        if let std::collections::hash_map::Entry::Vacant(e) =
            series_representatives.entry(series_uid.clone())
        {
            series_descriptions.insert(
                series_uid.clone(),
                get_tag_string(dcm_object, tags::SERIES_DESCRIPTION),
            );
            series_numbers.insert(series_uid, get_tag_string(dcm_object, tags::SERIES_NUMBER));
            e.insert(idx);
        }
    }

    println!(
        "Found {} series across {} files\n",
        series_representatives.len(),
        objects.len()
    );

    // Pass 2: Extract derivation info per series
    let mut derivation_infos: Vec<SeriesDerivationInfo> = Vec::new();

    for (series_uid, &obj_idx) in &series_representatives {
        let dcm_object = &objects[obj_idx];

        let image_type = get_tag_string(dcm_object, tags::IMAGE_TYPE);
        let is_derived = image_type
            .split('\\')
            .next()
            .map(|s| s.trim().eq_ignore_ascii_case("DERIVED"))
            .unwrap_or(false);

        let derivation_description = get_tag_string(dcm_object, tags::DERIVATION_DESCRIPTION);
        let derivation_codes = extract_derivation_codes(dcm_object);
        let frame_of_reference_uid = get_tag_string(dcm_object, tags::FRAME_OF_REFERENCE_UID);
        let modality = get_tag_string(dcm_object, tags::MODALITY);

        let referenced_series_uids = extract_all_from_sq(
            dcm_object,
            tags::REFERENCED_SERIES_SEQUENCE,
            tags::SERIES_INSTANCE_UID,
        );

        let source_sop_uids = extract_all_from_sq(
            dcm_object,
            tags::SOURCE_IMAGE_SEQUENCE,
            tags::REFERENCED_SOP_INSTANCE_UID,
        );

        let referenced_image_sop_uids = extract_all_from_sq(
            dcm_object,
            tags::REFERENCED_IMAGE_SEQUENCE,
            tags::REFERENCED_SOP_INSTANCE_UID,
        );

        derivation_infos.push(SeriesDerivationInfo {
            series_instance_uid: series_uid.clone(),
            series_description: series_descriptions
                .get(series_uid)
                .cloned()
                .unwrap_or_default(),
            series_number: series_numbers.get(series_uid).cloned().unwrap_or_default(),
            modality,
            image_type,
            is_derived,
            file_count: *series_file_counts.get(series_uid).unwrap_or(&0),
            derivation_description,
            derivation_codes,
            frame_of_reference_uid,
            referenced_series_uids,
            source_sop_uids,
            referenced_image_sop_uids,
        });
    }

    // Sort by series number
    derivation_infos
        .sort_by_key(|info| info.series_number.trim().parse::<i32>().unwrap_or(i32::MAX));

    // Phase 3: Build derivation graph
    let mut derivation_graph: HashMap<String, Vec<DerivationEdge>> = HashMap::new();
    let mut unresolved_sops: HashMap<String, Vec<String>> = HashMap::new();

    for info in &derivation_infos {
        let mut edges: Vec<DerivationEdge> = Vec::new();

        // Direct series references
        for ref_series_uid in &info.referenced_series_uids {
            if ref_series_uid != &info.series_instance_uid {
                edges.push(DerivationEdge {
                    source_series_uid: ref_series_uid.clone(),
                    evidence_tags: vec!["ReferencedSeriesSequence (0008,1115)"],
                });
            }
        }

        // SOP → Series for SourceImageSequence
        let mut source_series: HashSet<String> = HashSet::new();
        let mut unresolved: Vec<String> = Vec::new();
        for sop_uid in &info.source_sop_uids {
            if let Some(series_uid) = sop_to_series.get(sop_uid) {
                if series_uid != &info.series_instance_uid {
                    source_series.insert(series_uid.clone());
                }
            } else {
                unresolved.push(sop_uid.clone());
            }
        }
        for series_uid in source_series {
            if let Some(existing) = edges.iter_mut().find(|e| e.source_series_uid == series_uid) {
                if !existing
                    .evidence_tags
                    .contains(&"SourceImageSequence (0008,2112)")
                {
                    existing
                        .evidence_tags
                        .push("SourceImageSequence (0008,2112)");
                }
            } else {
                edges.push(DerivationEdge {
                    source_series_uid: series_uid,
                    evidence_tags: vec!["SourceImageSequence (0008,2112)"],
                });
            }
        }

        // SOP → Series for ReferencedImageSequence
        let mut ref_img_series: HashSet<String> = HashSet::new();
        for sop_uid in &info.referenced_image_sop_uids {
            if let Some(series_uid) = sop_to_series.get(sop_uid) {
                if series_uid != &info.series_instance_uid {
                    ref_img_series.insert(series_uid.clone());
                }
            } else {
                unresolved.push(sop_uid.clone());
            }
        }
        for series_uid in ref_img_series {
            if let Some(existing) = edges.iter_mut().find(|e| e.source_series_uid == series_uid) {
                if !existing
                    .evidence_tags
                    .contains(&"ReferencedImageSequence (0008,1140)")
                {
                    existing
                        .evidence_tags
                        .push("ReferencedImageSequence (0008,1140)");
                }
            } else {
                edges.push(DerivationEdge {
                    source_series_uid: series_uid,
                    evidence_tags: vec!["ReferencedImageSequence (0008,1140)"],
                });
            }
        }

        if !unresolved.is_empty() {
            unresolved_sops.insert(info.series_instance_uid.clone(), unresolved);
        }

        if !edges.is_empty() {
            derivation_graph.insert(info.series_instance_uid.clone(), edges);
        }
    }

    // Phase 4: Output report

    println!("=== DICOM Derivation Analysis ===\n");

    // Series overview
    println!("--- Series Overview ---");
    for info in &derivation_infos {
        let tag = if info.is_derived { "D" } else { "O" };
        let label = format_series_label(
            &info.series_instance_uid,
            &series_descriptions,
            &series_numbers,
        );
        println!(
            "  [{}] {}  ({} files, {})",
            tag, label, info.file_count, info.modality
        );
        println!("      ImageType: {}", info.image_type);
        if info.derivation_description != "N/A" && !info.derivation_description.is_empty() {
            println!(
                "      DerivationDescription (0008,2111): {}",
                info.derivation_description
            );
        }
        for code in &info.derivation_codes {
            println!(
                "      DerivationCode (0008,9215): {} ({}) \"{}\"",
                code.code_value, code.coding_scheme, code.code_meaning
            );
        }
        if !info.source_sop_uids.is_empty() {
            println!(
                "      SourceImageSequence (0008,2112): {} SOP references",
                info.source_sop_uids.len()
            );
        }
        if !info.referenced_image_sop_uids.is_empty() {
            println!(
                "      ReferencedImageSequence (0008,1140): {} SOP references",
                info.referenced_image_sop_uids.len()
            );
        }
        if !info.referenced_series_uids.is_empty() {
            println!(
                "      ReferencedSeriesSequence (0008,1115): {:?}",
                info.referenced_series_uids
            );
        }
        if info.frame_of_reference_uid != "N/A" {
            let short_for = info
                .frame_of_reference_uid
                .split('.')
                .next_back()
                .unwrap_or(&info.frame_of_reference_uid);
            println!("      FrameOfReference: ...{}", short_for);
        }
        println!();
    }

    // Derivation graph (derived → sources)
    let has_any_derived = derivation_infos.iter().any(|i| i.is_derived);
    if has_any_derived {
        println!("--- Derivation Map (derived <- source) ---");
        for info in &derivation_infos {
            if !info.is_derived {
                continue;
            }
            let label = format_series_label(
                &info.series_instance_uid,
                &series_descriptions,
                &series_numbers,
            );
            println!("  {}", label);

            if let Some(edges) = derivation_graph.get(&info.series_instance_uid) {
                for edge in edges {
                    let source_label = format_series_label(
                        &edge.source_series_uid,
                        &series_descriptions,
                        &series_numbers,
                    );
                    println!("    <- {}", source_label);
                    println!("       via: {}", edge.evidence_tags.join(", "));
                }
            }

            if let Some(unresolved) = unresolved_sops.get(&info.series_instance_uid) {
                println!(
                    "    <- {} unresolved external SOP reference(s) not in this dataset",
                    unresolved.len()
                );
                println!(
                    "       via: SourceImageSequence (0008,2112) / ReferencedImageSequence (0008,1140)"
                );
            }

            let total_refs = info.source_sop_uids.len() + info.referenced_image_sop_uids.len();
            if total_refs == 0
                && info.referenced_series_uids.is_empty()
                && !unresolved_sops.contains_key(&info.series_instance_uid)
            {
                println!("    (no explicit derivation references found in DICOM tags)");
                println!("    ImageType indicates DERIVED — likely inline reconstruction");
            }

            println!();
        }
    }

    // Dependency tree (source → derived)
    println!("--- Dependency Tree (source -> derived) ---");
    let mut source_to_derived: HashMap<String, Vec<(String, Vec<&'static str>)>> = HashMap::new();
    for (derived_uid, edges) in &derivation_graph {
        for edge in edges {
            source_to_derived
                .entry(edge.source_series_uid.clone())
                .or_default()
                .push((derived_uid.clone(), edge.evidence_tags.clone()));
        }
    }

    for info in &derivation_infos {
        let has_sources = derivation_graph.contains_key(&info.series_instance_uid);
        let has_derivatives = source_to_derived.contains_key(&info.series_instance_uid);

        if has_derivatives {
            let label = format_series_label(
                &info.series_instance_uid,
                &series_descriptions,
                &series_numbers,
            );
            let tag = if info.is_derived {
                "DERIVED"
            } else {
                "ORIGINAL"
            };
            println!("  {} [{}]", label, tag);
            if let Some(derivatives) = source_to_derived.get(&info.series_instance_uid) {
                let mut sorted_derivs = derivatives.clone();
                sorted_derivs.sort_by_key(|(uid, _)| {
                    series_numbers
                        .get(uid)
                        .and_then(|n| n.trim().parse::<i32>().ok())
                        .unwrap_or(i32::MAX)
                });
                for (derived_uid, evidence) in &sorted_derivs {
                    let derived_label =
                        format_series_label(derived_uid, &series_descriptions, &series_numbers);
                    println!("    -> {} via {}", derived_label, evidence.join(", "));
                }
            }
        } else if !has_sources {
            let label = format_series_label(
                &info.series_instance_uid,
                &series_descriptions,
                &series_numbers,
            );
            let tag = if info.is_derived {
                "DERIVED"
            } else {
                "ORIGINAL"
            };
            println!("  {} [{}] (no derivation links)", label, tag);
        }
    }

    // Frame of Reference groups
    println!("\n--- Frame of Reference Groups ---");
    let mut for_groups: HashMap<String, Vec<String>> = HashMap::new();
    for info in &derivation_infos {
        if info.frame_of_reference_uid != "N/A" {
            for_groups
                .entry(info.frame_of_reference_uid.clone())
                .or_default()
                .push(info.series_instance_uid.clone());
        }
    }
    for (for_uid, series_uids) in &for_groups {
        let short_for = for_uid.split('.').next_back().unwrap_or(for_uid);
        let labels: Vec<String> = series_uids
            .iter()
            .map(|uid| {
                let num = series_numbers.get(uid).map(|s| s.as_str()).unwrap_or("?");
                format!("S{}", num.trim())
            })
            .collect();
        println!("  FoR ...{}: {}", short_for, labels.join(", "));
    }

    println!();
    Ok(())
}

pub fn extract_xprotocol_from_zip(
    zip_bytes: &[u8],
    output_dir: &Path,
) -> Result<usize, Box<dyn std::error::Error>> {
    use std::collections::HashSet;
    use std::fs;
    use std::io::Write;

    fs::create_dir_all(output_dir)?;

    let mut archive = ZipArchive::new(Cursor::new(zip_bytes))?;
    let mut seen_series: HashSet<String> = HashSet::new();
    let mut extracted_count = 0;

    for i in 0..archive.len() {
        let mut file = archive.by_index(i)?;
        if file.is_dir() {
            continue;
        }

        let mut buf = Vec::new();
        if file.read_to_end(&mut buf).is_err() {
            continue;
        }

        let dcm_object = match OpenFileOptions::new()
            .read_until(tags::PIXEL_DATA)
            .from_reader(Cursor::new(&buf))
        {
            Ok(obj) => obj,
            Err(_) => continue,
        };

        let manufacturer = dcm_object
            .element(tags::MANUFACTURER)
            .map_or(String::new(), |e| {
                e.value().to_str().unwrap_or_default().to_uppercase()
            });

        if !manufacturer.contains("SIEMENS") {
            continue;
        }

        let series_uid = dcm_object
            .element(tags::SERIES_INSTANCE_UID)
            .map_or("unknown".to_string(), |e| {
                e.value().to_str().unwrap_or_default().to_string()
            });

        if seen_series.contains(&series_uid) {
            continue;
        }

        let tag_1019 = extract_tag_recursive(&dcm_object, Tag(0x0021, 0x1019));
        let tag_10fe = extract_tag_recursive(&dcm_object, Tag(0x0021, 0x10fe));

        if tag_1019.is_none() && tag_10fe.is_none() {
            seen_series.insert(series_uid);
            continue;
        }

        let series_desc =
            dcm_object
                .element(tags::SERIES_DESCRIPTION)
                .map_or("Unknown".to_string(), |e| {
                    let s = e.value().to_str().unwrap_or_default().trim().to_string();
                    if s.is_empty() {
                        "Unknown".to_string()
                    } else {
                        s
                    }
                });

        let series_num = dcm_object
            .element(tags::SERIES_NUMBER)
            .map_or(String::new(), |e| {
                e.value().to_str().unwrap_or_default().trim().to_string()
            });

        let short_uid = &series_uid.split('.').next_back().unwrap_or(&series_uid)[..8.min(
            series_uid
                .split('.')
                .next_back()
                .unwrap_or(&series_uid)
                .len(),
        )];

        let filename = if series_num.is_empty() {
            format!("{}_{}.xprot", sanitize_filename(&series_desc), short_uid)
        } else {
            format!(
                "{:04}_{}_{}.xprot",
                series_num.parse::<i32>().unwrap_or(0),
                sanitize_filename(&series_desc),
                short_uid
            )
        };

        let output_path = output_dir.join(&filename);
        let mut out_file = fs::File::create(&output_path)?;

        if let Some(value) = &tag_1019 {
            write!(out_file, "{}", value)?;
        }
        if let Some(value) = &tag_10fe {
            if tag_1019.is_some() {
                writeln!(out_file)?;
                writeln!(out_file)?;
            }
            write!(out_file, "{}", value)?;
        }

        println!(
            "  Extracted XProtocol: {} (series: {})",
            filename, series_desc
        );
        extracted_count += 1;
        seen_series.insert(series_uid);
    }

    Ok(extracted_count)
}

pub fn deep_scan_dicom_candidates_parallel(
    zip_bytes: &[u8],
    suppress_output: bool,
) -> Result<Vec<DeepDicomCandidate>, Box<dyn std::error::Error>> {
    let mut archive = ZipArchive::new(Cursor::new(zip_bytes))?;
    let mut all_candidates = Vec::new();

    for i in 0..archive.len() {
        let file = match archive.by_index(i) {
            Ok(f) => f,
            Err(_) => continue,
        };

        if file.size() < 132 {
            continue;
        }

        let name = file.name().to_string();
        let compressed_size = file.compressed_size();
        let uncompressed_size = file.size();

        // If this succeeds, we have a DICOM file, I suppose
        let dcm_result = OpenFileOptions::new()
            .read_until(tags::PIXEL_DATA)
            .from_reader(file);

        if dcm_result.is_err() {
            continue;
        }

        let dcm_object = dcm_result.unwrap();

        let study_instance_uid = get_tag_string(&dcm_object, tags::STUDY_INSTANCE_UID);
        let series_instance_uid = get_tag_string(&dcm_object, tags::SERIES_INSTANCE_UID);
        let patient_id = get_tag_string(&dcm_object, PATIENT_ID);
        let sop_instance_uid = get_tag_string(&dcm_object, tags::SOP_INSTANCE_UID);
        let modality = get_tag_string(&dcm_object, tags::MODALITY);

        let protocol_name = get_tag_string(&dcm_object, tags::PROTOCOL_NAME);
        let study_description = get_tag_string(&dcm_object, tags::STUDY_DESCRIPTION);
        let series_description = get_tag_string(&dcm_object, tags::SERIES_DESCRIPTION);
        let _series_date = get_tag_string(&dcm_object, tags::SERIES_DATE);
        let series_number = get_tag_string(&dcm_object, tags::SERIES_NUMBER);
        let _series_time = get_tag_string(&dcm_object, tags::SERIES_TIME);
        let manufacturer = get_tag_string(&dcm_object, tags::MANUFACTURER);

        let acquisition_type = get_tag_string(&dcm_object, tags::MR_ACQUISITION_TYPE);
        let acquisition_time = get_tag_string(&dcm_object, tags::ACQUISITION_TIME);
        let pixel_spacing = get_tag_string(&dcm_object, tags::PIXEL_SPACING);
        let slice_thickness = get_tag_string(&dcm_object, tags::SLICE_THICKNESS);
        let rows = get_tag_string(&dcm_object, tags::ROWS);
        let columns = get_tag_string(&dcm_object, tags::COLUMNS);
        let repetition_time = get_tag_string(&dcm_object, tags::REPETITION_TIME);
        let echo_time = get_tag_string(&dcm_object, tags::ECHO_TIME);
        let inversion_time = get_tag_string(&dcm_object, tags::INVERSION_TIME);
        let derivation_description = get_tag_string(&dcm_object, tags::DERIVATION_DESCRIPTION);

        // Check for referenced series (source series for derived images)
        let referenced_series_uid =
            if let Ok(ref_series_seq) = dcm_object.element(tags::REFERENCED_SERIES_SEQUENCE) {
                // Get the first item in the sequence
                if let dicom::core::value::Value::Sequence(seq) = ref_series_seq.value() {
                    if let Some(first_item) = seq.items().first() {
                        // Get the Series Instance UID from the referenced series
                        first_item
                            .element(tags::SERIES_INSTANCE_UID)
                            .map_or("N/A".to_string(), |e| {
                                e.value().to_str().unwrap().to_string()
                            })
                    } else {
                        "N/A".to_string()
                    }
                } else {
                    "N/A".to_string()
                }
            } else {
                "N/A".to_string()
            };

        let mut acquisition_duration = get_tag_string(&dcm_object, Tag(0x0018, 0x9073));
        let flip_angle = get_tag_string(&dcm_object, tags::FLIP_ANGLE);
        let number_of_averages = get_tag_string(&dcm_object, tags::NUMBER_OF_AVERAGES);
        let echo_train_length = get_tag_string(&dcm_object, tags::ECHO_TRAIN_LENGTH);

        // Parallel imaging factor - try vendor-specific tags in priority order
        let mut parallel_imaging_factor = {
            let val = get_tag_string(&dcm_object, Tag(0x0018, 0x9069));
            if val != "N/A" {
                val
            } else {
                let val = get_tag_string(&dcm_object, Tag(0x0051, 0x1011));
                if val != "N/A" {
                    val
                } else {
                    get_tag_string(&dcm_object, Tag(0x0018, 0x9078))
                }
            }
        };

        let magnetic_field_strength = get_tag_string(&dcm_object, tags::MAGNETIC_FIELD_STRENGTH);
        let spacing_between_slices = get_tag_string(&dcm_object, tags::SPACING_BETWEEN_SLICES);
        let image_type = get_tag_string(&dcm_object, tags::IMAGE_TYPE);
        let sop_class_uid = get_tag_string(&dcm_object, tags::SOP_CLASS_UID);

        if !suppress_output {
            println!("sop_class_uid: {}", sop_class_uid);
        }

        if sop_class_uid == *"1.2.840.10008.5.1.4.1.1.4.1" && modality == *"MR" {
            if !suppress_output {
                println!("This is an enhanced MR image DICOM file");
            }

            let _acquisition_number = get_tag_string(&dcm_object, tags::ACQUISITION_NUMBER);
            let _acquisiton_date_time = get_tag_string(&dcm_object, tags::ACQUISITION_DATE_TIME);
            let _content_qualification = get_tag_string(&dcm_object, tags::CONTENT_QUALIFICATION);
            let _resonant_nucleus = get_tag_string(&dcm_object, tags::RESONANT_NUCLEUS);
            let _kspace_filtering = get_tag_string(&dcm_object, tags::K_SPACE_FILTERING);
            let _magnetic_field_strength =
                get_tag_string(&dcm_object, tags::MAGNETIC_FIELD_STRENGTH);
            let _applicable_safety_standard_agency =
                get_tag_string(&dcm_object, tags::APPLICABLE_SAFETY_STANDARD_AGENCY);
            let _applicable_safety_standard_description =
                get_tag_string(&dcm_object, tags::APPLICABLE_SAFETY_STANDARD_DESCRIPTION);
            let _image_comments = get_tag_string(&dcm_object, tags::IMAGE_COMMENTS);
            let _isocenter_position = get_tag_string(&dcm_object, tags::ISOCENTER_POSITION);
            let _b1rms = get_tag_string(&dcm_object, tags::B1RMS);
            let _acquisition_contrast = get_tag_string(&dcm_object, tags::ACQUISITION_CONTRAST);
            let mr_fov_geometry_sequence =
                get_tag_string(&dcm_object, tags::MRFOV_GEOMETRY_SEQUENCE);
            let _inplane_phase_encoding_direction =
                get_tag_string(&dcm_object, tags::IN_PLANE_PHASE_ENCODING_DIRECTION);
            let mr_acquisition_frequency_encoding_steps =
                get_tag_string(&dcm_object, tags::MR_ACQUISITION_FREQUENCY_ENCODING_STEPS);
            let mr_acquisition_phase_encoding_steps_inplane = get_tag_string(
                &dcm_object,
                tags::MR_ACQUISITION_PHASE_ENCODING_STEPS_IN_PLANE,
            );
            let mr_acquisition_phase_encoding_steps_outofplane = get_tag_string(
                &dcm_object,
                tags::MR_ACQUISITION_PHASE_ENCODING_STEPS_OUT_OF_PLANE,
            );
            let _percent_sampling = get_tag_string(&dcm_object, tags::PERCENT_SAMPLING);
            let _percent_phase_field_of_view =
                get_tag_string(&dcm_object, tags::PERCENT_PHASE_FIELD_OF_VIEW);

            if !suppress_output {
                println!(
                    "MRFOV_GEOMETRY_SEQUENCE: {} freq: {} phas: {} kz: {}",
                    mr_fov_geometry_sequence,
                    mr_acquisition_frequency_encoding_steps,
                    mr_acquisition_phase_encoding_steps_inplane,
                    mr_acquisition_phase_encoding_steps_outofplane
                );
            }

            continue;
        }

        // If the Modality is "MR", get some additional information
        if modality == "MR" {
            let te = get_tag_string(&dcm_object, tags::ECHO_TIME);
            let tr = get_tag_string(&dcm_object, tags::REPETITION_TIME);
            let sar = get_tag_string(&dcm_object, tags::SAR);
            let _db_dt = get_tag_string(&dcm_object, tags::D_BDT);
            let _isocenter_position = get_tag_string(&dcm_object, tags::ISOCENTER_POSITION);
            let receive_coil_name = get_tag_string(&dcm_object, tags::RECEIVE_COIL_NAME);
            let pixel_bandwidth = get_tag_string(&dcm_object, tags::PIXEL_BANDWIDTH);
            let _number_pe = get_tag_string(&dcm_object, tags::NUMBER_OF_PHASE_ENCODING_STEPS);
            let acq_matrix = get_tag_string(&dcm_object, tags::ACQUISITION_MATRIX);
            let phase_encoding_direction =
                get_tag_string(&dcm_object, tags::IN_PLANE_PHASE_ENCODING_DIRECTION);
            let reconstruction_diameter =
                get_tag_string(&dcm_object, tags::RECONSTRUCTION_DIAMETER);
            let pixel_spacing = get_tag_string(&dcm_object, tags::PIXEL_SPACING);
            let rows = get_tag_string(&dcm_object, tags::ROWS);
            let columns = get_tag_string(&dcm_object, tags::COLUMNS);

            let _b1_rms = get_tag_string(&dcm_object, tags::B1RMS);
            let _bits_allocated = get_tag_string(&dcm_object, tags::BITS_ALLOCATED);
            let _bits_stored = get_tag_string(&dcm_object, tags::BITS_STORED);
            let _high_bit = get_tag_string(&dcm_object, tags::HIGH_BIT);
            let scanning_sequence = get_tag_string(&dcm_object, tags::SCANNING_SEQUENCE);
            let sequence_variant = get_tag_string(&dcm_object, tags::SEQUENCE_VARIANT);
            let scan_options = get_tag_string(&dcm_object, tags::SCAN_OPTIONS);
            let mr_acquisition_type = get_tag_string(&dcm_object, tags::MR_ACQUISITION_TYPE);
            let _inversion_time = get_tag_string(&dcm_object, tags::INVERSION_TIME);

            let _sequence_name = get_tag_string(&dcm_object, tags::SEQUENCE_NAME);
            let center_to_center_slice_gap =
                get_tag_string(&dcm_object, tags::SPACING_BETWEEN_SLICES);
            let percent_sampling = get_tag_string(&dcm_object, tags::PERCENT_SAMPLING);
            let percent_phase_fov = get_tag_string(&dcm_object, tags::PERCENT_PHASE_FIELD_OF_VIEW);
            let flip_angle = get_tag_string(&dcm_object, tags::FLIP_ANGLE);

            let _variable_flip_flag = get_tag_string(&dcm_object, tags::VARIABLE_FLIP_ANGLE_FLAG);
            let slice_thickness = get_tag_string(&dcm_object, tags::SLICE_THICKNESS);

            // Seems like GE does not report:
            // - dbdt
            // - isocenter position
            // - number_pe
            // - sequence_name

            // image matrix

            if !suppress_output {
                println!(
                    "{} \"{}\" [{},{},{}] DIM: {}, SAR: {} RX Coil {} BW: {}Hz/px, TE: {}, TR: {}, FA: {}, AMTX: {} PE_dir: {} FOV: {} pFOV: {}%, samp: {}%, RES: {}, rows: {}, cols: {}, thick: {}, c2c: {}, res: {}",
                    series_number,
                    series_description,
                    scanning_sequence,
                    sequence_variant,
                    scan_options,
                    mr_acquisition_type,
                    sar,
                    receive_coil_name,
                    pixel_bandwidth,
                    te,
                    tr,
                    flip_angle,
                    acq_matrix,
                    phase_encoding_direction,
                    reconstruction_diameter,
                    percent_phase_fov,
                    percent_sampling,
                    pixel_spacing,
                    rows,
                    columns,
                    slice_thickness,
                    center_to_center_slice_gap,
                    calculate_acq_resolution(
                        acq_matrix.clone(),
                        rows.clone(),
                        columns.clone(),
                        pixel_spacing.clone(),
                    )
                );
            }

            if manufacturer == "GE MEDICAL SYSTEMS" {
                let internal_sequence_name = get_tag_string(&dcm_object, Tag(0x0019, 0x109E));

                // this tag is "FL" as VR (single float)
                let ge_acquisition_duration = dcm_object
                    .element(Tag(0x0019, 0x105A))
                    .map_or(f32::NAN, |e| e.value().to_float32().unwrap());

                // If we didn't get acquisition duration from standard tag, use GE-specific one
                if acquisition_duration == "N/A" && !ge_acquisition_duration.is_nan() {
                    acquisition_duration = format!("{:.2}", ge_acquisition_duration);
                }

                let number_of_echoes = get_tag_string(&dcm_object, Tag(0x0019, 0x107E));
                let _table_delta = get_tag_string(&dcm_object, Tag(0x0019, 0x107F));
                let _gehc_private_creator_id = get_tag_string(&dcm_object, Tag(0x0043, 0x0010));
                let _bitmap_of_prescan_options = get_tag_string(&dcm_object, Tag(0x0043, 0x1001));
                let _gradient_offset_x = get_tag_string(&dcm_object, Tag(0x0043, 0x1002));
                let _gradient_offset_y = get_tag_string(&dcm_object, Tag(0x0043, 0x1003));
                let _gradient_offset_z = get_tag_string(&dcm_object, Tag(0x0043, 0x1004));

                let _image_is_original = get_tag_string(&dcm_object, Tag(0x0043, 0x1005));
                let _number_of_epi_shots = get_tag_string(&dcm_object, Tag(0x0043, 0x1006));
                let _views_per_segment = get_tag_string(&dcm_object, Tag(0x0043, 0x1007));
                let _respiratory_rate_bpm = get_tag_string(&dcm_object, Tag(0x0043, 0x1008));
                let _respiratory_trigger_point = get_tag_string(&dcm_object, Tag(0x0043, 0x1009));
                let _type_of_receiver_used = get_tag_string(&dcm_object, Tag(0x0043, 0x100A));
                let _peak_dbdt = get_tag_string(&dcm_object, Tag(0x0043, 0x100B));
                let _dbdt_limits_percent = get_tag_string(&dcm_object, Tag(0x0043, 0x100C));
                let _psd_estimatated_limit = get_tag_string(&dcm_object, Tag(0x0043, 0x100D));
                let _psd_estimated_limit_tps = get_tag_string(&dcm_object, Tag(0x0043, 0x100E));
                let _sar_avg_head = get_tag_string(&dcm_object, Tag(0x0043, 0x100F));
                let _application_name = get_tag_string(&dcm_object, Tag(0x0043, 0x1077));
                let _application_version = get_tag_string(&dcm_object, Tag(0x0043, 0x1078));
                let _slices_per_volume = get_tag_string(&dcm_object, Tag(0x0043, 0x1079));
                let asset_r_factors = get_tag_string(&dcm_object, Tag(0x0043, 0x1083));

                // Use GE Asset R factors for parallel imaging factor if not already set
                if parallel_imaging_factor == "N/A" && asset_r_factors != "N/A" {
                    parallel_imaging_factor = asset_r_factors.clone();
                }

                // note the acquisition duration is in micro seconds
                if !suppress_output {
                    println!(
                        "{} {:#?} {} {}",
                        internal_sequence_name,
                        Duration::from_micros(ge_acquisition_duration as u64),
                        number_of_echoes,
                        asset_r_factors,
                    );
                }

                //scan_gems_parm_01(&dcm_object, suppress_output);
            }
        }
        all_candidates.push(DeepDicomCandidate {
            index: i,
            name,
            compressed_size,
            uncompressed_size,
            study_instance_uid,
            series_instance_uid,
            sop_instance_uid,
            modality,
            manufacturer,
            patient_id,
            study_description,
            series_description,
            series_number,
            protocol_name,
            acquisition_type,
            pixel_spacing,
            slice_thickness,
            acquisition_time,
            rows,
            columns,
            repetition_time,
            echo_time,
            inversion_time,
            derivation_description,
            referenced_series_uid,
            acquisition_duration,
            flip_angle,
            number_of_averages,
            echo_train_length,
            parallel_imaging_factor,
            magnetic_field_strength,
            spacing_between_slices,
            image_type,
        });
    }
    Ok(all_candidates)
}

pub fn scan_dicom_candidates_parallel(
    zip_bytes: &[u8],
) -> Result<Vec<DicomCandidate>, Box<dyn std::error::Error>> {
    let mut archive = ZipArchive::new(Cursor::new(zip_bytes))?;
    let mut all_candidates = Vec::new();

    for i in 0..archive.len() {
        let mut file = match archive.by_index(i) {
            Ok(f) => f,
            Err(_) => continue,
        };

        if file.size() < 132 {
            continue;
        }

        // get the first 132 bytes of the dcm_object
        let mut header = [0u8; 132];
        if file.read_exact(&mut header).is_ok() {
            all_candidates.push((
                i,
                file.name().to_string(),
                file.compressed_size(),
                file.size(),
                header,
            ));
        }
    }

    let results: Vec<_> = all_candidates
        .into_par_iter()
        .filter_map(|(i, name, compressed, uncompressed, header)| {
            if &header[128..132] == b"DICM" {
                Some(DicomCandidate {
                    index: i,
                    name,
                    compressed_size: compressed,
                    uncompressed_size: uncompressed,
                })
            } else {
                None
            }
        })
        .collect();

    Ok(results)
}

/// Export series metadata to CSV file
fn export_series_metadata_csv(
    deep_candidates: &[DeepDicomCandidate],
    output_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::Write;

    // Group files by series to get one row per series
    let mut series_map: std::collections::HashMap<(String, String), &DeepDicomCandidate> =
        std::collections::HashMap::new();

    // For each series, we'll use the first file's metadata as representative
    for candidate in deep_candidates {
        let key = (
            candidate.study_instance_uid.clone(),
            candidate.series_instance_uid.clone(),
        );
        series_map.entry(key).or_insert(candidate);
    }

    // Create CSV file
    let csv_path = output_dir.join("series_metadata.csv");
    let mut csv_file = File::create(&csv_path)?;

    // Write CSV header
    writeln!(
        csv_file,
        "SeriesDescription,ProtocolName,SeriesNumber,Modality,SeriesInstanceUID,AcquisitionType,PixelSpacing,SliceThickness,SpacingBetweenSlices,FOV,TR,TE,TI,FlipAngle,NumberOfAverages,EchoTrainLength,ParallelImagingFactor,MagneticFieldStrength,ImageType,AcquisitionTime,AcquisitionDuration,DerivationDescription,ReferencedSeriesUID,FileCount"
    )?;

    // Count files per series
    let mut series_file_count: std::collections::HashMap<(String, String), usize> =
        std::collections::HashMap::new();
    for candidate in deep_candidates {
        let key = (
            candidate.study_instance_uid.clone(),
            candidate.series_instance_uid.clone(),
        );
        *series_file_count.entry(key).or_insert(0) += 1;
    }

    // Sort series by study and series for consistent output
    let mut sorted_series: Vec<_> = series_map.into_iter().collect();
    sorted_series.sort_by(|a, b| match a.0.0.cmp(&b.0.0) {
        std::cmp::Ordering::Equal => a.0.1.cmp(&b.0.1),
        other => other,
    });

    // Write data rows
    for ((study_uid, series_uid), candidate) in sorted_series {
        let file_count = series_file_count
            .get(&(study_uid.clone(), series_uid.clone()))
            .unwrap_or(&0);

        // Escape fields that might contain commas or quotes
        let escape_csv_field = |field: &str| -> String {
            if field.contains(',') || field.contains('"') || field.contains('\n') {
                format!("\"{}\"", field.replace("\"", "\"\""))
            } else {
                field.to_string()
            }
        };

        // Calculate FOV if pixel spacing and matrix size are available
        let fov = if candidate.pixel_spacing != "N/A"
            && candidate.rows != "N/A"
            && candidate.columns != "N/A"
        {
            // Parse pixel spacing (format: "value1\\value2")
            let spacing_parts: Vec<&str> = candidate.pixel_spacing.split('\\').collect();
            if spacing_parts.len() >= 2 {
                if let (Ok(spacing_x), Ok(spacing_y), Ok(rows), Ok(cols)) = (
                    spacing_parts[0].parse::<f32>(),
                    spacing_parts[1].parse::<f32>(),
                    candidate.rows.parse::<f32>(),
                    candidate.columns.parse::<f32>(),
                ) {
                    format!("{:.1}x{:.1}", spacing_x * cols, spacing_y * rows)
                } else {
                    "N/A".to_string()
                }
            } else {
                "N/A".to_string()
            }
        } else {
            "N/A".to_string()
        };

        writeln!(
            csv_file,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            escape_csv_field(&candidate.series_description),
            escape_csv_field(&candidate.protocol_name),
            escape_csv_field(&candidate.series_number),
            escape_csv_field(&candidate.modality),
            escape_csv_field(&series_uid),
            escape_csv_field(&candidate.acquisition_type),
            escape_csv_field(&candidate.pixel_spacing),
            escape_csv_field(&candidate.slice_thickness),
            escape_csv_field(&candidate.spacing_between_slices),
            escape_csv_field(&fov),
            escape_csv_field(&candidate.repetition_time),
            escape_csv_field(&candidate.echo_time),
            escape_csv_field(&candidate.inversion_time),
            escape_csv_field(&candidate.flip_angle),
            escape_csv_field(&candidate.number_of_averages),
            escape_csv_field(&candidate.echo_train_length),
            escape_csv_field(&candidate.parallel_imaging_factor),
            escape_csv_field(&candidate.magnetic_field_strength),
            escape_csv_field(&candidate.image_type),
            escape_csv_field(&candidate.acquisition_time),
            escape_csv_field(&candidate.acquisition_duration),
            escape_csv_field(&candidate.derivation_description),
            escape_csv_field(&candidate.referenced_series_uid),
            file_count
        )?;
    }

    println!("\nMetadata exported to: {}", csv_path.display());
    Ok(())
}

/// Extract and organize DICOM files from ZIP archive by series
fn extract_and_organize_dicoms(
    zip_bytes: &[u8],
    deep_candidates: &[DeepDicomCandidate],
    output_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs;
    use std::io::Write;

    // Create output directory if it doesn't exist
    fs::create_dir_all(output_dir)?;

    // Create a map of file indices to their study/series info
    let mut file_organization = std::collections::HashMap::new();
    for candidate in deep_candidates {
        file_organization.insert(
            candidate.index,
            (
                candidate.study_instance_uid.clone(),
                candidate.series_instance_uid.clone(),
                candidate.name.clone(),
                candidate.study_description.clone(),
                candidate.series_description.clone(),
                candidate.series_number.clone(),
            ),
        );
    }

    // Open the ZIP archive
    let mut archive = ZipArchive::new(Cursor::new(zip_bytes))?;

    // Track created directories to avoid redundant filesystem calls
    let mut created_dirs = std::collections::HashSet::new();

    // Extract and organize files
    for (index, (study_uid, series_uid, original_name, study_desc, series_desc, series_num)) in
        file_organization
    {
        // Create directory structure with descriptive names
        // Format: Study_[Description]_[ShortUID] / Series_[Number]_[Description]_[ShortUID]
        let study_folder_name = if study_desc != "N/A" && !study_desc.is_empty() {
            format!(
                "Study_{}_{}",
                sanitize_filename(&study_desc),
                &study_uid.split('.').next_back().unwrap_or(&study_uid)
                    [..8.min(study_uid.split('.').next_back().unwrap_or(&study_uid).len())]
            )
        } else {
            format!(
                "Study_{}",
                &study_uid.split('.').next_back().unwrap_or(&study_uid)
                    [..16.min(study_uid.split('.').next_back().unwrap_or(&study_uid).len())]
            )
        };

        let series_folder_name = if series_desc != "N/A" && !series_desc.is_empty() {
            let series_prefix = if series_num != "N/A" && !series_num.is_empty() {
                format!(
                    "Series_{:04}_{}",
                    series_num.parse::<i32>().unwrap_or(0),
                    sanitize_filename(&series_desc)
                )
            } else {
                format!("Series_{}", sanitize_filename(&series_desc))
            };
            format!(
                "{}_{}",
                series_prefix,
                &series_uid.split('.').next_back().unwrap_or(&series_uid)[..8.min(
                    series_uid
                        .split('.')
                        .next_back()
                        .unwrap_or(&series_uid)
                        .len()
                )]
            )
        } else {
            format!(
                "Series_{}",
                &series_uid.split('.').next_back().unwrap_or(&series_uid)[..16.min(
                    series_uid
                        .split('.')
                        .next_back()
                        .unwrap_or(&series_uid)
                        .len()
                )]
            )
        };

        let study_dir = output_dir.join(study_folder_name);
        let series_dir = study_dir.join(series_folder_name);

        // Create directories if not already created
        if !created_dirs.contains(&series_dir) {
            fs::create_dir_all(&series_dir)?;
            created_dirs.insert(series_dir.clone());
        }

        // Extract file from ZIP
        let mut zip_file = archive.by_index(index)?;

        // Determine output filename (preserve original name if possible)
        let file_name = PathBuf::from(&original_name)
            .file_name()
            .and_then(|n| n.to_str())
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("dicom_{:04}.dcm", index));

        let output_path = series_dir.join(&file_name);

        // Read and write file
        let mut buffer = Vec::new();
        std::io::copy(&mut zip_file, &mut buffer)?;

        let mut output_file = fs::File::create(&output_path)?;
        output_file.write_all(&buffer)?;

        println!("Extracted: {} -> {}", original_name, output_path.display());
    }

    println!("\nExtraction complete!");
    println!("Files organized in: {}", output_dir.display());

    // Print summary of organization with descriptive names
    type SeriesMap = std::collections::HashMap<String, (String, String)>;
    type StudyMap = std::collections::HashMap<String, (String, SeriesMap)>;
    let mut study_info: StudyMap = StudyMap::new();

    for candidate in deep_candidates {
        let series_info = study_info
            .entry(candidate.study_instance_uid.clone())
            .or_insert_with(|| {
                (
                    candidate.study_description.clone(),
                    std::collections::HashMap::new(),
                )
            });

        series_info.1.insert(
            candidate.series_instance_uid.clone(),
            (
                candidate.series_description.clone(),
                candidate.series_number.clone(),
            ),
        );
    }

    println!("\nOrganization summary:");
    println!("  {} studies", study_info.len());
    for (study_uid, (study_desc, series_map)) in &study_info {
        let study_display = if study_desc != "N/A" && !study_desc.is_empty() {
            study_desc.clone()
        } else {
            format!(
                "UID: {}",
                &study_uid.split('.').next_back().unwrap_or(study_uid)
                    [..16.min(study_uid.split('.').next_back().unwrap_or(study_uid).len())]
            )
        };
        println!("    Study [{}]: {} series", study_display, series_map.len());
    }

    Ok(())
}

/// Sanitize a string to be used as a filename
fn sanitize_filename(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            c => c,
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let zip_path = args.file;

    // --derivations works with both ZIP files and directories
    if args.derivations {
        analyze_derivations(&zip_path)?;
        return Ok(());
    }

    let start = Instant::now();
    // Open the ZIP archive and load into memory
    let mut zip_data = Vec::new();
    std::fs::File::open(&zip_path)
        .expect("Failed to open ZIP file")
        .read_to_end(&mut zip_data)
        .expect("Failed to read ZIP file");
    let fr_duration = start.elapsed();

    let arc_data = Arc::new(zip_data);
    let candidates = scan_dicom_candidates_parallel(&arc_data)?;

    let duration = start.elapsed();

    let deep_candidates = deep_scan_dicom_candidates_parallel(&arc_data, args.mrn)?;
    let deep_duration = start.elapsed();

    // If --mrn flag is set, output only the MRN and exit
    if args.mrn {
        // Get unique patient IDs from all candidates
        let mut patient_ids = std::collections::HashSet::new();
        for cand in &deep_candidates {
            if cand.patient_id != "N/A" {
                patient_ids.insert(cand.patient_id.clone());
            }
        }

        // Print each unique patient ID
        for patient_id in patient_ids {
            println!("{}", patient_id);
        }

        return Ok(());
    }

    // Display results
    println!("Found {} DICOM files in archive:\n", candidates.len());
    for cand in &candidates {
        println!("{:<40}", cand.name);
    }
    println!("\nZIP file read took: {:?}", fr_duration);
    //println!("ZIP archive inspection took: {:?}", idx_duration);
    println!("DICOM detection took: {:?}", duration - fr_duration);
    println!("Total time: {:?}", duration);
    println!(
        "Effective scan rate of {:?} files/sec",
        (candidates.len() as f64 / duration.as_secs_f64()).round()
    );
    let total_compressed: u64 = candidates.iter().map(|c| c.compressed_size).sum();
    let total_uncompressed: u64 = candidates.iter().map(|c| c.uncompressed_size).sum();

    println!("Total compressed size:   {} bytes", total_compressed);
    println!("Total uncompressed size: {} bytes", total_uncompressed);

    if total_compressed > 0 {
        let ratio = total_uncompressed as f64 / total_compressed as f64;
        println!("Compression ratio: {:.3} %", ratio * 100.0);
    } else {
        println!("Compression ratio: N/A (no compressed data)");
    }

    println!("Using {} threads", rayon::current_num_threads());

    println!("Deep scan time = {:?}", deep_duration - duration);
    println!("Deep scan found {} DICOM files", deep_candidates.len());

    // For each unique study_instance_uid, find all the candidates that contain it
    let mut study_instance_uid_map = std::collections::HashMap::new();
    for cand in &deep_candidates {
        study_instance_uid_map
            .entry(cand.study_instance_uid.clone())
            .or_insert_with(std::collections::HashSet::new)
            .insert(cand.series_instance_uid.clone());
    }
    // Print the study_instance_uid and the number of candidates that contain it
    for (study_instance_uid, candidates) in study_instance_uid_map {
        println!(
            "Study Instance UID: {} has {} distinct series",
            study_instance_uid,
            candidates.len()
        );
    }

    for cand in &deep_candidates {
        println!(
            "{:<40} {} [{}] {}, {}",
            cand.name,
            cand.modality,
            cand.manufacturer,
            cand.study_instance_uid,
            cand.series_instance_uid
        );
    }

    // Extract XProtocol data if --xprot directory is specified
    if let Some(xprot_dir) = args.xprot {
        println!("\n--- Extracting XProtocol data from Siemens DICOM files ---");
        let xprot_count = extract_xprotocol_from_zip(&arc_data, &xprot_dir)?;
        println!(
            "Extracted {} XProtocol file(s) to {}",
            xprot_count,
            xprot_dir.display()
        );
    }

    // Extract and organize files if output directory is specified
    if let Some(output_dir) = args.output {
        println!("\n--- Extracting and organizing DICOM files ---");
        extract_and_organize_dicoms(&arc_data, &deep_candidates, &output_dir)?;

        // Export metadata to CSV
        export_series_metadata_csv(&deep_candidates, &output_dir)?;
    }

    Ok(())
}

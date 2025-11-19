use std::{
    io::{Cursor, Read},
    path::PathBuf,
    str::FromStr,
    sync::Arc,
    time::{Duration, Instant},
};

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

        // get the study instance UID
        let study_instance_uid = dcm_object
            .element(tags::STUDY_INSTANCE_UID)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        // get the series instance UID
        let series_instance_uid = dcm_object
            .element(tags::SERIES_INSTANCE_UID)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        // get the patient ID (MRN)
        let patient_id = dcm_object
            .element(PATIENT_ID)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        // get the sop instance UID
        let sop_instance_uid = dcm_object
            .element(tags::SOP_INSTANCE_UID)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        // get the Modality
        let modality = dcm_object
            .element(tags::MODALITY)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        let series_description = dcm_object
            .element(tags::SERIES_DESCRIPTION)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        let _series_date = dcm_object
            .element(tags::SERIES_DATE)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        let series_number = dcm_object
            .element(tags::SERIES_NUMBER)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        let _series_time = dcm_object
            .element(tags::SERIES_TIME)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        // get the Manufacturer
        let manufacturer = dcm_object
            .element(tags::MANUFACTURER)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        // get the SOP Class UID
        let sop_class_uid = dcm_object
            .element(tags::SOP_CLASS_UID)
            .map_or("N/A".to_string(), |e| {
                e.value().to_str().unwrap().to_string()
            });

        if !suppress_output {
            println!("sop_class_uid: {}", sop_class_uid);
        }

        if sop_class_uid == *"1.2.840.10008.5.1.4.1.1.4.1" && modality == *"MR" {
            if !suppress_output {
                println!("This is an enhanced MR image DICOM file");
            }

            let _acquisition_number = dcm_object
                .element(tags::ACQUISITION_NUMBER)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _acquisiton_date_time = dcm_object
                .element(tags::ACQUISITION_DATE_TIME)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _content_qualification = dcm_object
                .element(tags::CONTENT_QUALIFICATION)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _resonant_nucleus = dcm_object
                .element(tags::RESONANT_NUCLEUS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _kspace_filtering = dcm_object
                .element(tags::K_SPACE_FILTERING)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _magnetic_field_strength = dcm_object
                .element(tags::MAGNETIC_FIELD_STRENGTH)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _applicable_safety_standard_agency = dcm_object
                .element(tags::APPLICABLE_SAFETY_STANDARD_AGENCY)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _applicable_safety_standard_description = dcm_object
                .element(tags::APPLICABLE_SAFETY_STANDARD_DESCRIPTION)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _image_comments = dcm_object
                .element(tags::IMAGE_COMMENTS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _isocenter_position = dcm_object
                .element(tags::ISOCENTER_POSITION)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _b1rms = dcm_object
                .element(tags::B1RMS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _acquisition_contrast = dcm_object
                .element(tags::ACQUISITION_CONTRAST)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            // check the geometry stuff
            let mr_fov_geometry_sequence = dcm_object
                .element(tags::MRFOV_GEOMETRY_SEQUENCE)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _inplane_phase_encoding_direction = dcm_object
                .element(tags::IN_PLANE_PHASE_ENCODING_DIRECTION)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let mr_acquisition_frequency_encoding_steps = dcm_object
                .element(tags::MR_ACQUISITION_FREQUENCY_ENCODING_STEPS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let mr_acquisition_phase_encoding_steps_inplane = dcm_object
                .element(tags::MR_ACQUISITION_PHASE_ENCODING_STEPS_IN_PLANE)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let mr_acquisition_phase_encoding_steps_outofplane = dcm_object
                .element(tags::MR_ACQUISITION_PHASE_ENCODING_STEPS_OUT_OF_PLANE)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _percent_sampling = dcm_object
                .element(tags::PERCENT_SAMPLING)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _percent_phase_field_of_view = dcm_object
                .element(tags::PERCENT_PHASE_FIELD_OF_VIEW)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

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
            // get the TE (echo time)
            let te = dcm_object
                .element(tags::ECHO_TIME)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });
            // get the TR (repetition time)
            let tr = dcm_object
                .element(tags::REPETITION_TIME)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            // get the matrix size

            // get the slice thickness

            // get the pixel spacing

            // get the number of slices

            // get the pixel bit depth

            // get the SAR value
            let sar = dcm_object
                .element(tags::SAR)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            // get the dB/dt value
            let _db_dt = dcm_object
                .element(tags::D_BDT)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            // get the isocenter position
            let _isocenter_position = dcm_object
                .element(tags::ISOCENTER_POSITION)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            // get the receive coil name
            let receive_coil_name = dcm_object
                .element(tags::RECEIVE_COIL_NAME)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            // get the pixel bandwidth
            let pixel_bandwidth = dcm_object
                .element(tags::PIXEL_BANDWIDTH)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _number_pe = dcm_object
                .element(tags::NUMBER_OF_PHASE_ENCODING_STEPS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let acq_matrix = dcm_object
                .element(tags::ACQUISITION_MATRIX)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            // phase encoding direction, redundant with acq_matrix
            let phase_encoding_direction = dcm_object
                .element(tags::IN_PLANE_PHASE_ENCODING_DIRECTION)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let reconstruction_diameter = dcm_object
                .element(tags::RECONSTRUCTION_DIAMETER)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let pixel_spacing = dcm_object
                .element(tags::PIXEL_SPACING)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let rows = dcm_object
                .element(tags::ROWS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let columns = dcm_object
                .element(tags::COLUMNS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _b1_rms = dcm_object
                .element(tags::B1RMS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _bits_allocated = dcm_object
                .element(tags::BITS_ALLOCATED)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _bits_stored = dcm_object
                .element(tags::BITS_STORED)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _high_bit = dcm_object
                .element(tags::HIGH_BIT)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let scanning_sequence = dcm_object
                .element(tags::SCANNING_SEQUENCE)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let sequence_variant = dcm_object
                .element(tags::SEQUENCE_VARIANT)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let scan_options = dcm_object
                .element(tags::SCAN_OPTIONS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let mr_acquisition_type = dcm_object
                .element(tags::MR_ACQUISITION_TYPE)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _inversion_time = dcm_object
                .element(tags::INVERSION_TIME)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _sequence_name = dcm_object
                .element(tags::SEQUENCE_NAME)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let center_to_center_slice_gap = dcm_object
                .element(tags::SPACING_BETWEEN_SLICES)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let percent_sampling = dcm_object
                .element(tags::PERCENT_SAMPLING)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let percent_phase_fov = dcm_object
                .element(tags::PERCENT_PHASE_FIELD_OF_VIEW)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let flip_angle = dcm_object
                .element(tags::FLIP_ANGLE)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let _variable_flip_flag = dcm_object
                .element(tags::VARIABLE_FLIP_ANGLE_FLAG)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

            let slice_thickness = dcm_object
                .element(tags::SLICE_THICKNESS)
                .map_or("N/A".to_string(), |e| {
                    e.value().to_str().unwrap().to_string()
                });

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
                let internal_sequence_name = dcm_object
                    .element(Tag(0x0019, 0x109E))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                // this tag is "FL" as VR (single float)
                let acquisition_duration = dcm_object
                    .element(Tag(0x0019, 0x105A))
                    .map_or(f32::NAN, |e| e.value().to_float32().unwrap());

                let number_of_echoes = dcm_object
                    .element(Tag(0x0019, 0x107E))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });
                let _table_delta = dcm_object
                    .element(Tag(0x0019, 0x107F))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _gehc_private_creator_id = dcm_object
                    .element(Tag(0x0043, 0x0010))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _bitmap_of_prescan_options = dcm_object
                    .element(Tag(0x0043, 0x1001))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _gradient_offset_x = dcm_object
                    .element(Tag(0x0043, 0x1002))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _gradient_offset_y = dcm_object
                    .element(Tag(0x0043, 0x1003))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _gradient_offset_z = dcm_object
                    .element(Tag(0x0043, 0x1004))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _image_is_original = dcm_object
                    .element(Tag(0x0043, 0x1005))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _number_of_epi_shots = dcm_object
                    .element(Tag(0x0043, 0x1006))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _views_per_segment = dcm_object
                    .element(Tag(0x0043, 0x1007))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _respiratory_rate_bpm = dcm_object
                    .element(Tag(0x0043, 0x1008))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _respiratory_trigger_point = dcm_object
                    .element(Tag(0x0043, 0x1009))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _type_of_receiver_used = dcm_object
                    .element(Tag(0x0043, 0x100A))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _peak_dbdt = dcm_object
                    .element(Tag(0x0043, 0x100B))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _dbdt_limits_percent = dcm_object
                    .element(Tag(0x0043, 0x100C))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _psd_estimatated_limit = dcm_object
                    .element(Tag(0x0043, 0x100D))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _psd_estimated_limit_tps = dcm_object
                    .element(Tag(0x0043, 0x100E))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _sar_avg_head = dcm_object
                    .element(Tag(0x0043, 0x100F))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _application_name = dcm_object
                    .element(Tag(0x0043, 0x1077))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _application_version = dcm_object
                    .element(Tag(0x0043, 0x1078))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let _slices_per_volume = dcm_object
                    .element(Tag(0x0043, 0x1079))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                let asset_r_factors = dcm_object
                    .element(Tag(0x0043, 0x1083))
                    .map_or("N/A".to_string(), |e| {
                        e.value().to_str().unwrap().to_string()
                    });

                // note the acquisition duration is in micro seconds
                if !suppress_output {
                    println!(
                        "{} {:#?} {} {}",
                        internal_sequence_name,
                        Duration::from_micros(acquisition_duration as u64),
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let zip_path = args.file;

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
    Ok(())
}

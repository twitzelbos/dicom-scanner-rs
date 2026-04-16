use dicom::core::{Tag, VR};
use dicom::object::{open_file, InMemDicomObject};
use dicom::core::value::Value;
use std::env;
use std::process;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::convert::TryInto;

/// Parse DICOM sequence structure and extract nested content
fn parse_dicom_sequence(data: &[u8]) -> Option<String> {
    let mut offset = 0;

    // Check for sequence item delimiter (0xFFFE, 0xE000)
    if data.len() < 8 || data[0] != 0xFE || data[1] != 0xFF || data[2] != 0x00 || data[3] != 0xE0 {
        return None;
    }
    offset += 4;

    // Read item length (4 bytes, little-endian)
    let item_length = u32::from_le_bytes(data[offset..offset+4].try_into().ok()?) as usize;
    offset += 4;

    // Parse tags within the sequence item
    let item_end = if item_length == 0xFFFFFFFF {
        // Undefined length - search for item delimiter
        data.len()
    } else {
        offset + item_length
    };

    while offset + 8 <= item_end && offset + 8 <= data.len() {
        // Read tag (group, element)
        let group = u16::from_le_bytes(data[offset..offset+2].try_into().ok()?);
        let element = u16::from_le_bytes(data[offset+2..offset+4].try_into().ok()?);
        offset += 4;

        // Check for item/sequence delimiters
        if group == 0xFFFE {
            if element == 0xE00D || element == 0xE0DD {
                // Item or Sequence delimitation
                offset += 4; // Skip length field
                continue;
            }
        }

        // Read VR (2 bytes)
        let vr_bytes = if offset + 2 <= data.len() {
            [data[offset], data[offset + 1]]
        } else {
            break;
        };
        let vr = String::from_utf8_lossy(&vr_bytes).to_string();
        offset += 2;


        // Check if this looks like a valid VR (uppercase letters)
        let valid_vr = vr_bytes[0].is_ascii_uppercase() && vr_bytes[1].is_ascii_uppercase();

        let value_length = if !valid_vr {
            // Implicit VR - the "VR" bytes are actually the first 2 bytes of a 4-byte length
            offset -= 2; // Back up
            if offset + 4 > data.len() { break; }
            let len = u32::from_le_bytes(data[offset..offset+4].try_into().ok()?) as usize;
            offset += 4;
            len
        } else if vr == "OB" || vr == "OW" || vr == "SQ" || vr == "UN" || vr == "UT" {
            // These VRs have 2 reserved bytes then 4-byte length
            offset += 2; // Skip reserved bytes
            if offset + 4 > data.len() { break; }
            u32::from_le_bytes(data[offset..offset+4].try_into().ok()?) as usize
        } else {
            // Other VRs have 2-byte length
            if offset + 2 > data.len() { break; }
            u16::from_le_bytes(data[offset..offset+2].try_into().ok()?) as usize
        };

        if !valid_vr {
            // Already advanced offset above
        } else if vr == "OB" || vr == "OW" || vr == "SQ" || vr == "UN" || vr == "UT" {
            offset += 4;
        } else {
            offset += 2;
        }

        // Only print details for significant tags
        if group == 0x0021 && (element == 0x0010 || element == 0x1019) {
            println!("  Found tag ({:04x},{:04x}) Length={}", group, element, value_length);
        }

        // Check if this tag contains the XProtocol
        if offset + value_length <= data.len() {
            let tag_data = &data[offset..offset+value_length];

            // Look for XProtocol marker in this tag's data
            if let Some(xprot_pos) = tag_data.windows(11)
                .position(|window| window.starts_with(b"<XProtocol>")) {
                // Extract content and trim trailing null bytes
                let content = &tag_data[xprot_pos..];
                let trimmed_len = content.iter().rposition(|&b| b != 0).map(|i| i + 1).unwrap_or(0);
                return Some(String::from_utf8_lossy(&content[..trimmed_len]).to_string());
            }

            // For private creator tags, show the content
            if group == 0x0021 && element == 0x0010 {
                let creator = String::from_utf8_lossy(tag_data).trim_end_matches('\0').to_string();
                println!("  Private Creator: {}", creator);
            }
        }

        offset += value_length;
    }

    None
}

/// Check if data starts with CSA header signature "SV10"
fn is_csa_header(data: &[u8]) -> bool {
    data.len() >= 4 && &data[0..4] == b"SV10"
}

/// Parse CSA header and extract XProtocol content
fn parse_csa_header(data: &[u8]) -> Option<String> {
    if !is_csa_header(data) {
        return None;
    }

    // Skip signature (4 bytes) and unused bytes (4 bytes)
    let mut offset = 8;

    if data.len() < offset + 4 {
        return None;
    }

    // Read number of elements (little-endian u32)
    let num_elements = u32::from_le_bytes(data[offset..offset+4].try_into().ok()?) as usize;
    offset += 4;

    // Parse each CSA element
    for _ in 0..num_elements {
        if data.len() < offset + 64 + 4 + 4 + 4 + 4 {
            break;
        }

        // Read element name (64 bytes, null-terminated)
        let name_bytes = &data[offset..offset+64];
        let name = String::from_utf8_lossy(name_bytes)
            .trim_end_matches('\0')
            .to_string();
        offset += 64;

        // Skip VM (4 bytes)
        offset += 4;

        // Read VR (4 bytes) - usually like "OB\0\0" or "UN\0\0"
        let _vr_bytes = &data[offset..offset+4];
        offset += 4;

        // Skip sync pattern (4 bytes)
        offset += 4;

        // Read length (little-endian u32)
        let length = u32::from_le_bytes(data[offset..offset+4].try_into().ok()?) as usize;
        offset += 4;

        if data.len() < offset + length {
            break;
        }

        // Check if this element contains XProtocol data
        // Common names: "MrPhoenixProtocol", "MrProtocol", "MrEvaProtocol"
        if name.contains("Protocol") || name == "MrPhoenixProtocol" {
            let content_data = &data[offset..offset+length];

            // The content might still have some binary padding, so let's try to find
            // the start of the actual XProtocol text (usually starts with "### ASCCONV")
            if let Some(start_pos) = content_data.windows(10)
                .position(|window| window.starts_with(b"### ASCCONV") || window.starts_with(b"<XProtocol>")) {
                let content = &content_data[start_pos..];
                let trimmed_len = content.iter().rposition(|&b| b != 0).map(|i| i + 1).unwrap_or(0);
                return Some(String::from_utf8_lossy(&content[..trimmed_len]).to_string());
            } else {
                // If we can't find markers, try to extract as much readable text as possible
                let trimmed_len = content_data.iter().rposition(|&b| b != 0).map(|i| i + 1).unwrap_or(0);
                return Some(String::from_utf8_lossy(&content_data[..trimmed_len]).to_string());
            }
        }

        offset += length;
    }

    None
}

fn extract_tag(obj: &InMemDicomObject, tag_to_find: Tag, path: &str) -> Option<String> {
    // First try to find tag in the root
    if let Ok(tag_element) = obj.element(tag_to_find) {
        if let Ok(bytes) = tag_element.to_bytes() {
            // Special handling for tag (0021,10fe) which may contain CSA header or DICOM sequence
            if tag_to_find == Tag(0x0021, 0x10fe) {
                // Check for DICOM sequence format
                if bytes.len() >= 4 && bytes[0] == 0xFE && bytes[1] == 0xFF && bytes[2] == 0x00 && bytes[3] == 0xE0 {
                    if let Some(parsed) = parse_dicom_sequence(&bytes) {
                        return Some(parsed);
                    }
                } else if is_csa_header(&bytes) {
                    if let Some(parsed) = parse_csa_header(&bytes) {
                        return Some(parsed);
                    }
                }
            }
            // Default: return as string, trimming trailing null bytes
            let trimmed_len = bytes.iter().rposition(|&b| b != 0).map(|i| i + 1).unwrap_or(0);
            return Some(String::from_utf8_lossy(&bytes[..trimmed_len]).to_string());
        }
    }

    // Search through all elements for sequences
    for element in obj.iter() {
        if element.vr() == VR::SQ {
            let value = element.value();
            if let Value::Sequence(sequence) = value {
                for (i, item) in sequence.items().iter().enumerate() {
                    let item_path = format!("{}/{}[{}]", path, element.header().tag, i);
                    if let Some(tag_value) = extract_tag(item, tag_to_find, &item_path) {
                        return Some(tag_value);
                    }
                }
            }
        }
    }

    None
}

fn get_series_description(obj: &InMemDicomObject) -> String {
    // Try to get Series Description (0008,103e)
    if let Ok(tag) = obj.element(Tag(0x0008, 0x103e)) {
        if let Ok(bytes) = tag.to_bytes() {
            let desc = String::from_utf8_lossy(&bytes).trim().to_string();
            if !desc.is_empty() {
                return sanitize_filename(&desc);
            }
        }
    }

    // Fallback to Series Number (0020,0011)
    if let Ok(tag) = obj.element(Tag(0x0020, 0x0011)) {
        if let Ok(bytes) = tag.to_bytes() {
            let series_num = String::from_utf8_lossy(&bytes).trim().to_string();
            if !series_num.is_empty() {
                return format!("Series_{}", sanitize_filename(&series_num));
            }
        }
    }

    "Unknown_Series".to_string()
}

fn sanitize_filename(name: &str) -> String {
    name.chars()
        .map(|c| match c {
            '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            c if c.is_control() => '_',
            c => c,
        })
        .collect::<String>()
        .trim()
        .to_string()
}

fn process_file(file_path: &Path, output_dir: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let obj = open_file(file_path)?;

    let tag_1019_value = extract_tag(&obj, Tag(0x0021, 0x1019), "root");
    let tag_10fe_value = extract_tag(&obj, Tag(0x0021, 0x10fe), "root");

    if tag_1019_value.is_some() || tag_10fe_value.is_some() {
        let series_name = get_series_description(&obj);
        let output_filename = format!("{}.xprot", series_name);
        let output_path = output_dir.join(output_filename);

        let mut file = File::create(&output_path)?;

        // Write tag (0021,1019) if present
        if let Some(value) = &tag_1019_value {
            write!(file, "{}", value)?;
            println!("Extracted (0021,1019) from {}", file_path.display());
        }

        // Write tag (0021,10fe) if present
        if let Some(value) = &tag_10fe_value {
            if tag_1019_value.is_some() {
                writeln!(file)?; // Add empty line as separator
                writeln!(file)?; // Add second line for clarity when both tags present
            }
            write!(file, "{}", value)?;
            println!("Extracted (0021,10fe) from {}", file_path.display());
        }

        println!("Saved to {}", output_path.display());
    } else {
        println!("No (0021,1019) or (0021,10fe) tags found in {}", file_path.display());
    }

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.len() > 3 {
        eprintln!("Usage: {} <INPUT_PATH> [OUTPUT_DIR]", args[0]);
        eprintln!("  INPUT_PATH: DICOM file or directory containing DICOM files");
        eprintln!("  OUTPUT_DIR: Directory to save .xprot files (default: current directory)");
        process::exit(1);
    }

    let input_path = Path::new(&args[1]);
    let output_dir = if args.len() == 3 {
        Path::new(&args[2])
    } else {
        Path::new(".")
    };

    // Create output directory if it doesn't exist
    if let Err(e) = fs::create_dir_all(output_dir) {
        eprintln!("Error creating output directory '{}': {}", output_dir.display(), e);
        process::exit(1);
    }

    if input_path.is_file() {
        // Process single file
        if let Err(e) = process_file(input_path, output_dir) {
            eprintln!("Error processing file '{}': {}", input_path.display(), e);
            process::exit(1);
        }
    } else if input_path.is_dir() {
        // Process all files in directory
        let entries = match fs::read_dir(input_path) {
            Ok(entries) => entries,
            Err(e) => {
                eprintln!("Error reading directory '{}': {}", input_path.display(), e);
                process::exit(1);
            }
        };

        let mut processed_count = 0;
        let mut error_count = 0;

        for entry in entries {
            let entry = match entry {
                Ok(entry) => entry,
                Err(e) => {
                    eprintln!("Error reading directory entry: {}", e);
                    error_count += 1;
                    continue;
                }
            };

            let file_path = entry.path();
            if file_path.is_file() {
                match process_file(&file_path, output_dir) {
                    Ok(()) => processed_count += 1,
                    Err(e) => {
                        eprintln!("Error processing file '{}': {}", file_path.display(), e);
                        error_count += 1;
                    }
                }
            }
        }

        println!("Processing complete. Files processed: {}, Errors: {}", processed_count, error_count);
    } else {
        eprintln!("Input path '{}' is neither a file nor a directory", input_path.display());
        process::exit(1);
    }
}

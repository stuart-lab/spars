use crate::io_utils;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::fs;

pub fn subset_matrix(
    input: Option<String>,
    output: &str,
    rows_file: Option<String>,
    cols_file: Option<String>,
    no_reindex: bool,
) -> Result<(), Box<dyn Error>> {

    // Check if both rows_file and cols_file are None
    if rows_file.is_none() && cols_file.is_none() {
        return Err("At least one of --rows or --cols must be specified.".into());
    }

    // check if input is a file or directory
    let is_directory = input.as_ref().map_or(false, |path| Path::new(path).is_dir());

    // Read the indices to retain
    let rows_to_retain = if let Some(ref file_path) = rows_file {
        read_indices(file_path)?
    } else {
        Vec::new() // Empty vector indicates all rows are retained
    };

    let cols_to_retain = if let Some(ref file_path) = cols_file {
        read_indices(file_path)?
    } else {
        Vec::new() // Empty vector indicates all columns are retained
    };

    let row_mapping = if !rows_to_retain.is_empty() && !no_reindex {
        Some(create_index_mapping(&rows_to_retain))
    } else {
        None
    };

    let col_mapping = if !cols_to_retain.is_empty() && !no_reindex {
        Some(create_index_mapping(&cols_to_retain))
    } else {
        None
    };

    let rows_set: HashSet<usize> = if !rows_to_retain.is_empty() {
        rows_to_retain.clone().into_iter().collect()
    } else {
        HashSet::new() // Empty set indicates all rows are retained
    };

    let cols_set: HashSet<usize> = if !cols_to_retain.is_empty() {
        cols_to_retain.clone().into_iter().collect()
    } else {
        HashSet::new() // Empty set indicates all columns are retained
    };

    if is_directory {
        fs::create_dir_all(output)?;
        let dir_path = Path::new(input.as_ref().unwrap());
        let output_path = Path::new(output);
        let matrix_file = dir_path.join("matrix.mtx.gz").to_string_lossy().into_owned();
        let output_matrix = output_path.join("matrix.mtx.gz").to_string_lossy().into_owned();
        let features_file = dir_path.join("features.tsv.gz").to_string_lossy().into_owned();
        let barcodes_file = if dir_path.join("barcodes.tsv.gz").exists() {
            dir_path.join("barcodes.tsv.gz").to_string_lossy().into_owned()
        } else {
            dir_path.join("barcodes.tsv").to_string_lossy().into_owned()
        };

        let output_features = output_path.join("features.tsv.gz").to_string_lossy().into_owned();
        let output_barcodes = output_path.join("barcodes.tsv.gz").to_string_lossy().into_owned();

        subset_mtx_file(&matrix_file, &output_matrix, &row_mapping, &col_mapping, &rows_set, &cols_set)?;
        subset_dimnames_file(&features_file, &output_features, &rows_to_retain, &row_mapping)?;
        subset_dimnames_file(&barcodes_file, &output_barcodes, &cols_to_retain, &col_mapping)?;

    } else {
        let matrix_file = input.as_ref().unwrap();
        subset_mtx_file(&matrix_file, output, &row_mapping, &col_mapping, &rows_set, &cols_set)?;
    }

    Ok(())
}

fn read_indices(file_path: &str) -> Result<Vec<usize>, Box<dyn Error>> {
    let reader = io_utils::get_reader(file_path)?;
    let mut indices = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed_line = line.trim();

        if trimmed_line.is_empty() {
            continue;
        }

        match trimmed_line.parse::<usize>() {
            Ok(idx) => indices.push(idx),
            Err(_) => {
                eprintln!(
                    "Warning: Failed to parse index on line {} in '{}'",
                    line_number + 1,
                    file_path
                );
            }
        }
    }

    Ok(indices)
}

fn subset_mtx_file(
    input_file: &str,
    output_file: &str,
    row_mapping: &Option<HashMap<usize, usize>>,
    col_mapping: &Option<HashMap<usize, usize>>,
    rows_set: &HashSet<usize>,
    cols_set: &HashSet<usize>,
) -> Result<(), Box<dyn Error>> {
    let temp_data_file = "temp_data.mtx";

    // Open the input and temporary data files
    let mut reader = io_utils::get_reader(input_file)?;
    let mut temp_writer = io_utils::get_writer(&temp_data_file)?;

    let mut header_lines = Vec::new();
    let mut n_nonzeros = 0;
    let mut data_started = false;

    let mut buffer = String::new();

    let mut nz_elements = 0;

    while reader.read_line(&mut buffer)? > 0 {
        let trimmed_line = buffer.trim();

        if trimmed_line.starts_with('%') || trimmed_line.is_empty() {
            header_lines.push(buffer.clone());
            buffer.clear();
            continue;
        }

        if !data_started {
            // This is the line with dimensions and n_nonzeros, skip for now
            data_started = true;
            // Store the original dimensions for later use
            header_lines.push(buffer.clone());
            buffer.clear();
            continue;
        }

        nz_elements += 1;
        if nz_elements % 1_000_000 == 0 {
            print!("\rProcessed {} M elements", nz_elements / 1_000_000);
            std::io::stdout().flush().expect("Can't flush output");
        }

        // Parse data lines
        let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
        if parts.len() >= 2 {
            let row_idx: usize = parts[0].parse().unwrap();
            let col_idx: usize = parts[1].parse().unwrap();
            let value = if parts.len() >= 3 {
                parts[2]
            } else {
                "1"
            }; // Default value for pattern matrices

            let row_included = if rows_set.is_empty() {
                true // All rows are retained
            } else {
                rows_set.contains(&row_idx)
            };

            let col_included = if cols_set.is_empty() {
                true // All columns are retained
            } else {
                cols_set.contains(&col_idx)
            };

            if row_included && col_included {
                n_nonzeros += 1;

                let new_row = if let Some(ref mapping) = row_mapping {
                    mapping.get(&row_idx).unwrap_or(&row_idx)
                } else {
                    &row_idx
                };

                let new_col = if let Some(ref mapping) = col_mapping {
                    mapping.get(&col_idx).unwrap_or(&col_idx)
                } else {
                    &col_idx
                };

                writeln!(temp_writer, "{} {} {}", new_row, new_col, value)?;
            }
        }

        buffer.clear();
    }

    temp_writer.flush()?;

    // Clear the progress line after processing is complete
    println!("");

    // Now write the header to the output file
    let mut output_writer = io_utils::get_writer(output_file)?;

    // Write the header lines excluding the last one
    for (index, line) in header_lines.iter().enumerate() {
        if index < header_lines.len() - 1 {
            writeln!(output_writer, "{}", line.trim_end())?;
        }
    }

    // Calculate new dimensions
    let (orig_n_rows, orig_n_cols, _orig_n_nonzeros) = parse_header_line(&header_lines)?;

    let n_rows = if let Some(ref mapping) = row_mapping {
        mapping.len()
    } else if !rows_set.is_empty() {
        rows_set.len()
    } else {
        orig_n_rows // All rows retained
    };

    let n_cols = if let Some(ref mapping) = col_mapping {
        mapping.len()
    } else if !cols_set.is_empty() {
        cols_set.len()
    } else {
        orig_n_cols // All columns retained
    };

    // record program information in the header
    writeln!(output_writer, "%metadata json: {{\"software_version\": \"spars-1.0.0\", \"command\": \"spars subset\"}}")?;
    writeln!(output_writer, "{} {} {}", n_rows, n_cols, n_nonzeros)?;

    output_writer.flush()?;

    // Concatenate the temp data file to the output file
    let mut temp_reader = io_utils::get_reader(&temp_data_file)?;
    io::copy(&mut temp_reader, &mut output_writer)?;

    output_writer.flush()?;

    // Clean up the temporary data file
    std::fs::remove_file(temp_data_file)?;

    Ok(())
}


fn subset_dimnames_file(
    input_file: &str,
    output_file: &str,
    indices_to_retain: &[usize],
    mapping: &Option<HashMap<usize, usize>>,
) -> Result<(), Box<dyn Error>> {
    let reader = io_utils::get_reader(input_file)?;
    let mut writer = io_utils::get_writer(output_file)?;

    let indices_set: HashSet<usize> = indices_to_retain.iter().cloned().collect();

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        if indices_set.is_empty() || indices_set.contains(&(index + 1)) {
            let _ = if let Some(ref mapping) = mapping {
                *mapping.get(&(index + 1)).unwrap_or(&(index + 1))
            } else {
                index + 1
            };
            writeln!(writer, "{}", line)?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn create_index_mapping(indices: &[usize]) -> HashMap<usize, usize> {
    indices
        .iter()
        .enumerate()
        .map(|(new_idx, &old_idx)| (old_idx, new_idx + 1))
        .collect()
}

fn parse_header_line(header_lines: &Vec<String>) -> Result<(usize, usize, usize), Box<dyn Error>> {
    // Find the last non-comment line in header_lines
    for line in header_lines.iter().rev() {
        let trimmed_line = line.trim();
        if trimmed_line.starts_with('%') || trimmed_line.is_empty() {
            continue;
        } else {
            let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
            if parts.len() >= 3 {
                let n_rows = parts[0].parse()?;
                let n_cols = parts[1].parse()?;
                let n_nonzeros = parts[2].parse()?;
                return Ok((n_rows, n_cols, n_nonzeros));
            } else {
                break;
            }
        }
    }
    Err("Error parsing header line for matrix dimensions.".into())
}
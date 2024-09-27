use crate::io_utils;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::io::{self, BufRead, Write};

pub fn subset_matrix(
    input_file: &str,
    output_file: &str,
    rows_file: Option<String>,
    cols_file: Option<String>,
    no_reindex: bool,
) -> Result<(), Box<dyn Error>> {
    // Read the indices to retain
    let rows_to_retain = if let Some(ref file_path) = rows_file {
        read_indices(file_path)?
    } else {
        Vec::new()
    };

    let cols_to_retain = if let Some(ref file_path) = cols_file {
        read_indices(file_path)?
    } else {
        Vec::new()
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
        rows_to_retain.into_iter().collect()
    } else {
        HashSet::new()
    };

    let cols_set: HashSet<usize> = if !cols_to_retain.is_empty() {
        cols_to_retain.into_iter().collect()
    } else {
        HashSet::new()
    };

    let temp_data_file = "temp_data.mtx";

    // Open the input and temporary data files
    let mut reader = io_utils::get_reader(input_file)?;
    let mut temp_writer = io_utils::get_writer(&temp_data_file)?;

    let mut header_lines = Vec::new();
    let mut n_nonzeros = 0;
    let mut data_started = false;

    // Read and process the input file line by line
    let mut buffer = String::new();
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
            buffer.clear();
            continue;
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

            let row_included = rows_set.is_empty() || rows_set.contains(&row_idx);
            let col_included = cols_set.is_empty() || cols_set.contains(&col_idx);

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

    // Now write the header to the output file
    let mut output_writer = io_utils::get_writer(output_file)?;

    // Write the header lines
    for line in &header_lines {
        writeln!(output_writer, "{}", line.trim_end())?;
    }

    // Write the new dimensions and n_nonzeros
    let n_rows = if let Some(ref mapping) = row_mapping {
        mapping.len()
    } else {
        // If no reindexing, use original dimensions
        header_lines
            .last()
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .parse()
            .unwrap()
    };

    let n_cols = if let Some(ref mapping) = col_mapping {
        mapping.len()
    } else {
        header_lines
            .last()
            .unwrap()
            .split_whitespace()
            .nth(1)
            .unwrap()
            .parse()
            .unwrap()
    };

    writeln!(output_writer, "{} {} {}", n_rows, n_cols, n_nonzeros)?;

    output_writer.flush()?;

    // Concatenate the temp data file to the output file
    let mut temp_reader = io_utils::get_reader(&temp_data_file)?;
    io::copy(&mut temp_reader, &mut output_writer)?;

    output_writer.flush()?;

    // Clean up the temporary data file
    std::fs::remove_file(temp_data_file)?;

    println!("Subsetted matrix has been written to '{}'", output_file);

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

fn create_index_mapping(indices: &[usize]) -> HashMap<usize, usize> {
    indices
        .iter()
        .enumerate()
        .map(|(new_idx, &old_idx)| (old_idx, new_idx + 1))
        .collect()
}
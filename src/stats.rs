use crate::io_utils;
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, Write};

pub fn compute_stats(input_file: &str, output_prefix: &str) -> Result<(), Box<dyn Error>> {
    let row_output_file = format!("{}_row.tsv", output_prefix);
    let col_output_file = format!("{}_col.tsv", output_prefix);

    let mut reader = io_utils::get_reader(input_file)?;

    let mut n_rows = 0;
    let mut n_cols = 0;

    let mut line = String::new();
    let mut header_lines = Vec::new();

    // Read header
    loop {
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        let trimmed_line = line.trim();

        if trimmed_line.starts_with('%') || trimmed_line.is_empty() {
            header_lines.push(line.clone());
            line.clear();
            continue; // Skip comments and empty lines
        } else {
            // Parse dimensions
            let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
            if parts.len() >= 3 {
                n_rows = parts[0].parse()?;
                n_cols = parts[1].parse()?;
                header_lines.push(line.clone());
                line.clear();
                break;
            } else {
                eprintln!("Error: Invalid header line: '{}'", line);
                return Ok(());
            }
        }
    }

    let mut row_stats: HashMap<usize, Stats> = HashMap::new();
    let mut col_stats: HashMap<usize, Stats> = HashMap::new();
    let mut row_nonzero_counts: HashMap<usize, usize> = HashMap::new();
    let mut col_nonzero_counts: HashMap<usize, usize> = HashMap::new();

    // Process data entries
    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed_line = line.trim();

        if trimmed_line.starts_with('%') || trimmed_line.is_empty() {
            continue; // Skip comments and empty lines
        }

        if line_number % 1_000_000 == 0 {
            print!("\rProcessed {} M elements", line_number / 1_000_000);
            std::io::stdout().flush().expect("Can't flush output");
        }

        let parts: Vec<&str> = trimmed_line.split_whitespace().collect();

        if parts.len() < 3 {
            eprintln!(
                "Warning: Line {} is malformed (expected at least 3 columns): '{}'",
                line_number + 1,
                line
            );
            continue;
        }

        let row_idx = match parts[0].parse::<usize>() {
            Ok(idx) => idx,
            Err(_) => {
                eprintln!(
                    "Warning: Failed to parse row index on line {}: '{}'",
                    line_number + 1,
                    line
                );
                continue;
            }
        };

        let col_idx = match parts[1].parse::<usize>() {
            Ok(idx) => idx,
            Err(_) => {
                eprintln!(
                    "Warning: Failed to parse column index on line {}: '{}'",
                    line_number + 1,
                    line
                );
                continue;
            }
        };

        let value = match parts[2].parse::<f64>() {
            Ok(val) => val,
            Err(_) => {
                eprintln!(
                    "Warning: Failed to parse value on line {}: '{}'",
                    line_number + 1,
                    line
                );
                continue;
            }
        };

        // Update row statistics
        row_stats
            .entry(row_idx)
            .and_modify(|stats| stats.update(value))
            .or_insert_with(|| Stats::new(value));

        *row_nonzero_counts.entry(row_idx).or_insert(0) += 1;

        // Update column statistics
        col_stats
            .entry(col_idx)
            .and_modify(|stats| stats.update(value))
            .or_insert_with(|| Stats::new(value));

        *col_nonzero_counts.entry(col_idx).or_insert(0) += 1;
    }

    // Finalize stats for all rows
    for row_idx in 1..=n_rows {
        let nonzero_count = row_nonzero_counts.get(&row_idx).cloned().unwrap_or(0);
        row_stats
            .entry(row_idx)
            .and_modify(|stats| stats.finalize(n_cols, nonzero_count))
            .or_insert_with(|| {
                let mut stats = Stats::new(0.0);
                stats.count = n_cols;
                stats.nonzero_count = 0;
                stats
            });
    }

    // Finalize stats for all columns
    for col_idx in 1..=n_cols {
        let nonzero_count = col_nonzero_counts.get(&col_idx).cloned().unwrap_or(0);
        col_stats
            .entry(col_idx)
            .and_modify(|stats| stats.finalize(n_rows, nonzero_count))
            .or_insert_with(|| {
                let mut stats = Stats::new(0.0);
                stats.count = n_rows;
                stats.nonzero_count = 0;
                stats
            });
    }

    // Write row statistics to TSV file
    let mut row_file = io_utils::get_writer(&row_output_file)?;
    writeln!(row_file, "Row\tNonZeroCount\tSum\tMean\tVariance\tStdDev\tMin\tMax")?;
    let mut sorted_rows: Vec<_> = row_stats.iter().collect();
    sorted_rows.sort_by_key(|&(idx, _)| *idx);
    for (&idx, stats) in sorted_rows {
        writeln!(
            row_file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            idx,
            stats.nonzero_count,
            stats.sum,
            format!("{:.4}", stats.mean()),
            format!("{:.4}", stats.variance()),
            format!("{:.4}", stats.std_dev()),
            stats.min,
            stats.max
        )?;
    }

    // Write column statistics to TSV file
    let mut col_file = io_utils::get_writer(&col_output_file)?;
    writeln!(
        col_file,
        "Column\tNonZeroCount\tSum\tMean\tVariance\tStdDev\tMin\tMax"
    )?;
    let mut sorted_cols: Vec<_> = col_stats.iter().collect();
    sorted_cols.sort_by_key(|&(idx, _)| *idx);
    for (&idx, stats) in sorted_cols {
        writeln!(
            col_file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            idx,
            stats.nonzero_count,
            stats.sum,
            format!("{:.4}", stats.mean()),
            format!("{:.4}", stats.variance()),
            format!("{:.4}", stats.std_dev()),
            stats.min,
            stats.max
        )?;
    }

    Ok(())
}

struct Stats {
    count: usize,
    nonzero_count: usize,
    sum: f64,
    sum_of_squares: f64,
    min: f64,
    max: f64,
}

impl Stats {
    fn new(value: f64) -> Self {
        Stats {
            count: 0, // Will be set during finalization
            nonzero_count: 1,
            sum: value,
            sum_of_squares: value * value,
            min: value,
            max: value,
        }
    }

    fn update(&mut self, value: f64) {
        self.nonzero_count += 1;
        self.sum += value;
        self.sum_of_squares += value * value;

        // Update min and max
        if value < self.min {
            self.min = value;
        }
        if value > self.max {
            self.max = value;
        }
    }

    fn finalize(&mut self, total_count: usize, nonzero_count: usize) {
        // Adjust the count to include zeros
        self.count = total_count;
        self.nonzero_count = nonzero_count;

        // Zeros contribute zero to sum and sum_of_squares
        // No adjustment needed
    }

    fn mean(&self) -> f64 {
        self.sum / self.count as f64
    }

    fn variance(&self) -> f64 {
        // Variance = [sum_of_squares - n * mean^2] / n
        let mean = self.mean();
        (self.sum_of_squares - self.count as f64 * mean * mean) / self.count as f64
    }

    fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }
}
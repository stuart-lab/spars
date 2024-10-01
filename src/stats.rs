use crate::io_utils;
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, Write};
use std::cmp::Ordering;
use std::str::FromStr;

#[derive(Clone, Copy)]
enum SortColumn {
    Index,
    NonZeroCount,
    Sum,
    Mean,
    Variance,
    StdDev,
    Min,
    Max,
}

impl FromStr for SortColumn {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "index" => Ok(SortColumn::Index),
            "nonzerocount" => Ok(SortColumn::NonZeroCount),
            "sum" => Ok(SortColumn::Sum),
            "mean" => Ok(SortColumn::Mean),
            "variance" => Ok(SortColumn::Variance),
            "stddev" => Ok(SortColumn::StdDev),
            "min" => Ok(SortColumn::Min),
            "max" => Ok(SortColumn::Max),
            _ => Err(format!("Invalid sort column: {}", s)),
        }
    }
}

pub fn compute_stats(input_file: &str, output_prefix: &str, sort_by: Option<String>) -> Result<(), Box<dyn Error>> {

    let valid_sort_columns = vec!["NonZeroCount", "Sum", "Mean", "Variance", "StdDev", "Min", "Max"];
    if let Some(ref sort_by) = sort_by {
        if !valid_sort_columns.contains(&sort_by.as_str()) {
            return Err(format!("Invalid sort column: {:?}", sort_by).into());
        }
    }

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

    // Clear the progress line after processing is complete
    println!("");

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

    // Sort and write row statistics
    let row_output_file = format!("{}_row.tsv", output_prefix);
    write_stats(&row_output_file, &row_stats, sort_by.as_deref(), "Row")?;

    // Sort and write column statistics
    let col_output_file = format!("{}_col.tsv", output_prefix);
    write_stats(&col_output_file, &col_stats, sort_by.as_deref(), "Column")?;

    Ok(())
}

fn write_stats(
    output_file: &str,
    stats: &HashMap<usize, Stats>,
    sort_by: Option<&str>,
    index_name: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = io_utils::get_writer(output_file)?;
    writeln!(file, "{}\tNonZeroCount\tSum\tMean\tVariance\tStdDev\tMin\tMax", index_name)?;

    let mut sorted_stats: Vec<(usize, &Stats)> = stats.iter().map(|(&k, v)| (k, v)).collect();
    
    if let Some(sort_column) = sort_by.and_then(|s| s.parse::<SortColumn>().ok()) {
        sort_stats(&mut sorted_stats, sort_column);
    } else {
        sorted_stats.sort_by_key(|&(idx, _)| idx);
    }

    for (idx, stats) in &sorted_stats {
        writeln!(
            file,
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

fn sort_stats(stats: &mut [(usize, &Stats)], sort_column: SortColumn) {
    stats.sort_by(|&(idx_a, a), &(idx_b, b)| {
        let cmp = match sort_column {
            SortColumn::Index => idx_a.cmp(&idx_b),
            SortColumn::NonZeroCount => a.nonzero_count.cmp(&b.nonzero_count),
            SortColumn::Sum => a.sum.partial_cmp(&b.sum).unwrap_or(Ordering::Equal),
            SortColumn::Mean => a.mean().partial_cmp(&b.mean()).unwrap_or(Ordering::Equal),
            SortColumn::Variance => a.variance().partial_cmp(&b.variance()).unwrap_or(Ordering::Equal),
            SortColumn::StdDev => a.std_dev().partial_cmp(&b.std_dev()).unwrap_or(Ordering::Equal),
            SortColumn::Min => a.min.partial_cmp(&b.min).unwrap_or(Ordering::Equal),
            SortColumn::Max => a.max.partial_cmp(&b.max).unwrap_or(Ordering::Equal),
        };
        cmp.reverse()
    });
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
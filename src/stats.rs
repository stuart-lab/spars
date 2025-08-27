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
    PearsonResidualVar,
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
            "pearsonresidualvar" => Ok(SortColumn::PearsonResidualVar),
            _ => Err(format!("Invalid sort column: {}", s)),
        }
    }
}

pub fn compute_stats(input_file: &str, output_prefix: &str, sort_by: Option<String>, theta: Option<f64>) -> Result<(), Box<dyn Error>> {
    let default_theta = 100.0;
    let theta = theta.unwrap_or(default_theta);

    let valid_sort_columns = vec!["NonZeroCount", "Sum", "Mean", "Variance", "StdDev", "Min", "Max", "PearsonResidualVar"];
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

    // Second pass to compute residuals
    let mut reader = io_utils::get_reader(input_file)?;
    
    // Skip header
    let mut line = String::new();
    loop {
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        let trimmed_line = line.trim();

        if !trimmed_line.starts_with('%') && !trimmed_line.is_empty() {
            let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
            if parts.len() >= 3 {
                line.clear();
                break;
            }
        }
        line.clear();
    }

    let mut row_pearson_residual_sums: HashMap<usize, f64> = HashMap::new();
    let mut row_pearson_residual_squares: HashMap<usize, f64> = HashMap::new();
    
    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed_line = line.trim();

        if trimmed_line.starts_with('%') || trimmed_line.is_empty() {
            continue; // Skip comments and empty lines
        }

        if line_number % 1_000_000 == 0 {
            print!("\rComputing residuals: {} M elements", line_number / 1_000_000);
            std::io::stdout().flush().expect("Can't flush output");
        }

        let parts: Vec<&str> = trimmed_line.split_whitespace().collect();
        if parts.len() < 3 {
            continue;
        }

        let row_idx = match parts[0].parse::<usize>() {
            Ok(idx) => idx,
            Err(_) => continue,
        };

        let value = match parts[2].parse::<f64>() {
            Ok(val) => val,
            Err(_) => continue,
        };

        // Get the row mean from previously computed stats
        if let Some(row_stat) = row_stats.get(&row_idx) {
            let mean = row_stat.mean();
            
            // Calculate Pearson residual: (X - μ) / sqrt(μ + μ²/θ)
            let denominator = (mean + mean * mean / theta).sqrt();
            let pearson_residual = if denominator != 0.0 { (value - mean) / denominator } else { 0.0 };
            
            // Clip Pearson residuals by sqrt(n_cols)
            let clip_threshold = (n_cols as f64).sqrt();
            let clipped_pearson_residual = pearson_residual.max(-clip_threshold).min(clip_threshold);
            
            // Update Pearson residual sums for variance calculation
            *row_pearson_residual_sums.entry(row_idx).or_insert(0.0) += clipped_pearson_residual;
            *row_pearson_residual_squares.entry(row_idx).or_insert(0.0) += clipped_pearson_residual * clipped_pearson_residual;
        }
    }
    
    println!("");  // Clear the progress line

    // Add contributions from implicit zeros to Pearson residual variance
    for (row_idx, stats) in &row_stats {
        let num_zeros = stats.count - stats.nonzero_count;
        if num_zeros > 0 {
            let mean = stats.mean();
            if mean > 0.0 {  // Only process if mean > 0 (avoid division by zero)
                let denominator = (mean + mean * mean / theta).sqrt();
                let zero_residual = -mean / denominator;
                
                // Clip by sqrt(n_cols) as done for non-zeros
                let clip_threshold = (n_cols as f64).sqrt();
                let clipped_residual = zero_residual.max(-clip_threshold).min(clip_threshold);
                
                // Add contribution of all zeros to the sums
                let zero_contribution = num_zeros as f64 * clipped_residual;
                let zero_squares_contribution = num_zeros as f64 * clipped_residual * clipped_residual;
                
                *row_pearson_residual_sums.entry(*row_idx).or_insert(0.0) += zero_contribution;
                *row_pearson_residual_squares.entry(*row_idx).or_insert(0.0) += zero_squares_contribution;
            }
        }
    }

    // Calculate residual variances for each row
    for (row_idx, stats) in &mut row_stats {
        // Calculate Pearson residual variance
        if let (Some(sum), Some(sum_squares)) = (row_pearson_residual_sums.get(row_idx), row_pearson_residual_squares.get(row_idx)) {
            if stats.count > 1 {
                let mean = sum / stats.count as f64;
                stats.pearson_residual_var = (sum_squares - stats.count as f64 * mean * mean) / (stats.count as f64 - 1.0);
            } else {
                stats.pearson_residual_var = 0.0;
            }
        } else {
            stats.pearson_residual_var = 0.0;
        }
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
    
    // Different headers for rows and columns
    if index_name == "Row" {
        writeln!(file, "{}\tNonZeroCount\tSum\tMean\tVariance\tStdDev\tMin\tMax\tPearsonResidualVar", index_name)?;
    } else {
        writeln!(file, "{}\tNonZeroCount\tSum\tMean\tVariance\tStdDev\tMin\tMax", index_name)?;
    }

    let mut sorted_stats: Vec<(usize, &Stats)> = stats.iter().map(|(&k, v)| (k, v)).collect();
    
    if let Some(sort_column) = sort_by.and_then(|s| s.parse::<SortColumn>().ok()) {
        sort_stats(&mut sorted_stats, sort_column);
    } else {
        sorted_stats.sort_by_key(|&(idx, _)| idx);
    }

    for (idx, stats) in &sorted_stats {
        if index_name == "Row" {
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                idx,
                stats.nonzero_count,
                stats.sum,
                format!("{:.4}", stats.mean()),
                format!("{:.4}", stats.variance()),
                format!("{:.4}", stats.std_dev()),
                stats.min,
                stats.max,
                format!("{:.4}", stats.pearson_residual_var),
            )?;
        } else {
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
            SortColumn::PearsonResidualVar => a.pearson_residual_var.partial_cmp(&b.pearson_residual_var).unwrap_or(Ordering::Equal),
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
    pearson_residual_var: f64,
}

impl Stats {
    fn new(value: f64) -> Self {
        Stats {
            count: 1,
            nonzero_count: if value != 0.0 { 1 } else { 0 },
            sum: value,
            sum_of_squares: value * value,
            min: value,
            max: value,
            pearson_residual_var: 0.0,
        }
    }

    fn update(&mut self, value: f64) {
        self.count += 1;
        if value != 0.0 {
            self.nonzero_count += 1;
        }
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

    // fn finalize(&mut self, total_count: usize, nonzero_count: usize) {
    //     // Adjust the count and nonzero_count
    //     self.count = total_count;
    //     self.nonzero_count = nonzero_count;

    //     // No need to adjust sum or sum_of_squares as zeros don't contribute

    //     // min is nonzero minumum
    // }

    fn finalize(&mut self, total_count: usize, nonzero_count: usize) {
    // Adjust the count and nonzero_count
    self.count = total_count;
    self.nonzero_count = nonzero_count;

    // If there are any zeros (implicit), min should be 0
    if nonzero_count < total_count {
        self.min = 0.0;
    }
    // Note: max stays as is - the maximum non-zero value is correct
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_stats_handles_implicit_zeros() {
        // Create a Stats struct with explicit values only
        let mut stats = Stats::new(5.0);
        stats.update(10.0);
        
        // Before finalize: only tracked 2 values
        assert_eq!(stats.count, 2);
        assert_eq!(stats.min, 5.0);  // Wrong! Missing zeros
        
        // After finalize with 10 total elements (8 implicit zeros)
        stats.finalize(10, 2);
        
        assert_eq!(stats.count, 10);
        assert_eq!(stats.nonzero_count, 2);
        assert_eq!(stats.min, 0.0);  // Should now be 0
        assert_eq!(stats.max, 10.0);  // Max unchanged
        
        // Mean should be (5+10)/10 = 1.5
        assert!((stats.mean() - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_compute_stats_small_matrix() {
        // Create a small test matrix
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, 
            "%%MatrixMarket matrix coordinate real general\n\
             3 4 3\n\
             1 1 6.0\n\
             1 2 6.0\n\
             2 3 12.0").unwrap();
        
        let temp_path = temp_file.path().to_str().unwrap();
        
        // Run stats computation
        compute_stats(temp_path, "test_output", Some("PearsonResidualVar".to_string()), Some(100.0))
            .expect("Stats computation failed");
        
        // Check the output files exist
        assert!(std::path::Path::new("test_output_row.tsv").exists());
        assert!(std::path::Path::new("test_output_col.tsv").exists());
        
        // Clean up
        std::fs::remove_file("test_output_row.tsv").ok();
        std::fs::remove_file("test_output_col.tsv").ok();
    }

    #[test]
    fn test_pearson_residual_calculation() {
        // Test the math for Pearson residuals with zeros
        let mean: f64 = 1.2;
        let theta: f64 = 100.0;
        let n_cols: f64 = 10.0;  // Make this f64 too
        
        let denominator = (mean + mean * mean / theta).sqrt();
        let zero_residual = -mean / denominator;
        let clip_threshold = n_cols.sqrt();  // Now works since n_cols is f64
        
        // Zero residual should be negative
        assert!(zero_residual < 0.0);
        
        // Should be within clipping threshold
        assert!(zero_residual.abs() <= clip_threshold);
    }
}
use crate::io_utils;
use crate::loess;
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, Write};
use rand::thread_rng;
use rand::seq::SliceRandom;


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

    // Collect means and variances for all rows
    let mut row_indices = Vec::new();
    let mut row_means = Vec::new();
    let mut row_variances = Vec::new();
    let mut row_means_log = Vec::new();
    let mut row_variances_log = Vec::new();

    for row_idx in 1..=n_rows {
        if let Some(stats) = row_stats.get(&row_idx) {
            row_indices.push(row_idx);
            let mean = stats.mean();
            let variance = stats.variance();
            
            // Log-transform mean and variance, handling zero values
            let log_mean = if mean > 0.0 { mean.ln() } else { f64::NEG_INFINITY };
            let log_variance = if variance > 0.0 { variance.ln() } else { f64::NEG_INFINITY };
            
            row_means.push(mean);
            row_variances.push(variance);
            row_means_log.push(log_mean);
            row_variances_log.push(log_variance);
        }
    }

    // Compute residual variances with LOESS
    println!("");
    println!("Fitting Mean-Variance model");
    let span = 0.3;
    let residual_variances = compute_residual_variance(&row_means_log, &row_variances_log, span)?;

    // Update the Stats structs with residual variances
    for (&row_idx, &res_var) in row_indices.iter().zip(residual_variances.iter()) {
        if let Some(stats) = row_stats.get_mut(&row_idx) {
            stats.residual_variance = Some(res_var);
        }
    }

    // Write row statistics to TSV file
    let mut row_file = io_utils::get_writer(&row_output_file)?;
    writeln!(
        row_file,
        "Row\tNonZeroCount\tSum\tMean\tVariance\tStdDev\tMin\tMax\tResidualVariance"
    )?;
    let mut sorted_rows: Vec<_> = row_stats.iter().collect();
    sorted_rows.sort_by_key(|&(idx, _)| *idx);
    for (&idx, stats) in sorted_rows {
        writeln!(
            row_file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            idx,
            stats.nonzero_count,
            stats.sum,
            format!("{:.4}", stats.mean()),
            format!("{:.4}", stats.variance()),
            format!("{:.4}", stats.std_dev()),
            stats.min,
            stats.max,
            match stats.residual_variance {
                Some(rv) => format!("{:.4}", rv),
                None => "NA".to_string(),
            }
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
    residual_variance: Option<f64>,
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
            residual_variance: None,
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

fn compute_residual_variance(
    means: &Vec<f64>,
    variances: &Vec<f64>,
    span: f64,
) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    if means.len() != variances.len() {
        return Err("Means and variances vectors must have the same length.".into());
    }

    // Filter out -inf values
    let filtered: Vec<(f64, f64, usize)> = means.iter()
        .zip(variances.iter())
        .enumerate()
        .filter(|(_, (&m, &v))| m.is_finite() && v.is_finite())
        .map(|(i, (&m, &v))| (m, v, i))
        .collect();

    let mut filtered_means = Vec::new();
    let mut filtered_variances = Vec::new();
    let mut indices = Vec::new();
    for (mean, variance, index) in filtered {
        filtered_means.push(mean);
        filtered_variances.push(variance);
        indices.push(index);
    }

    // Check if we should use binning or fit the model to all values
    let use_binning = false;

    if !use_binning {
        // If not using binning, directly apply LOESS to all filtered data
        let x: Vec<f64> = filtered_means.clone();
        let y: Vec<f64> = filtered_variances.clone();
        let fitted = loess::loess(&x, &y, span, 2, &x);

        // Create the result vector with the same length as the original input
        let mut result = vec![f64::NAN; means.len()];
        for (&idx, &fitted_value) in indices.iter().zip(fitted.iter()) {
            result[idx] = fitted_value;
        }

        return Ok(result);
    }

    // Compute min and max of means
    let min_mean = filtered_means.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_mean = filtered_means.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Define the number of bins
    let num_bins = 200;
    let bin_width = (max_mean - min_mean) / num_bins as f64;

    println!("{}", bin_width);

    // Create bins
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); num_bins];

    for (idx, &mean_value) in filtered_means.iter().enumerate() {
        let bin_idx = ((mean_value - min_mean) / bin_width).floor() as usize;
        let bin_idx = bin_idx.min(num_bins - 1); // Ensure the index is within bounds
        bins[bin_idx].push(idx);
    }

    // Sample up to 500 instances from each bin
    let mut sampled_indices = Vec::new();
    let mut rng = thread_rng();

    for bin in bins {
        let n_in_bin = bin.len();
        if n_in_bin == 0 {
            continue; // Skip empty bins
        }
        if n_in_bin <= 500 {
            sampled_indices.extend(bin);
        } else {
            let sampled = bin.choose_multiple(&mut rng, 500).cloned().collect::<Vec<usize>>();
            sampled_indices.extend(sampled);
        }
    }

    // Create sampled mean and variance vectors
    let sampled_means = sampled_indices.iter().map(|&i| filtered_means[i]).collect::<Vec<f64>>();
    let sampled_variances = sampled_indices.iter().map(|&i| filtered_variances[i]).collect::<Vec<f64>>();

    let degree = 2; // quadratic fitting
    let var_fitted = loess::loess(&sampled_means, &sampled_variances, span, degree, &filtered_means);

    // Compute residuals for filtered values
    let filtered_residuals: Vec<f64> = filtered_variances
        .iter()
        .zip(var_fitted.iter())
        .map(|(&yi, &y_fit)| yi - y_fit)
        .collect();

    // Create full residual variance vector, setting 0 for -inf values
    let mut residual_variance = vec![0.0; means.len()];
    for ((&r, &idx), &var) in filtered_residuals.iter().zip(indices.iter()).zip(filtered_variances.iter()) {
        residual_variance[idx] = if var > 0.0 { r * r / var } else { 0.0 };
    }

    Ok(residual_variance)
}
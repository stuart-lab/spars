use nalgebra::{DMatrix, DVector};


fn weighted_linear_regression(
    x: &[f64],
    y: &[f64],
    weights: &[f64],
) -> (f64, f64) {
    let n = x.len();

    // Construct weighted design matrix and response vector
    let mut w_sqrt = Vec::with_capacity(n);
    for &w in weights {
        w_sqrt.push(w.sqrt());
    }

    // Design matrix with ones (intercept) and x
    let mut x_matrix = DMatrix::from_element(n, 2, 1.0);
    for i in 0..n {
        x_matrix[(i, 1)] = x[i];
    }

    // Apply weights
    for i in 0..n {
        let w = w_sqrt[i];
        x_matrix.row_mut(i).scale_mut(w);
    }

    // Response vector
    let mut y_vector = DVector::from_column_slice(y);
    for i in 0..n {
        y_vector[i] *= w_sqrt[i];
    }

    // Solve normal equations using QR decomposition
    let qr = x_matrix.qr();
    let result = qr.solve(&y_vector);

    if let Some(coeffs) = result {
        let intercept = coeffs[0];
        let slope = coeffs[1];
        (intercept, slope)
    } else {
        // Handle singularity
        (0.0, 0.0)
    }
}

fn weighted_quadratic_regression(
    x: &[f64],
    y: &[f64],
    weights: &[f64],
) -> (f64, f64, f64) {
    let n = x.len();

    // Construct weighted design matrix and response vector
    let mut w_sqrt = Vec::with_capacity(n);
    for &w in weights {
        w_sqrt.push(w.sqrt());
    }

    // Design matrix with ones (intercept), x, and x^2
    let mut x_matrix = DMatrix::from_element(n, 3, 1.0);
    for i in 0..n {
        x_matrix[(i, 1)] = x[i];
        x_matrix[(i, 2)] = x[i] * x[i];
    }

    // Apply weights
    for i in 0..n {
        let w = w_sqrt[i];
        x_matrix.row_mut(i).scale_mut(w);
    }

    // Response vector
    let mut y_vector = DVector::from_column_slice(y);
    for i in 0..n {
        y_vector[i] *= w_sqrt[i];
    }

    // Solve least squares problem using SVD
    let svd = x_matrix.svd(true, true);
    let result = svd.solve(&y_vector, 1e-6);

    if let Ok(coeffs) = result {
        let intercept = coeffs[0];
        let linear = coeffs[1];
        let quadratic = coeffs[2];
        (intercept, linear, quadratic)
    } else {
        // Handle singularity
        (0.0, 0.0, 0.0)
    }
}

pub fn loess(
    x: &[f64],
    y: &[f64],
    span: f64,
    degree: usize, // 1 for linear, 2 for quadratic
    x_predict: &[f64], // Points at which to predict
) -> Vec<f64> {
    let n = x.len();
    let mut y_predict = Vec::with_capacity(x_predict.len());

    let k = ((span * n as f64).ceil() as usize).max(degree + 1); // Ensure enough points

    for &x0 in x_predict {
        // Find distances from x0 to all x
        let mut distances: Vec<(usize, f64)> = x
            .iter()
            .enumerate()
            .map(|(i, &xi)| (i, (xi - x0).abs()))
            .collect();

        // Sort by distance
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // Select k nearest neighbors
        let neighbors: Vec<usize> = distances.iter().take(k).map(|&(i, _)| i).collect();

        // Maximum distance in neighborhood
        let max_distance = distances[k - 1].1.max(1e-12); // Avoid division by zero

        // Compute weights
        let weights: Vec<f64> = neighbors
            .iter()
            .map(|&i| tricube_weight((x[i] - x0).abs(), max_distance))
            .collect();

        // Get local x and y
        let x_local: Vec<f64> = neighbors.iter().map(|&i| x[i]).collect();
        let y_local: Vec<f64> = neighbors.iter().map(|&i| y[i]).collect();

        // Fit local model and predict at x0
        let y0 = if degree == 1 {
            let (intercept, slope) = weighted_linear_regression(&x_local, &y_local, &weights);
            intercept + slope * x0
        } else {
            let (a0, a1, a2) = weighted_quadratic_regression(&x_local, &y_local, &weights);
            a0 + a1 * x0 + a2 * x0 * x0
        };

        y_predict.push(y0);
    }

    y_predict
}

fn tricube_weight(distance: f64, max_distance: f64) -> f64 {
    let u = distance / max_distance;
    if u >= 1.0 {
        0.0
    } else {
        let tmp = 1.0 - u * u * u;
        tmp * tmp * tmp
    }
}
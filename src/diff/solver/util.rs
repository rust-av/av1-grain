use v_frame::plane::Plane;

use crate::diff::BLOCK_SIZE;

/// Solves Ax = b, where x and b are column vectors of size nx1 and A is nxn
#[allow(clippy::many_single_char_names)]
pub(super) fn linsolve(
    n: usize,
    a: &mut [f64],
    stride: usize,
    b: &mut [f64],
    x: &mut [f64],
) -> bool {
    // Forward elimination
    for k in 0..(n - 1) {
        // Bring the largest magnitude to the diagonal position
        ((k + 1)..n).rev().for_each(|i| {
            if a[(i - 1) * stride + k].abs() < a[i * stride + k].abs() {
                (0..n).for_each(|j| {
                    a.swap(i * stride + j, (i - 1) * stride + j);
                });
                b.swap(i, i - 1);
            }
        });

        for i in k..(n - 1) {
            if a[k * stride + k].abs() < f64::EPSILON {
                return false;
            }
            let c = a[(i + 1) * stride + k] / a[k * stride + k];
            (0..n).for_each(|j| {
                a[(i + 1) * stride + j] -= c * a[k * stride + j];
            });
            b[i + 1] -= c * b[k];
        }
    }

    // Backward substitution
    for i in (0..n).rev() {
        if a[i * stride + i].abs() < f64::EPSILON {
            return false;
        }
        let mut c = 0.0f64;
        for j in (i + 1)..n {
            c += a[i * stride + j] * x[j];
        }
        x[i] = (b[i] - c) / a[i * stride + i];
    }

    true
}

pub(super) fn multiply_mat(
    m1: &[f64],
    m2: &[f64],
    res: &mut [f64],
    m1_rows: usize,
    inner_dim: usize,
    m2_cols: usize,
) {
    assert!(res.len() >= m1_rows * m2_cols);
    assert!(m1.len() >= m1_rows * inner_dim);
    assert!(m2.len() >= m2_cols * inner_dim);
    let mut idx = 0;
    for row in 0..m1_rows {
        for col in 0..m2_cols {
            let mut sum = 0f64;
            for inner in 0..inner_dim {
                // SAFETY: We do the bounds checks once at the top to improve performance.
                unsafe {
                    sum += m1.get_unchecked(row * inner_dim + inner)
                        * m2.get_unchecked(inner * m2_cols + col);
                }
            }
            // SAFETY: We do the bounds checks once at the top to improve performance.
            unsafe {
                *res.get_unchecked_mut(idx) = sum;
            }
            idx += 1;
        }
    }
}

#[must_use]
pub(super) fn normalized_cross_correlation(a: &[f64], b: &[f64], n: usize) -> f64 {
    let mut c = 0f64;
    let mut a_len = 0f64;
    let mut b_len = 0f64;
    for i in 0..n {
        a_len += a[i].powi(2);
        b_len += b[i].powi(2);
        c += a[i] * b[i];
    }
    c / (a_len.sqrt() * b_len.sqrt())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn extract_ar_row(
    coords: &[[isize; 2]],
    num_coords: usize,
    source: &Plane<u8>,
    denoised: &Plane<u8>,
    alt_source: Option<&Plane<u8>>,
    alt_denoised: Option<&Plane<u8>>,
    x: usize,
    y: usize,
    buffer: &mut [f64],
) -> f64 {
    let source_data = source.data_origin();
    let denoised_data = denoised.data_origin();
    let stride = source.cfg.stride;
    for i in 0..num_coords {
        let x_i = x as isize + coords[i][0];
        let y_i = y as isize + coords[i][1];
        debug_assert!(x_i >= 0);
        debug_assert!(y_i >= 0);
        buffer[i] = f64::from(source_data[y_i as usize * stride + x_i as usize])
            - f64::from(denoised_data[y_i as usize * stride + x_i as usize]);
    }
    let val = f64::from(source.data_origin()[y * stride + x])
        - f64::from(denoised.data_origin()[y * stride + x]);

    if let Some(alt_source) = alt_source {
        if let Some(alt_denoised) = alt_denoised {
            let alt_source_data = alt_source.data_origin();
            let alt_denoised_data = alt_denoised.data_origin();

            let alt_stride = alt_source.cfg.stride;
            let mut source_sum = 0u64;
            let mut denoised_sum = 0u64;
            let mut num_samples = 0usize;
            for dy_i in 0..(1 << source.cfg.ydec) {
                let y_up = (y << source.cfg.ydec) + dy_i;
                for dx_i in 0..(1 << source.cfg.xdec) {
                    let x_up = (x << source.cfg.xdec) + dx_i;
                    source_sum += u64::from(alt_source_data[y_up * alt_stride + x_up]);
                    denoised_sum += u64::from(alt_denoised_data[y_up * alt_stride + x_up]);
                    num_samples += 1;
                }
            }
            buffer[num_coords] = (source_sum as f64 - denoised_sum as f64) / num_samples as f64;
        }
    }

    val
}

#[must_use]
pub(super) fn get_block_mean(
    source: &Plane<u8>,
    frame_dims: (usize, usize),
    x_o: usize,
    y_o: usize,
) -> f64 {
    let max_h = (frame_dims.1 - y_o).min(BLOCK_SIZE);
    let max_w = (frame_dims.0 - x_o).min(BLOCK_SIZE);
    let data = source.data_origin();
    assert!(data.len() >= (y_o + max_h - 1) * source.cfg.stride + x_o + max_w - 1);

    let mut block_sum = 0u64;
    for y in 0..max_h {
        for x in 0..max_w {
            let index = (y_o + y) * source.cfg.stride + x_o + x;
            // SAFETY: We bounds check once before entering the loop
            // to improve vectorization
            block_sum += unsafe { u64::from(*data.get_unchecked(index)) };
        }
    }

    block_sum as f64 / (max_w * max_h) as f64
}

#[must_use]
pub(super) fn get_noise_var(
    source: &Plane<u8>,
    denoised: &Plane<u8>,
    frame_dims: (usize, usize),
    x_o: usize,
    y_o: usize,
    block_w: usize,
    block_h: usize,
) -> f64 {
    let max_h = (frame_dims.1 - y_o).min(block_h);
    let max_w = (frame_dims.0 - x_o).min(block_w);

    let source_data = source.data_origin();
    let denoised_data = denoised.data_origin();

    let mut noise_var_sum = 0u64;
    let mut noise_sum = 0i64;
    for y in 0..max_h {
        for x in 0..max_w {
            let index = (y_o + y) * source.cfg.stride + x_o + x;
            let noise = i64::from(source_data[index]) - i64::from(denoised_data[index]);
            noise_sum += noise;
            noise_var_sum += noise.pow(2) as u64;
        }
    }

    let noise_mean = noise_sum as f64 / (max_w * max_h) as f64;
    noise_var_sum as f64 / (max_w * max_h) as f64 - noise_mean.powi(2)
}

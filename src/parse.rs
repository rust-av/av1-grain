// Copyright (c) 2022-2022, The rav1e contributors. All rights reserved
//
// This source code is subject to the terms of the BSD 2 Clause License and
// the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
// was not distributed with this source code in the LICENSE file, you can
// obtain it at www.aomedia.org/license/software. If the Alliance for Open
// Media Patent License 1.0 was not distributed with this source code in the
// PATENTS file, you can obtain it at www.aomedia.org/license/patent.

use arrayvec::ArrayVec;

use crate::{GrainTableSegment, NUM_UV_COEFFS, NUM_UV_POINTS, NUM_Y_COEFFS, NUM_Y_POINTS};

/// This file has the implementation details of the grain table.
///
/// The file format is an ascii representation for readability and
/// editability. Array parameters are separated from the non-array
/// parameters and prefixed with a few characters to make for easy
/// localization with a parameter set. Each entry is prefixed with "E"
/// and the other parameters are only specified if "apply-grain" is
/// non-zero.
///
/// ```text
/// filmgrn1
/// E <start-time> <end-time> <apply-grain> <random-seed> <dynamic-grain>
///  p <ar_coeff_lag> <ar_coeff_shift> <grain_scale_shift> ...
///  sY <num_y_points> <point_0_x> <point_0_y> ...
///  sCb <num_cb_points> <point_0_x> <point_0_y> ...
///  sCr <num_cr_points> <point_0_x> <point_0_y> ...
///  cY <ar_coeff_y_0> ....
///  cCb <ar_coeff_cb_0> ....
///  cCr <ar_coeff_cr_0> ....
/// E <start-time> ...
/// ```
///
/// # Errors
///
/// - If the file cannot be opened
/// - If the file does not contain a properly formatted film grain table
#[inline]
pub fn parse_grain_table(input: &str) -> anyhow::Result<Vec<GrainTableSegment>> {
    let lines = non_empty_lines(input);
    let Some((header, body)) = lines.split_first() else {
        anyhow::bail!("Expected filmgrn1 header");
    };

    let (header_tag, header_values) = header.split()?;
    anyhow::ensure!(
        header_tag == "filmgrn1" && header_values.is_empty(),
        "line {}: expected filmgrn1 header",
        header.line_no
    );

    parse_segments(body)
}

#[derive(Debug)]
struct LineFields<'a> {
    line_no: usize,
    fields: Vec<&'a str>,
}

impl<'a> LineFields<'a> {
    fn split(&self) -> anyhow::Result<(&'a str, &[&'a str])> {
        self.fields
            .split_first()
            .map(|(tag, values)| (*tag, values))
            .ok_or_else(|| anyhow::anyhow!("line {}: expected a non-empty line", self.line_no))
    }
}

fn non_empty_lines(input: &str) -> Vec<LineFields<'_>> {
    input
        .lines()
        .enumerate()
        .filter_map(|(index, line)| {
            let fields: Vec<_> = line.split_whitespace().collect();
            (!fields.is_empty()).then_some(LineFields {
                line_no: index + 1,
                fields,
            })
        })
        .collect()
}

fn expect_tag<'a>(line: &'a LineFields<'_>, tag: &str) -> anyhow::Result<&'a [&'a str]> {
    let (actual, values) = line.split()?;
    anyhow::ensure!(
        actual == tag,
        "line {}: expected {tag} line, got {actual}",
        line.line_no
    );
    Ok(values)
}

fn parse_values<T>(line: &LineFields<'_>, values: &[&str], label: &str) -> anyhow::Result<Vec<T>>
where
    T: std::str::FromStr,
    T::Err: std::fmt::Display,
{
    values
        .iter()
        .map(|value| {
            value.parse::<T>().map_err(|error| {
                anyhow::anyhow!(
                    "line {}: failed to parse {label} value `{value}`: {error}",
                    line.line_no
                )
            })
        })
        .collect()
}

fn parse_segments(lines: &[LineFields<'_>]) -> anyhow::Result<Vec<GrainTableSegment>> {
    let mut index = 0;
    let mut saw_segment = false;
    let mut segments = Vec::new();

    while let Some(line) = lines.get(index) {
        let (tag, _) = line.split()?;
        if tag != "E" {
            if saw_segment {
                break;
            }
            anyhow::bail!("line {}: expected E line, got {tag}", line.line_no);
        }

        saw_segment = true;
        let e_params = e_params(line)?;
        index += 1;

        if !e_params.apply {
            continue;
        }

        let p_params = p_params(required_line(lines, index, "p")?)?;
        index += 1;
        let s_y_params =
            s_params::<NUM_Y_POINTS>(required_line(lines, index, "sY")?, "sY", "Y-plane")?;
        index += 1;
        let s_cb_params =
            s_params::<NUM_UV_POINTS>(required_line(lines, index, "sCb")?, "sCb", "Cb-plane")?;
        index += 1;
        let s_cr_params =
            s_params::<NUM_UV_POINTS>(required_line(lines, index, "sCr")?, "sCr", "Cr-plane")?;
        index += 1;

        let coeff_count = (2 * p_params.ar_coeff_lag * (p_params.ar_coeff_lag + 1)) as usize;
        let c_y_params = c_params::<NUM_Y_COEFFS>(
            required_line(lines, index, "cY")?,
            "cY",
            "Y-plane",
            coeff_count,
        )?;
        index += 1;
        let c_cb_params = c_params::<NUM_UV_COEFFS>(
            required_line(lines, index, "cCb")?,
            "cCb",
            "Cb-plane",
            coeff_count + 1,
        )?;
        index += 1;
        let c_cr_params = c_params::<NUM_UV_COEFFS>(
            required_line(lines, index, "cCr")?,
            "cCr",
            "Cr-plane",
            coeff_count + 1,
        )?;
        index += 1;

        segments.push(GrainTableSegment {
            start_time: e_params.start,
            end_time: e_params.end,
            scaling_points_y: s_y_params,
            scaling_points_cb: s_cb_params,
            scaling_points_cr: s_cr_params,
            scaling_shift: p_params.scaling_shift,
            ar_coeff_lag: p_params.ar_coeff_lag,
            ar_coeffs_y: c_y_params,
            ar_coeffs_cb: c_cb_params,
            ar_coeffs_cr: c_cr_params,
            ar_coeff_shift: p_params.ar_coeff_shift,
            cb_mult: p_params.cb_mult,
            cb_luma_mult: p_params.cb_luma_mult,
            cb_offset: p_params.cb_offset,
            cr_mult: p_params.cr_mult,
            cr_luma_mult: p_params.cr_luma_mult,
            cr_offset: p_params.cr_offset,
            overlap_flag: p_params.overlap_flag,
            chroma_scaling_from_luma: p_params.chroma_scaling_from_luma,
            grain_scale_shift: p_params.grain_scale_shift,
            random_seed: e_params.seed,
        });
    }

    anyhow::ensure!(saw_segment, "Expected at least one grain table segment");
    Ok(segments)
}

fn required_line<'a, 'input>(
    lines: &'a [LineFields<'input>],
    index: usize,
    expected_tag: &str,
) -> anyhow::Result<&'a LineFields<'input>> {
    let Some(line) = lines.get(index) else {
        anyhow::bail!("Expected {expected_tag} line, got end of input");
    };
    let (actual_tag, _) = line.split()?;
    anyhow::ensure!(
        actual_tag == expected_tag,
        "line {}: expected {expected_tag} line, got {actual_tag}",
        line.line_no
    );
    Ok(line)
}

#[derive(Debug, Clone, Copy)]
struct EParams {
    pub start: u64,
    pub end: u64,
    pub apply: bool,
    pub seed: u16,
}

fn e_params(line: &LineFields<'_>) -> anyhow::Result<EParams> {
    let values = expect_tag(line, "E")?;
    let [start, end, apply, seed, _dynamic_grain] = values else {
        anyhow::bail!(
            "line {}: expected 5 values on E line, got {}",
            line.line_no,
            values.len()
        );
    };

    let start = start.parse::<u64>().map_err(|error| {
        anyhow::anyhow!("line {}: failed to parse start_time: {error}", line.line_no)
    })?;
    let end = end.parse::<u64>().map_err(|error| {
        anyhow::anyhow!("line {}: failed to parse end_time: {error}", line.line_no)
    })?;
    let apply = apply.parse::<u8>().map_err(|error| {
        anyhow::anyhow!(
            "line {}: failed to parse apply_grain: {error}",
            line.line_no
        )
    })? > 0;
    let seed = seed.parse::<u16>().map_err(|error| {
        anyhow::anyhow!(
            "line {}: failed to parse random_seed: {error}",
            line.line_no
        )
    })?;

    anyhow::ensure!(
        end >= start,
        "line {}: start time must be before end time",
        line.line_no
    );

    Ok(EParams {
        start,
        end,
        apply,
        seed,
    })
}

#[derive(Debug, Clone, Copy)]
struct PParams {
    ar_coeff_lag: u8,
    ar_coeff_shift: u8,
    grain_scale_shift: u8,
    scaling_shift: u8,
    chroma_scaling_from_luma: bool,
    overlap_flag: bool,
    cb_mult: u8,
    cb_luma_mult: u8,
    cb_offset: u16,
    cr_mult: u8,
    cr_luma_mult: u8,
    cr_offset: u16,
}

#[allow(clippy::too_many_lines)]
fn p_params(line: &LineFields<'_>) -> anyhow::Result<PParams> {
    let values = expect_tag(line, "p")?;
    let [
        ar_coeff_lag,
        ar_coeff_shift,
        grain_scale_shift,
        scaling_shift,
        chroma_scaling_from_luma,
        overlap_flag,
        cb_mult,
        cb_luma_mult,
        cb_offset,
        cr_mult,
        cr_luma_mult,
        cr_offset,
    ] = values
    else {
        anyhow::bail!(
            "line {}: expected 12 values on p line, got {}",
            line.line_no,
            values.len()
        );
    };

    let params = PParams {
        ar_coeff_lag: ar_coeff_lag.parse().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse ar_coeff_lag: {error}",
                line.line_no
            )
        })?,
        ar_coeff_shift: ar_coeff_shift.parse().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse ar_coeff_shift: {error}",
                line.line_no
            )
        })?,
        grain_scale_shift: grain_scale_shift.parse().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse grain_scale_shift: {error}",
                line.line_no
            )
        })?,
        scaling_shift: scaling_shift.parse().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse scaling_shift: {error}",
                line.line_no
            )
        })?,
        chroma_scaling_from_luma: chroma_scaling_from_luma.parse::<u8>().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse chroma_scaling_from_luma: {error}",
                line.line_no
            )
        })? > 0,
        overlap_flag: overlap_flag.parse::<u8>().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse overlap_flag: {error}",
                line.line_no
            )
        })? > 0,
        cb_mult: cb_mult.parse().map_err(|error| {
            anyhow::anyhow!("line {}: failed to parse cb_mult: {error}", line.line_no)
        })?,
        cb_luma_mult: cb_luma_mult.parse().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse cb_luma_mult: {error}",
                line.line_no
            )
        })?,
        cb_offset: cb_offset.parse().map_err(|error| {
            anyhow::anyhow!("line {}: failed to parse cb_offset: {error}", line.line_no)
        })?,
        cr_mult: cr_mult.parse().map_err(|error| {
            anyhow::anyhow!("line {}: failed to parse cr_mult: {error}", line.line_no)
        })?,
        cr_luma_mult: cr_luma_mult.parse().map_err(|error| {
            anyhow::anyhow!(
                "line {}: failed to parse cr_luma_mult: {error}",
                line.line_no
            )
        })?,
        cr_offset: cr_offset.parse().map_err(|error| {
            anyhow::anyhow!("line {}: failed to parse cr_offset: {error}", line.line_no)
        })?,
    };

    anyhow::ensure!(
        (8..=11).contains(&params.scaling_shift),
        "line {}: scaling_shift must be between 8 and 11",
        line.line_no
    );
    anyhow::ensure!(
        params.ar_coeff_lag <= 3,
        "line {}: ar_coeff_lag must be between 0 and 3",
        line.line_no
    );
    anyhow::ensure!(
        (6..=9).contains(&params.ar_coeff_shift),
        "line {}: ar_coeff_shift must be between 6 and 9",
        line.line_no
    );

    Ok(params)
}

fn s_params<const N: usize>(
    line: &LineFields<'_>,
    tag: &str,
    label: &str,
) -> anyhow::Result<ArrayVec<[u8; 2], N>> {
    let values = expect_tag(line, tag)?;
    let parsed = parse_values::<u8>(line, values, label)?;
    let Some((count, points)) = parsed.split_first() else {
        anyhow::bail!("line {}: expected {label} point count", line.line_no);
    };

    let count = *count as usize;
    let expected_values = count * 2;
    anyhow::ensure!(
        points.len() == expected_values,
        "line {}: expected {expected_values} {label} point values, got {}",
        line.line_no,
        points.len()
    );
    anyhow::ensure!(
        count <= N,
        "line {}: expected at most {N} {label} points, got {count}",
        line.line_no
    );

    let mut output = ArrayVec::new();
    for point in points.chunks_exact(2) {
        let [x, y] = point else {
            unreachable!("chunks_exact(2) only yields two-value chunks");
        };
        output.push([*x, *y]);
    }
    Ok(output)
}

fn c_params<const N: usize>(
    line: &LineFields<'_>,
    tag: &str,
    label: &str,
    expected_count: usize,
) -> anyhow::Result<ArrayVec<i8, N>> {
    let values = expect_tag(line, tag)?;
    let parsed = parse_values::<i8>(line, values, label)?;
    anyhow::ensure!(
        parsed.len() == expected_count,
        "line {}: expected {expected_count} {label} coeffs, got {}",
        line.line_no,
        parsed.len()
    );
    anyhow::ensure!(
        parsed.len() <= N,
        "line {}: expected at most {N} {label} coeffs, got {}",
        line.line_no,
        parsed.len()
    );

    let mut output = ArrayVec::new();
    for value in parsed {
        output.push(value);
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_luma_only_table() {
        // This is the luma-only table format generated by
        // both aomenc's photon noise utility and by av1an.
        let input = r#"filmgrn1
E 0 9223372036854775807 1 7391 1
  p 0 6 0 8 0 1 0 0 0 0 0 0
  sY 14  0 20 20 5 39 4 59 3 78 3 98 3 118 3 137 3 157 3 177 3 196 3 216 4 235 4 255 4
  sCb 0
  sCr 0
  cY
  cCb 0
  cCr 0
"#;
        let expected = GrainTableSegment {
            start_time: 0,
            end_time: 9_223_372_036_854_775_807,
            scaling_points_y: ArrayVec::from([
                [0, 20],
                [20, 5],
                [39, 4],
                [59, 3],
                [78, 3],
                [98, 3],
                [118, 3],
                [137, 3],
                [157, 3],
                [177, 3],
                [196, 3],
                [216, 4],
                [235, 4],
                [255, 4],
            ]),
            scaling_points_cb: ArrayVec::new(),
            scaling_points_cr: ArrayVec::new(),
            scaling_shift: 8,
            ar_coeff_lag: 0,
            ar_coeffs_y: ArrayVec::new(),
            ar_coeffs_cb: ArrayVec::try_from([0].as_slice()).expect("Arrayvec has capacity"),
            ar_coeffs_cr: ArrayVec::try_from([0].as_slice()).expect("Arrayvec has capacity"),
            ar_coeff_shift: 6,
            cb_mult: 0,
            cb_luma_mult: 0,
            cb_offset: 0,
            cr_mult: 0,
            cr_luma_mult: 0,
            cr_offset: 0,
            overlap_flag: true,
            chroma_scaling_from_luma: false,
            grain_scale_shift: 0,
            random_seed: 7391,
        };
        let output = parse_grain_table(input).expect("Test failed");
        assert_eq!(vec![expected], output);
    }

    #[test]
    fn parse_luma_chroma_table() {
        // This is the luma+chroma table format generated by
        // both aomenc's photon noise utility and by av1an.
        let input = r#"filmgrn1
E 0 9223372036854775807 1 7391 1
  p 0 6 0 8 0 1 128 192 256 128 192 256
  sY 14  0 0 20 4 39 3 59 3 78 3 98 3 118 4 137 4 157 4 177 4 196 4 216 5 235 5 255 5
  sCb 10 0 0 28 0 57 0 85 0 113 0 142 0 170 0 198 0 227 0 255 1
  sCr 10 0 0 28 0 57 0 85 0 113 0 142 0 170 0 198 0 227 0 255 1
  cY
  cCb 0
  cCr 0
"#;
        let expected = GrainTableSegment {
            start_time: 0,
            end_time: 9_223_372_036_854_775_807,
            scaling_points_y: ArrayVec::from([
                [0, 0],
                [20, 4],
                [39, 3],
                [59, 3],
                [78, 3],
                [98, 3],
                [118, 4],
                [137, 4],
                [157, 4],
                [177, 4],
                [196, 4],
                [216, 5],
                [235, 5],
                [255, 5],
            ]),
            scaling_points_cb: ArrayVec::from([
                [0, 0],
                [28, 0],
                [57, 0],
                [85, 0],
                [113, 0],
                [142, 0],
                [170, 0],
                [198, 0],
                [227, 0],
                [255, 1],
            ]),
            scaling_points_cr: ArrayVec::from([
                [0, 0],
                [28, 0],
                [57, 0],
                [85, 0],
                [113, 0],
                [142, 0],
                [170, 0],
                [198, 0],
                [227, 0],
                [255, 1],
            ]),
            scaling_shift: 8,
            ar_coeff_lag: 0,
            ar_coeffs_y: ArrayVec::new(),
            ar_coeffs_cb: ArrayVec::try_from([0].as_slice()).expect("Arrayvec has capacity"),
            ar_coeffs_cr: ArrayVec::try_from([0].as_slice()).expect("Arrayvec has capacity"),
            ar_coeff_shift: 6,
            cb_mult: 128,
            cb_luma_mult: 192,
            cb_offset: 256,
            cr_mult: 128,
            cr_luma_mult: 192,
            cr_offset: 256,
            overlap_flag: true,
            chroma_scaling_from_luma: false,
            grain_scale_shift: 0,
            random_seed: 7391,
        };
        let output = parse_grain_table(input).expect("Test failed");
        assert_eq!(vec![expected], output);
    }

    #[test]
    fn parse_complex_table() {
        let input = r#"filmgrn1
E 0 417083 1 7391 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 6  0 53 13 53 40 64 94 49 121 46 255 46
	sCb 2 0 14 255 13
	sCr 2 0 12 255 14
	cY 1 -4 1 4 8 3 -2 -6 9 14 -27 -25 -2 4 5 15 -80 94 28 -3 -2 6 -47 121
	cCb -3 1 -4 6 -1 2 -2 1 11 -10 -2 -16 -1 3 -2 -14 -26 65 19 -3 -5 2 -6 75 -1
	cCr 0 0 -4 8 -1 0 1 2 -1 -9 4 -7 -5 -2 -5 -14 0 45 18 3 -3 4 8 49 5
E 417083 7090416 1 0 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 4  0 46 40 54 108 39 255 38
	sCb 2 0 14 255 14
	sCr 2 0 12 255 14
	cY 1 -4 1 5 8 4 -2 -6 9 13 -28 -28 -5 5 5 13 -76 91 32 -1 -3 7 -50 124
	cCb -2 1 -3 3 -2 1 -1 2 8 -10 0 -12 -2 2 -1 -14 -20 61 18 -1 -4 -2 -1 70 -1
	cCr 0 0 -3 6 -1 -1 0 1 -2 -8 6 -4 -5 -2 -6 -12 4 41 17 4 -2 3 13 44 5
E 7090416 7507500 1 0 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 4  0 54 40 64 108 46 255 44
	sCb 2 0 14 255 13
	sCr 2 0 12 255 14
	cY 1 -4 2 3 7 3 -2 -6 9 14 -26 -25 -3 5 6 15 -81 95 27 -3 -3 5 -46 121
	cCb -2 1 -4 4 -2 1 -1 2 9 -12 3 -13 -1 2 -2 -16 -26 66 17 -2 -5 -1 1 73 0
	cCr 1 -1 -5 8 -1 -1 1 1 -3 -9 9 -5 -6 -2 -7 -14 1 44 17 3 -3 5 15 46 4
E 7507500 10010000 1 0 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 4  0 49 40 59 108 43 255 41
	sCb 2 0 14 255 14
	sCr 2 0 13 255 15
	cY 1 -4 0 6 8 3 -2 -5 8 14 -29 -26 -3 4 3 15 -76 92 29 -2 -3 8 -49 121
	cCb -3 0 -3 6 0 1 -2 1 10 -9 -4 -15 -1 2 -1 -13 -22 62 20 -3 -4 2 -7 73 -1
	cCr -1 0 -3 6 0 0 0 2 0 -9 2 -7 -5 -1 -4 -14 0 45 19 2 -2 3 7 50 4
E 10010000 13346666 1 0 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 6  0 33 27 39 40 53 54 55 108 52 255 52
	sCb 2 0 16 255 14
	sCr 2 0 11 255 12
	cY 1 -4 1 5 9 4 -2 -7 12 11 -27 -30 -5 5 6 10 -73 89 35 -1 -3 6 -49 124
	cCb -2 0 -2 1 -2 1 -2 0 9 -9 -2 -14 -1 2 0 -11 -26 65 18 -2 -4 -2 -8 75 -5
	cCr 0 0 -4 5 -2 0 1 3 -1 -9 6 -5 -5 -1 -6 -14 1 43 18 4 -3 3 13 49 3
E 13346666 16683333 1 0 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 6  0 36 27 42 40 58 54 60 108 57 255 57
	sCb 2 0 15 255 14
	sCr 4 0 11 40 17 94 13 255 13
	cY 1 -4 1 5 8 3 -2 -6 10 12 -27 -27 -4 4 5 12 -73 90 32 -2 -3 6 -47 121
	cCb -2 0 -3 4 -1 1 -2 0 10 -9 -2 -14 1 3 -1 -10 -24 62 16 -2 -4 0 -6 72 -7
	cCr 0 0 -3 6 -1 0 1 3 1 -9 3 -7 -5 -1 -5 -14 -2 46 19 2 -3 3 7 54 3
E 16683333 17100416 1 0 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 7  0 41 13 41 27 49 40 66 54 68 108 65 255 65
	sCb 2 0 18 255 14
	sCr 4 0 11 40 18 67 14 255 13
	cY 0 -3 1 4 7 3 -2 -5 7 13 -27 -23 -3 4 5 15 -79 94 26 -3 -2 5 -45 120
	cCb -1 -2 -1 1 0 0 -3 -2 12 -6 -3 -15 3 2 2 -8 -42 75 12 -3 -4 -2 -8 82 -3
	cCr 0 0 -5 7 -2 0 1 3 0 -11 6 -7 -5 -1 -6 -15 -5 48 18 2 -3 3 10 55 2
E 17100416 20020000 1 0 1
	p 3 7 0 11 0 1 128 192 256 128 192 256
	sY 6  0 37 27 44 40 61 54 63 108 60 255 60
	sCb 2 0 14 255 14
	sCr 4 0 11 40 18 94 13 255 13
	cY 1 -3 0 6 7 2 -1 -5 7 13 -28 -25 -2 3 3 13 -73 91 29 -2 -2 7 -47 119
	cCb -2 -1 -3 4 0 1 -2 -1 11 -7 -6 -15 1 2 -1 -9 -25 63 16 -3 -4 2 -11 73 -8
	cCr -1 1 -2 6 0 1 0 2 3 -9 -2 -10 -4 0 -3 -14 -6 50 20 0 -3 3 -1 59 3
E 20020000 9223372036854775807 1 0 1
	p 3 6 0 11 0 1 128 192 256 128 192 256
	sY 6  0 32 27 37 40 50 54 52 121 49 255 49
	sCb 4 0 21 40 23 81 17 255 15
	sCr 2 0 11 255 12
	cY 1 -3 1 2 5 3 -2 -6 8 6 -12 -18 -2 3 5 7 -42 44 21 -3 -1 4 -29 67
	cCb -1 0 1 0 -1 0 -1 0 5 -4 -3 -9 1 1 2 -4 -21 39 10 -2 -3 -2 -7 44 1
	cCr 1 0 -3 2 -3 -1 0 1 -1 -4 5 -2 -1 -1 -5 -6 3 20 10 4 -2 0 9 23 -1"#;
        let output = parse_grain_table(input);
        assert!(output.is_ok());
    }

    #[test]
    fn parse_rejects_missing_required_line() {
        let input = r#"filmgrn1
E 0 100 1 7391 1
  p 0 6 0 8 0 1 0 0 0 0 0 0
  sY 0
  sCb 0
  sCr 0
  cY
  cCb 0
"#;

        let err = parse_grain_table(input).expect_err("missing cCr line should fail");
        assert!(err.to_string().contains("Expected cCr line"));
    }

    #[test]
    fn parse_rejects_too_many_luma_scaling_points() {
        let input = r#"filmgrn1
E 0 100 1 7391 1
  p 0 6 0 8 0 1 0 0 0 0 0 0
  sY 15 0 0 1 0 2 0 3 0 4 0 5 0 6 0 7 0 8 0 9 0 10 0 11 0 12 0 13 0 14 0
  sCb 0
  sCr 0
  cY
  cCb 0
  cCr 0
"#;

        let err = parse_grain_table(input).expect_err("too many Y points should fail");
        assert!(
            err.to_string()
                .contains("expected at most 14 Y-plane points")
        );
    }
}

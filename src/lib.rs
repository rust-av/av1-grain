// Copyright (c) 2022-2022, The rav1e contributors. All rights reserved
//
// This source code is subject to the terms of the BSD 2 Clause License and
// the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
// was not distributed with this source code in the LICENSE file, you can
// obtain it at www.aomedia.org/license/software. If the Alliance for Open
// Media Patent License 1.0 was not distributed with this source code in the
// PATENTS file, you can obtain it at www.aomedia.org/license/patent.

#![deny(clippy::all)]
#![warn(clippy::nursery)]
#![warn(clippy::pedantic)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::default_trait_access)]
#![allow(clippy::inconsistent_struct_constructor)]
#![allow(clippy::inline_always)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::redundant_closure_for_method_calls)]
#![allow(clippy::similar_names)]
#![allow(clippy::struct_excessive_bools)]
#![allow(clippy::use_self)]
#![warn(clippy::clone_on_ref_ptr)]
#![warn(clippy::create_dir)]
#![warn(clippy::dbg_macro)]
#![warn(clippy::default_numeric_fallback)]
#![warn(clippy::exit)]
#![warn(clippy::filetype_is_file)]
#![warn(clippy::float_cmp_const)]
#![warn(clippy::if_then_some_else_none)]
#![warn(clippy::lossy_float_literal)]
#![warn(clippy::map_err_ignore)]
#![warn(clippy::mem_forget)]
#![warn(clippy::mod_module_files)]
#![warn(clippy::multiple_inherent_impl)]
#![warn(clippy::pattern_type_mismatch)]
#![warn(clippy::rc_buffer)]
#![warn(clippy::rc_mutex)]
#![warn(clippy::rest_pat_in_fully_bound_structs)]
#![warn(clippy::same_name_method)]
#![warn(clippy::str_to_string)]
#![warn(clippy::string_to_string)]
#![warn(clippy::undocumented_unsafe_blocks)]
#![warn(clippy::unnecessary_self_imports)]
#![warn(clippy::unneeded_field_pattern)]
#![warn(clippy::use_debug)]
#![warn(clippy::verbose_file_reads)]

#[cfg(feature = "create")]
mod create;
#[cfg(feature = "diff")]
mod diff;
#[cfg(all(feature = "estimate", feature = "unstable"))]
mod estimate;
#[cfg(feature = "parse")]
mod parse;
mod util;

use arrayvec::ArrayVec;
#[cfg(feature = "create")]
pub use create::*;
#[cfg(feature = "diff")]
pub use diff::*;
#[cfg(all(feature = "estimate", feature = "unstable"))]
pub use estimate::*;
#[cfg(feature = "parse")]
pub use parse::*;
pub use v_frame;

/// The max number of luma scaling points for grain synthesis
pub const NUM_Y_POINTS: usize = 14;
/// The max number of scaling points per chroma plane for grain synthesis
pub const NUM_UV_POINTS: usize = 10;
/// The max number of luma coefficients for grain synthesis
pub const NUM_Y_COEFFS: usize = 24;
/// The max number of coefficients per chroma plane for grain synthesis
pub const NUM_UV_COEFFS: usize = 25;

/// A randomly generated u16 to be used as a starting random seed
/// for grain synthesis. The idea behind using a constant random seed
/// is so that encodes are deterministic and reproducible.
pub const DEFAULT_GRAIN_SEED: u16 = 10956;

pub type ScalingPoints = ArrayVec<[u8; 2], NUM_Y_POINTS>;

/// Specifies parameters for enabling decoder-side grain synthesis for
/// a segment of video from `start_time` to `end_time`.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serialize", derive(serde::Deserialize, serde::Serialize))]
pub struct GrainTableSegment {
    /// The beginning timestamp of this segment, in 10,000,000ths of a second.
    pub start_time: u64,
    /// The ending timestamp of this segment, not inclusive, in 10,000,000ths of
    /// a second.
    pub end_time: u64,

    /// Values for the cutoffs and scale factors for luma scaling points
    pub scaling_points_y: ArrayVec<[u8; 2], NUM_Y_POINTS>,
    /// Values for the cutoffs and scale factors for Cb scaling points
    pub scaling_points_cb: ArrayVec<[u8; 2], NUM_UV_POINTS>,
    /// Values for the cutoffs and scale factors for Cr scaling points
    pub scaling_points_cr: ArrayVec<[u8; 2], NUM_UV_POINTS>,

    /// Determines the range and quantization step of the standard deviation
    /// of film grain.
    ///
    /// Accepts values between `8..=11`.
    pub scaling_shift: u8,

    /// A factor specifying how many AR coefficients are provided,
    /// based on the forumla `coeffs_len = (2 * ar_coeff_lag * (ar_coeff_lag +
    /// 1))`.
    ///
    /// Accepts values between `0..=3`.
    pub ar_coeff_lag: u8,
    /// Values for the AR coefficients for luma scaling points
    pub ar_coeffs_y: ArrayVec<i8, NUM_Y_COEFFS>,
    /// Values for the AR coefficients for Cb scaling points
    pub ar_coeffs_cb: ArrayVec<i8, NUM_UV_COEFFS>,
    /// Values for the AR coefficients for Cr scaling points
    pub ar_coeffs_cr: ArrayVec<i8, NUM_UV_COEFFS>,
    /// Shift value: Specifies the range of acceptable AR coefficients
    /// 6: [-2, 2)
    /// 7: [-1, 1)
    /// 8: [-0.5, 0.5)
    /// 9: [-0.25, 0.25)
    pub ar_coeff_shift: u8,
    /// Multiplier to the grain strength of the Cb plane
    pub cb_mult: u8,
    /// Multiplier to the grain strength of the Cb plane inherited from the luma
    /// plane
    pub cb_luma_mult: u8,
    /// A base value for the Cb plane grain
    pub cb_offset: u16,
    /// Multiplier to the grain strength of the Cr plane
    pub cr_mult: u8,
    /// Multiplier to the grain strength of the Cr plane inherited from the luma
    /// plane
    pub cr_luma_mult: u8,
    /// A base value for the Cr plane grain
    pub cr_offset: u16,

    /// Whether film grain blocks should overlap or not
    pub overlap_flag: bool,
    /// Scale chroma grain from luma instead of providing chroma scaling points
    pub chroma_scaling_from_luma: bool,
    /// Specifies how much the Gaussian random numbers should be scaled down
    /// during the grain synthesis process.
    pub grain_scale_shift: u8,
    /// Random seed used for generating grain
    pub random_seed: u16,
}

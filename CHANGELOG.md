# CHANGELOG

## Version 0.5.0

- [Breaking] Bump MSRV to Rust 1.95.
- [Breaking] Update `v_frame` dependency to 0.7.
- [Breaking] Update package licensing metadata to use the bundled `LICENSE` file, covering BSD-2-Clause and AOMPL-1.0 terms.
- Remove the `nom`, `cfg-if`, `quickcheck`, and `quickcheck_macros` dependencies.
- Remove default features from the optional `num-rational` dependency.
- Rewrite grain table parsing without `nom`.
- Fix noise estimation for padded `v_frame` planes.
- Remove unused legacy FFmpeg and VapourSynth diff code.
- Consolidate CI workflows into a single workflow and add `just precommit`.
- Simplify tests, clippy lint configuration, and formatting configuration.

## Version 0.4.2

- perf/safety: use get_dbg utility, eliminate current unchecked accesses
- fix: operator precedence in solver could lead to OOB access

## Version 0.4.1

- Further improvements to photon noise generation based on [https://github.com/juliobbv-p/svt-av1-hdr/pull/32](https://github.com/juliobbv-p/svt-av1-hdr/pull/32)

## Version 0.4.0

- Improve limited range photon noise generation.
  - [Breaking] Introduces a `full_range` field on `NoiseGenArgs`
- Bump to Rust edition 2024

## Version 0.3.0

- [Breaking] Update `v_frame` dependency to 0.5

## Version 0.2.5

- Bump `nom` dependency to 8.0

## Version 0.2.4

- Fix a compilation issue with `--no-default-features`

## Version 0.2.3

- Many speed optimizations to diff

## Version 0.2.2

- Fix issue where `NoiseModel` may fail in certain circumstances.
- Considerably speed up `NoiseModel` calculations.

## Version 0.2.1

- Bump `v_frame` to 0.3
- Fix a clippy warning

## Version 0.2.0

- [Breaking] Change the name of `generate_grain_params` to `generate_photon_noise_params`. This was done to support the future `generate_film_grain_params` feature.
- [Feature] Add the `diff` module which contains the `DiffGenerator` struct. This takes in a series of source frames and denoised frames and generates a grain table based on the difference. This feature is enabled by default.

## Version 0.1.4

- Fix a bug that prevented `generate_luma_noise_points` from generating any luma noise points.
  - ALL previous versions have been yanked because of the severity of this bug. Please update to this one.

## Version 0.1.3

- Be more consistent in using `anyhow::Result`

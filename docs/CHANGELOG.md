# Changelog

All notable changes to this project will be documented in this file.

## [4.0.0] - 2025-08-28

### Highlights (4.0.0)

- Finalized non-integer (beta) base support in all core functions.
- Vectorized multi-n exponent expansion via shared lattice (`ExpoExpand`).
- Combined multi-n grid rendering with aligned integer columns.
- High-resolution (256) colormap usage; removed low-base banding.
- Figure reuse & handle struct in `GridDigitVis` (naming/tagging for exports).
- `ColorScale` (linear/exponential) residual mapping; unified integer colorbar ticks.
- Simplified titles and consistent colorbar semantics.

### Migration (4.0.0)

- Replace legacy script patterns with function pipeline (see README Quick Start).
- Any direct `GridDigitViz` calls now routed through updated `GridDigitVis` (alias retained).

## [3.1.0] - 2025-08-28

### Added (3.1.0)

- `GridDigitVis`: Figure / handle management (`Figure`, `FigureName`, `Tag`) and second output struct with handles.
- `GridDigitVis`: `ColorScale` option (`linear` | `exponential`) for residual mode color distribution.
- Unified colorbar tick strategy (integer ticks 0..floor(baseRadix)-1) consistent across modes.

### Changed (3.1.0)

- Always use full-resolution (256) base colormap (parula) to avoid low-base banding; discrete digit colors derived from high-res map.
- Titles simplified (plain root symbol) and row count removed from multi-row grid title for clarity.
- Residual scaling uses base units (0..baseRadix) internally while displaying integer ticks.

### Fixed (3.1.0)

- Colormap truncation artifacts for small `baseRadix` (<3) by eliminating `parula(digitCount)` down-sampling.
- Off-by-one colorbar tick issue (extra tick at end) corrected.

## [3.0.0] - 2025-08-26

### Added (3.0.0)

- Modular `src/` folder with core functions.
- Batch exponent expansion (`ExpoExpand`).
- Grid digit image export (`GridDigitViz`).
- Dynamic video coloring & sliding window masking.
- Documentation (README, CITATION, MIT License).
- Basic unit tests.

### Changed (3.0.0)

- Refactored monolithic logic into testable functions.
- Generalized base handling integrated from earlier prototypes.

### Removed

- Implicit decimal-only assumptions.

## [2.0.0] - 2025 (Earlier)

### Added (2.0.0)

- Vectorized digit extraction across uniform exponent grid.
- Base generalization (parameterized baseRadix).
- Repetition guides & improved polar layout.

### Changed (2.0.0)

- Performance improvements over v1 loops.

## [1.0.0] - 2025 (Early Prototype)

### Added (1.0.0)

- Initial concatenated sequence visualization (decimal only).
- Per-layer looped digit extraction & video sweep.

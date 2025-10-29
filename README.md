# SchizoVis (v4.0)

**Author**: Maziar Moussavi
**Current Version**: 4.0.1 (2025-08-28)
**Languages**: MATLAB

## Summary

SchizoVis is a MATLAB project for visualizing the mathematical structure of the "schizophrenic number" sequence. This sequence is defined by the recurrence relation:

f(0) = 0, f(k) = b * f(k-1) + k

where `b` is a user-selectable `baseRadix` (>= 2). In base 10, this is equivalent to concatenating the integers 1, 2, 3, ... k. SchizoVis generalizes this to any base, including non-integer "beta" bases.

The project provides tools to:
*   Generate the sequence and its square root using symbolic math to maintain precision.
*   Expand the square root into its digit or residual representation in the chosen base.
*   Visualize the resulting digit patterns as stacked polar plots or grid images.

## What's New in v4.0

This major release focuses on generalized bases, performance, and unified multi-n visualization:

*   **Full Non-Integer (Beta) Base Support**: End-to-end support for non-integer bases in generation, expansion, and visualization.
*   **Vectorized Expansion**: The `ExpoExpand` function efficiently expands multiple square root values over a shared exponent lattice.
*   **Unified Multi-n Grid View**: A combined grid view that displays multiple `n` values as rows with aligned integer parts, available in both 'digits' and 'residual' modes.
*   **High-Resolution Colormap Pipeline**: Ensures stable and consistent color appearance, especially for small bases.
*   **Enhanced Figure Management**: `GridDigitVis` now offers better control over figures with options like `Figure`, `FigureName`, `Tag`, and returns a handle structure.
*   **Exponential Color Scaling**: A `ColorScale` option (`linear` or `exponential`) for emphasizing residuals.
*   **Unified Colorbar Ticks**: Integer colorbar ticks (0..floor(baseRadix)-1) across all modes and simplified plain-text titles.

## Repository Structure

```text
.
├── README.md                # Project overview and usage guide
├── LICENSE                  # Project license
├── .gitignore               # Files and directories to be ignored by Git
├── setup.m                  # MATLAB script to set up the project path
├── src/                     # Core library functions
│   ├── SchizoGen.m          # Generates the schizophrenic number sequence
│   ├── ExpoExpand.m         # Expands symbolic values into digits/residuals
│   ├── PolarDigitVis.m      # Creates stacked polar digit plots
│   └── GridDigitVis.m       # Creates grid-based digit visualizations
├── examples/                # Runnable example and driver scripts
│   └── SchizoViz_demo_v3.m
├── legacy/                  # Historical versions of the project
├── docs/                    # Additional documentation and changelog
├── outputs/                 # (Git-ignored) Default directory for figures and videos
└── tests/                   # MATLAB unit tests
```

## Quick Start

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/MahZadYar/SchizoVis.git
    ```
2.  **Run the setup script in MATLAB:**
    From the repository's root directory in MATLAB, run:
    ```matlab
    setup
    ```
    This will add the `src` directory to your MATLAB path.
3.  **Run the demo:**
    Open and run the main demonstration script to see SchizoVis in action:
    ```matlab
    examples/SchizoViz_demo_v3.m
    ```
4.  **Customize and Explore:**
    Adjust the parameters at the top of the demo script (e.g., `n` range, `baseRadix`, `precisionOrder`) to explore different visualizations. The generated images and videos will be saved in the `outputs` directory.

## Core Concepts

### 1. Sequence Construction (`SchizoGen`)

The `SchizoGen` function is the starting point of the workflow. It takes a vector of `n` values and a `baseRadix` and generates the corresponding f(n) and sqrt(f(n)) values. By using MATLAB's Symbolic Math Toolbox, it avoids numerical precision loss, which is crucial for the deep analysis of the digit patterns.

### 2. Exponent Lattice

To analyze the digits of the square root values, a shared grid of exponents is created. This "exponent lattice" is a descending vector of exponents that defines the positions (powers of the base) at which the digits or residuals are sampled. This ensures that the expansions for different `n` values are aligned and comparable.

### 3. Batch Expansion (`ExpoExpand`)

`ExpoExpand` is the computational core of the project. It takes the symbolic square root values from `SchizoGen` and expands them into a matrix of digits or residuals. It is highly vectorized to process multiple values in a batch, making it efficient. The `mode` parameter controls whether the output is discrete digits (`'digits'`) or continuous, normalized residuals (`'residual'`).

### 4. Visualization

SchizoVis offers two primary visualization modalities:

*   **Polar Mapping (`PolarDigitVis`)**: This function creates a "clock-like" radial plot. The radius corresponds to the exponent, and the angle corresponds to the digit value. When visualizing multiple `n` values, the plots are stacked in 3D. In residual mode, the color of each point represents the residual magnitude.

*   **Grid Digit Image (`GridDigitVis`)**: This function generates a rectangular image where each pixel corresponds to a digit. When multiple `n` values are provided, they are rendered as rows in the image, with their integer parts aligned. This is useful for observing how the digit patterns evolve as `n` increases.

## Function Overview

The `src` directory contains the core functions of SchizoVis:

*   `SchizoGen(n, baseRadix)`: Generates the symbolic f(n) and sqrt(f(n)) values.
*   `ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode)`: Expands the symbolic square root values into a matrix of digits or residuals.
*   `PolarDigitVis(exponents, digitsMat, baseRadix, varargin)`: Creates a stacked polar visualization of the digit matrix.
*   `GridDigitVis(symVal, baseRadix, outFile, nVal, varargin)`: Creates a grid-based image visualization of the digit matrix.

## Minimal Example

Here is a minimal script to generate and visualize a SchizoVis plot:

```matlab
% --- Parameters ---
nMin = 1;
nMax = 101;
nStepSize = 2;
baseRadix = 12;         % The numeral base for the sequence
precisionOrder = 1500;  % The depth of fractional exponents to sample
mode = "digits";        % "digits" or "residual"

% --- Generate the sequence ---
nList = nMin:nStepSize:nMax;
[fVals, sVals] = SchizoGen(nList, baseRadix);

% --- Create the exponent lattice ---
p_min = -abs(precisionOrder);
global_p_max = ceil(double(log(vpa(sVals(end)))/log(baseRadix)));
exponents = global_p_max:-1:p_min;

% --- Expand the values ---
digitsMat = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);

% --- Visualize the result ---
PolarDigitVis(exponents, digitsMat, baseRadix, 'Title', 'SchizoVis Sample');
```

## Parameter Notes

*   `nMin`, `nMax`, `nStepSize`: These parameters define the range of `n` values to be sampled. Larger `n` values require more memory and computation time.
*   `baseRadix`: The numeral base, `b`. This can be any number greater than 1, including non-integers. It controls the digit alphabet size and the angular segmentation in the polar plot.
*   `precisionOrder`: This determines the negative depth of the exponent lattice (i.e., how many fractional digits are computed). Larger values provide a more detailed view of the fractional part but increase memory usage.
*   `mode`: This can be `'digits'` or `'residual'`.
    *   `'digits'`: The expansion returns integer digits, resulting in discrete color values in the visualizations.
    *   `'residual'`: The expansion returns normalized residual values in the range [0, 1), resulting in smooth color gradients.

## Performance Tips

*   The size of `digitsMat` (L x M) can become very large with a high `precisionOrder` or a large number of `n` values. Monitor memory usage accordingly.
*   `ExpoExpand` includes a validation step that may issue warnings if precision mismatches are detected. If this occurs, you may need to adjust the precision parameters.
*   Video exports can generate very large files. For smaller file sizes, consider using a compressed format like MPEG-4 if lossless quality is not required.

## Testing

To run the suite of unit tests, execute the following in MATLAB:

```matlab
results = runtests('tests');
table(results)
```

## Contributing

Contributions are welcome! Please follow these steps:

1.  Open an issue to describe the bug or enhancement.
2.  Fork the repository and create a feature branch.
3.  Add or extend tests to cover your changes.
4.  Submit a pull request that references the issue.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

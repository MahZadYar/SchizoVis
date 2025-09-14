SchizoVis (v4.0)
=================

Author: Maziar Moussavi
Current Version: 4.0.1 (2025-08-28)
Languages: MATLAB (core algorithms + visualization)

## Purpose

SchizoVis is a MATLAB project for visualizing the mathematical structure of the "schizophrenic number" sequence. This sequence is defined by the recurrence relation:

  f(0) = 0,  f(k) = b * f(k-1) + k

where `b` is a user-selectable numeral base (`baseRadix`), which can be an integer or a non-integer value greater than 1. For base 10, this process is equivalent to concatenating the integers from 1 to k (e.g., f(3) = 123).

The primary goal of this project is to explore the intricate patterns that emerge from the square roots of this sequence, `sqrt(f(n))`. It extracts the digits (or continuous residuals) of these square roots across a specified range of exponents and visualizes them in two main ways:

1.  **Polar Visualization (`PolarDigitVis`)**: Maps digits to polar coordinates, with the digit value determining the angle and the exponent determining the radius. When visualizing multiple `n` values, these polar plots are stacked along the z-axis, creating a 3D structure.
2.  **Grid Visualization (`GridDigitVis`)**: Represents the digits as colored pixels in a grid, which can be arranged as a square for a single number or as a multi-row image for comparing multiple numbers.

The project is designed to be modular and extensible, separating the data generation, expansion, and visualization stages.

## Setup and Installation

Follow these steps to get the project running on your local machine.

### 1. Clone the Repository

First, clone the repository to your local machine using Git:

```bash
git clone https://github.com/MahZadYar/SchizoVis.git
cd SchizoVis
```

### 2. Add Source Files to MATLAB Path

Before you can run the examples or use the core functions, you need to add the `src` directory to your MATLAB path. A convenience script, `setup.m`, is provided for this purpose.

From the root directory of the project in your MATLAB session, run:

```matlab
>> setup
```

This command will add the `src` folder to the MATLAB path for the current session, making all the project's functions available. You only need to do this once per session.

## Usage

The main workflow of SchizoVis involves three steps:
1.  **Generate** the symbolic sequence values.
2.  **Expand** the symbolic values into a matrix of digits or residuals.
3.  **Visualize** the resulting matrix.

The `examples/SchizoViz_demo_v3.m` script provides a comprehensive demonstration of how to use the various functions and features of the project.

### Minimal Example

Here is a minimal example of how to use the core functions to generate and visualize a sequence:

```matlab
% 1. Define parameters
nMin = 1; nMax = 101; nStepSize = 2;
baseRadix = 12;
precisionOrder = 1500;
mode = "digits"; % or "residual"

% 2. Generate the symbolic sequence
nList = nMin:nStepSize:nMax;
[fVals, sVals] = SchizoGen(nList, baseRadix);

% 3. Define the exponent range for expansion
p_min = -abs(precisionOrder);
% Use vpa for high-precision calculation of the max exponent
s_max_vpa = vpa(sVals(end));
global_p_max = ceil(double(log(s_max_vpa) / log(baseRadix)));
exponents = global_p_max:-1:p_min;

% 4. Expand the symbolic values into a digit matrix
digitsMat = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);

% 5. Visualize the result
figure;
PolarDigitVis(exponents, digitsMat, baseRadix, 'Title', 'SchizoVis Sample');
```

### Core Functions

*   `SchizoGen(n, baseRadix)`: Generates the symbolic `f(n)` and `sqrt(f(n))` values.
*   `ExpoExpand(sVals, exponents, ...)`: Performs the base-b expansion on the symbolic values.
*   `PolarDigitVis(exponents, digitsMat, ...)`: Creates the stacked polar plot.
*   `GridDigitVis(sVals, baseRadix, ...)`: Creates the grid image visualization.

Each function is thoroughly documented. Use the `help` command in MATLAB (e.g., `help SchizoGen`) to see detailed information about each function's parameters and options.

## Repository Structure

```text
README.md                Project overview & usage
LICENSE                  License (MIT by default – adjust if needed)
.gitignore               Ignore rules (figures, videos, temp, artifacts)
setup.m                  Convenience path setup helper
src/                     Core library functions
  SchizoGen.m
  ExpoExpand.m
  PolarDigitVis.m
  GridDigitVis.m
examples/                Runnable example / driver scripts
  SchizoViz_demo_v3.m
legacy/                  Historical versions
  SchizoViz-v1.m
  SchizoViz-v2.m
docs/                    Additional documentation & change log
  CHANGELOG.md
outputs/                 (gitignored) Figures, videos, exported images
tests/                   MATLAB unit tests
  testSchizoGen.m
  testExpoExpand.m
```

## Testing

To run the test suite, execute the following command in the MATLAB command window from the project's root directory:

```matlab
results = runtests('tests');
table(results)
```

## Contributing

Contributions are welcome! Please follow these steps:

1.  Open an issue to describe the enhancement or bug.
2.  Fork the repository and create a feature branch.
3.  Add or extend tests for your changes.
4.  Submit a pull request referencing the issue.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

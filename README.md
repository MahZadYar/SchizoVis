SchizoVis (v3)
=================

Author: Maziar Moussavi  
Current Version: 3 (2025-08)  
Languages: MATLAB (core algorithms + visualization)

Summary
-------
SchizoVis (renamed from SchizoViz) visualizes the evolving square root of the "schizophrenic number" sequence:

  f(0) = 0,  f(k) = b * f(k-1) + k

Where b = baseRadix (user selectable, b >= 2). For base 10 this is ordinary decimal concatenation of 1..k. For other bases it generalizes the process into base-b digit append operations. The project extracts layered digit (or residual) structure of sqrt(f(n)) across a controlled exponent window, mapping digits to polar coordinates and stacking layers (n values) in Z. It also supports exporting progressive videos and grid digit images.

What's New in v3
----------------
v3 refactors the exploratory v1/v2 scripts into modular, testable components:

* Modular functions in `src/`: `SchizoGen`, `ExpoExpand`, `PolarDigitVis`, `GridDigitViz`.
* Efficient vectorized expansion of multiple sqrt values across a shared exponent lattice (`ExpoExpand`).
* Generalized numeral base (any b >= 2) instead of fixed decimal logic.
* Two visualization modalities: stacked polar digit/residual plot and square grid digit image export.
* Optional dynamic coloring & sliding window video rendering.
* Cleaner separation between generation, expansion, and rendering stages.
* Basic MATLAB unit tests under `tests/`.

Repository Structure
--------------------

```text
README.md                Project overview & usage
LICENSE                  License (MIT by default – adjust if needed)
.gitignore               Ignore rules (figures, videos, temp, artifacts)
setup.m                  Convenience path setup helper
src/                     Core library functions
  SchizoGen.m
  ExpoExpand.m
  PolarDigitVis.m
  GridDigitViz.m
examples/                Runnable example / driver scripts
  SchizoViz_demo_v3.m
legacy/                  Historical versions (unaltered except header notes)
  SchizoViz-v1.m
  SchizoViz-v2.m
docs/                    Additional documentation & change log
  CHANGELOG.md
outputs/                 (gitignored) Figures, videos, exported images
tests/                   MATLAB unit tests (run with runtests)
  testSchizoGen.m
  testExpoExpand.m
```

Quick Start
-----------

1. Clone the repository:

  ```bash
  git clone https://github.com/YOUR_GITHUB_USERNAME/SchizoVis.git
  ```
2. In MATLAB, from the repo root run:

  ```matlab
  setup
  ```
3. Open and run the example:

  ```matlab
  examples/SchizoViz_demo_v3.m
  ```
4. Adjust parameters near the top (n range, baseRadix, precisionOrder, export toggles).

Core Concepts
-------------

1. Sequence Construction (`SchizoGen`): Builds symbolic f(n) and sqrt(f(n)) without premature numeric rounding.
2. Exponent Lattice: A shared descending exponent vector is constructed from the largest sqrt value; defines digit sampling positions.
3. Batch Expansion (`ExpoExpand`): Vectorized extraction of either discrete digits or residual values for all layers.
4. Polar Mapping (`PolarDigitVis`): Maps each digit to an angular slot (0..baseRadix-1) and exponent to radial coordinate; layers stack along Z.
5. Grid Digit Image (`GridDigitViz`): Renders a single sqrt(f(n)) into a square (or rectangular) colored digit raster (integer + fractional digits).

Example (Minimal)
-----------------

```matlab
nMin = 1; nMax = 101; nStepSize = 2;
baseRadix = 12; precisionOrder = 1500; mode = "digits"; % or "residual"
nList = nMin:nStepSize:nMax;
[fVals, sVals] = SchizoGen(nList, baseRadix);
p_min = -abs(precisionOrder);
global_p_max = ceil(double(log(vpa(sVals(end)))/log(baseRadix)));
exponents = global_p_max:-1:p_min;
digitsMat = ExpoExpand(sVals, exponents, precisionOrder, baseRadix, mode);
PolarDigitVis(exponents, digitsMat, baseRadix, 'Title', 'SchizoVis Sample');
```

Parameter Notes
---------------

* nMin, nMax, nStepSize: Define sampled n values.
* baseRadix (b): Controls the digit alphabet size and angular segmentation.
* precisionOrder: The negative depth (|p_min|) of fractional exponent sampling.
* mode: "digits" uses floor at each position; "residual" carries full residual values for smooth color/position (less common here).
* exportFigure / exportVideo / exportNumber (in example script) toggle resource-heavy exports.

Performance Tips
----------------

* Large precisionOrder or large nMax -> high memory; monitor `digitsMat` size (L x M).
* `ExpoExpand` auto-inflates decimal precision with guard digits; if validation reveals mismatches, increase guard or reduce precision scope.
* Video exports (lossless AVI) can become very large; switch to MPEG-4 when acceptable.

Testing
-------

Run all tests:

```matlab
results = runtests('tests');
table(results)
```

Add additional tests following MATLAB's `matlab.unittest` patterns.

Outputs
-------

All exported PNG/AVI/MP4 files should be directed to `outputs/` (default in example script uses an absolute path; change it to `fullfile(pwd,'outputs')` for portability).

Version History (Summary)
-------------------------

See `docs/CHANGELOG.md` for detailed entries.

* v1: Initial monolithic script (decimal only, per-layer loops, early 3D scatter & video sweep).
* v2: Vectorized digit extraction, generalized base parameter, improved guide overlays.
* v3: Modularization, batch exponent expansion, grid digit image export, dynamic video coloring, repository hygiene.

Planned / Ideas
---------------

* Optional GPU acceleration (prototype flag in v2; not active in v3).
* Interactive UI (App Designer) for parameter sweeping.
* Export to WebGL (e.g., Plotly) for browser sharing.
* Python translation (SymPy + Matplotlib / Plotly) for broader accessibility.

Citation
--------

If you use SchizoVis (formerly SchizoViz) in research or publications, please cite (CITATION.cff provided):

Moussavi, M. (2025). SchizoVis (Version 3). GitHub repository: [https://github.com/YOUR_GITHUB_USERNAME/SchizoVis](https://github.com/YOUR_GITHUB_USERNAME/SchizoVis)

License
-------

MIT License (modifiable). See `LICENSE`.

Contributing
------------

1. Open an issue describing enhancement / bug.
2. Fork and create a feature branch.
3. Add/extend tests for changes.
4. Submit PR referencing the issue.

Acknowledgements
----------------

Developed by Maziar Moussavi. Assisted organization & documentation support by AI pairing (GitHub Copilot style).

---
Enjoy exploring digit structures across numeral bases!

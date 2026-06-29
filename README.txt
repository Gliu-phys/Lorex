LOREX: LORENTZIAN EXTRACTOR FOR MATLAB
======================================

Lorex.m is a robust MATLAB pipeline designed to automate the fitting and extraction of overlapping Lorentzian peaks in 1D spectra and/or 2D datasets.

Unlike standard curve-fitting tools that struggle with multi-peak arrays, Lorex.m uses a dynamic 2nd-derivative 5-point stencil windowing technique to isolate peak cores, automatically rejects noise, and utilizes a constrained multi-peak refitting algorithm to untangle overlapping tails.

FEATURES:
- Automated Peak Finding: No need to manually click or guess peak locations.
- Overlap Handling: Drops degenerate peaks and refits residuals on the fly.
- Dynamic Bounding: Automatically calculates physical constraints based on the user's FWHM guess.
- Visual Diagnostics: Built-in visualization tools for step-by-step peak rejection and fit analysis.

PREREQUISITES:
- MATLAB R2018b or newer
- MATLAB Optimization Toolbox
- MATLAB Signal Processing Toolbox

USAGE EXAMPLE:

% 1. Load your raw data (must contain frequency vector 'x' and amplitude data 'y')
load('sampledata.mat'); 

% 2. Run Lorex.m with custom parameters
[y_rec, IntermediateResults, FinalResults] = lorex(frequencydata, spectra, ...
    'FWHM_guess',    0.5, ...                    % Guessed peak width in Hz
    'MinProminence', 0.01, ...                   % 1% threshold for peak finding
    'DevTol',        0.5, ...                    % Allowed peak drift in Hz
    'ResMult',       100, ...                    % Output resolution is 100x denser than input
    'CheckPlots',    true, ...                   % Show step-by-step plots
    'ShowHeatmap',   true);                      % Show 2D before-and-after heatmaps


CITATION:
If you use this code in your research, please cite it via Zenodo:

Liu, G. (2026). Lorex: Lorentzian Extractor for MATLAB (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.21019679


ACKNOWLEDGMENTS:
This function utilizes the highly robust lorentzfit.m algorithm developed by Jered R. Wells, retrieved from the MATLAB Central File Exchange. The lorentzfit code is included within this package under the terms of the BSD License.

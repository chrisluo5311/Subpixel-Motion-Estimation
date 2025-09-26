# Subpixel Motion Estimation Without Interpolation

## Introduction
This project reproduces a hybrid subpixel motion estimation pipeline for motion deblurring that combines coarse block matching with a first-order Taylor refinement so that fractional-pixel accuracy can be reached without explicitly interpolating intermediate frames. The accompanying experiments compare the hybrid strategy against classical full-search and logarithmic-search baselines to highlight the trade-offs between reconstruction quality, computational cost, and suitability for real-time restoration workflows.

## Features
- End-to-end demo that loads consecutive frames, estimates bidirectional motion, reconstructs a compensated frame, and reports MSE/PSNR together with visualizations of the motion vectors and residuals.
- Bidirectional motion estimation wrapper that fuses forward and backward block matching results to obtain robust motion-vector fields.
- Modular block-matching core with interchangeable search strategies (full search or logarithmic search) and the choice between Taylor-series refinement or bilinear interpolation for subpixel accuracy.
- Motion-compensation routine that upsamples the reference frame, warps pixels according to the estimated motion field, and resizes the result back to the original resolution for evaluation or deblurring.
- Optional motion-blur experiment that synthesizes blurred frames, estimates a PSF from the recovered motion vectors, and applies regularized deconvolution to illustrate downstream image restoration.

## Main Dependencies
- MATLAB or Octave-compatible environment with core image-processing functions such as `im2double`, `imread`, and `imshow` for data handling and visualization.
- Image Processing Toolbox (or equivalent) providing utilities like `fspecial`, `imfilter`, `imgaussfilt`, and `deconvreg` for blur synthesis, smoothing, and deconvolution experiments.

## Project Report
- [Subpixel Motion Estimation Without Interpolation with Two Blocking Methods](./Subpixel%20Motion%20Estimation%20Without%20Interpolation%20with%20Two%20Blocking%20Methods.pdf)

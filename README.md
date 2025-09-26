# Subpixel Motion Estimation Without Interpolation

## Introduction
This project reproduces a hybrid subpixel motion estimation pipeline for motion deblurring that combines coarse block matching with a first-order Taylor refinement so that fractional-pixel accuracy can be reached without explicitly interpolating intermediate frames. The accompanying experiments compare the hybrid strategy against classical full-search and logarithmic-search baselines to highlight the trade-offs between reconstruction quality, computational cost, and suitability for real-time restoration workflows.

## Features
- End-to-end demo that loads consecutive frames, estimates bidirectional motion, reconstructs a compensated frame, and reports MSE/PSNR together with visualizations of the motion vectors and residuals.【F:Demo.m†L2-L77】
- Bidirectional motion estimation wrapper that fuses forward and backward block matching results to obtain robust motion-vector fields.【F:Bidirectional_ME.m†L1-L40】
- Modular block-matching core with interchangeable search strategies (full search or logarithmic search) and the choice between Taylor-series refinement or bilinear interpolation for subpixel accuracy.【F:Motion_Est.m†L1-L110】【F:FullSearch.m†L1-L132】【F:LogSearch.m†L1-L178】【F:FullSearch_interpolation.m†L1-L113】【F:LogSearch_interpolation.m†L1-L160】
- Motion-compensation routine that upsamples the reference frame, warps pixels according to the estimated motion field, and resizes the result back to the original resolution for evaluation or deblurring.【F:reconstruct.m†L1-L83】
- Optional motion-blur experiment that synthesizes blurred frames, estimates a PSF from the recovered motion vectors, and applies regularized deconvolution to illustrate downstream image restoration.【F:deblurring_example.m†L1-L106】

## Main Dependencies
- MATLAB or Octave-compatible environment with core image-processing functions such as `im2double`, `imread`, and `imshow` for data handling and visualization.【F:Demo.m†L13-L77】
- Image Processing Toolbox (or equivalent) providing utilities like `fspecial`, `imfilter`, `imgaussfilt`, and `deconvreg` for blur synthesis, smoothing, and deconvolution experiments.【F:deblurring_example.m†L8-L53】

## Project Report
- [Subpixel Motion Estimation Without Interpolation with Two Blocking Methods](./Subpixel%20Motion%20Estimation%20Without%20Interpolation%20with%20Two%20Blocking%20Methods.pdf)

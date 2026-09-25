# Quantitative Signal-Processing Analysis of Color Doppler Twinkling

## Overview

This project investigates color Doppler twinkling across the ultrasound signal-processing chain using PMMA, low-twinkling PMMA, and metal targets.

## Processing Workflow

1. Channel RF acquisition
2. Beamforming
3. IQ demodulation
4. Wall filtering
5. Power Doppler calculation
6. x1 and x2 feature extraction
7. 3D joint feature-space analysis
8. Pulse-wise and frame-wise temporal analysis

## MATLAB Scripts

- `Beamforming_Twinkling_paper_code_batch_L7_4_final.m`  
  Beamforming, IQ conversion, wall filtering, Power Doppler calculation, and ROI extraction for L7-4 data.

- `Beamforming_Twinkling_paper_code_batch_L7_4_final.m`  
  Same processing workflow for L11-4v data.

- `Paramter_extraction_x1_x2_verification_Figure_for_paper_IUS.m`  
  Calculates x1, x2, log10(PD), BC3D, and generates 2D and 3D feature-space plots.

- `PostPorcess_Twinkle_Post_processing_Paper.m`  
  Generates RF, IQ, and Power Doppler summary plots and temporal variability results.

- `PostProcess_BF_RF_Variability_ACF_4panel.m`  
  Performs pulse-wise and frame-wise beamformed RF analysis and autocorrelation.

## Figures

![Signal-processing framework](figures/Picture5.svg)

![3D feature space](figures/Picture6.svg)


## Requirements

- MATLAB
- MATLAB Signal Processing Toolbox (MUST)
- MUST ultrasound toolbox

## Notes

The current scripts contain local file paths that must be updated before running on another computer.

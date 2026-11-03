# Radiomics Heterogeneity Module for 3D Slicer

A 3D volumetric radiomics module that computes heterogeneity features from
CT or PET volumes using a user-provided segmentation.

## Feature categories

| Category | Features |
|----------|----------|
| **First-order** (19) | `firstorder_Count`, `firstorder_Mean`, `firstorder_StdDev`, `firstorder_Variance`, `firstorder_Median`, `firstorder_Minimum`, `firstorder_Maximum`, `firstorder_Range`, `firstorder_InterquartileRange`, `firstorder_P5`, `firstorder_P25`, `firstorder_P75`, `firstorder_P95`, `firstorder_MeanAbsoluteDeviation`, `firstorder_Skewness`, `firstorder_Kurtosis`, `firstorder_Entropy`, `firstorder_CoV`, `firstorder_RobustSpreadEstimator` |
| **3D GLCM** (8) | `glcm_Contrast`, `glcm_Correlation`, `glcm_Energy`, `glcm_Homogeneity`, `glcm_Entropy`, `glcm_IdmVariance`, `glcm_InformationMeasureCorrelation1`, `glcm_InformationMeasureCorrelation2` |
| **3D GLRLM** (10) | `glrlm_RunLengthEmphasis`, `glrlm_ShortRunEmphasis`, `glrlm_LongRunEmphasis`, `glrlm_GrayLevelNonUniformity`, `glrlm_RunLengthNonUniformity`, `glrlm_RunLengthNonUniformityNormalized`, `glrlm_ShortRunLargeGrayEmphasis`, `glrlm_LongRunLargeGrayEmphasis`, `glrlm_RunLengthGrayLevelNonUniformity`, `glrlm_GrayLevelVariance` |
| **3D GLZDM** (5) | `gldm_SmallDependenceEmphasis`, `gldm_LargeDependenceEmphasis`, `gldm_DependencePercent`, `gldm_DependenceVariance`, `gldm_DependenceEntropy` |
| **NGTDM** (5) | `ngtdm_Coarseness`, `ngtdm_Contrast`, `ngtdm_Complexity`, `ngtdm_Strength`, `ngtdm_Granularity` |
| **Morphometry** (7) | `morphometry_VoxelCount`, `morphometry_Volume_mm3`, `morphometry_BoundingBox_X_mm`, `morphometry_BoundingBox_Y_mm`, `morphometry_BoundingBox_Z_mm`, `morphometry_Compactness`, `morphometry_Sphericity` |

**Total: 54 features** following the IBSI (Imaging Biomarker Standardization Initiative) naming convention.

## Installation

1. Copy the `RadiomicsHeterogeneity/` folder to your Slicer extensions directory:

    ```
    ~/.config/NA-MIC/Slicer/Extensions/lib/Slicer-5.4/RadiomicsHeterogeneity/
    ```

    (Adjust the Slicer version number as needed — e.g., `Slicer-5.7`.)

    Alternatively, place it at:

    ```
    <Slicer-install>/lib/Slicer-X.Y/RadiomicsHeterogeneity/
    ```

2. Restart Slicer.

3. The module appears under **Quantification → Radiomics Heterogeneity**.

## Usage

1. Load a CT or PET volume into the scene.
2. Create or load a segmentation with at least one segment (use the **Segment Editor** module).
3. Open **Quantification → Radiomics Heterogeneity**.
4. Select the **Volume** and **Segmentation** from the drop-down lists.
5. Select the **Gray Levels** (number of bins for intensity quantization) from the
   drop-down: **8, 16, 32, or 64**. Lower values speed up computation and
   reduce noise sensitivity.
6. Click **Compute Features**.
7. Results appear in the **Results** table. Click **Export CSV** to save.

## How it works

- The segment is converted to a binary labelmap in the volume's geometry.
- Only voxels inside the mask are used for computation.
- All texture features (GLCM, GLRLM, GLZDM) are computed in **3D** over all
  unique neighbour directions (7 positive-shift directions) and averaged for
  rotation-invariant results.
- Intensities are quantized to *N* gray levels (user-selected: 8, 16, 32, or 64)
  before texture computation.

## Dependencies

- 3D Slicer (>= 5.4)
- NumPy
- SciPy (for skewness/kurtosis)

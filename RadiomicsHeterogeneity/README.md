# Radiomics Heterogeneity Module for 3D Slicer

A 3D volumetric radiomics module that computes heterogeneity features from
CT or PET volumes inside a selected segment of a segmentation.

## Feature categories (87 features)

| Category | Features |
|----------|----------|
| **First-order** (25) | VoxelCount, Mean, StdDev, Variance, Median, Minimum, Maximum, Range, P5, P10, P25, P75, P90, P95, InterquartileRange, MeanAbsoluteDeviation, RobustMeanAbsoluteDeviation, MedianAbsoluteDeviation, RootMeanSquared, Skewness, Kurtosis (excess), Entropy, Uniformity, CoefficientOfVariation, QuartileCoefficientOfDispersion |
| **Shape** (9) | VoxelCount, VoxelVolume_mm3, MeshVolume_mm3, SurfaceArea_mm2, SurfaceVolumeRatio, Sphericity, BoundingBox_X/Y/Z_mm |
| **GLCM** (16) | JointEnergy, JointEntropy, Contrast, Correlation, Idm, Id, Idmn, Idn, InverseVariance, DifferenceEntropy, SumEntropy, ClusterShade, ClusterProminence, Imc1, Imc2, Autocorrelation |
| **GLRLM** (16) | Short/Long Run Emphasis, Gray Level / Run Length Non-Uniformity (+ normalized), Run Percentage, Gray Level / Run Variance, Run Entropy, Low/High Gray Level Emphasis, Short/Long Run Low/High Gray Level Emphasis |
| **GLSZM** (16) | Same set as GLRLM, computed on 26-connected zones |
| **NGTDM** (5) | Coarseness, Contrast, Busyness, Complexity, Strength |

Feature definitions follow IBSI / PyRadiomics. Values will not be identical to
PyRadiomics because no resampling is performed and GLCM/GLRLM matrices are
merged over directions (IBSI "3D merged") rather than averaged per direction.

## Installation

1. In Slicer: *Edit → Application Settings → Modules → Additional module paths*,
   add this `RadiomicsHeterogeneity/` folder.
2. Restart Slicer.
3. The module appears under **Quantification → Radiomics Heterogeneity**.

## Usage

1. Load a CT or PET volume and a segmentation (or create one in **Segment Editor**).
2. Open **Quantification → Radiomics Heterogeneity**.
3. Select the **Volume**, then the **Segmentation** and **Segment**.
4. Select **Gray levels** (8–128) for texture discretisation.
5. Click **Compute Features**. Results (with volume/segment names) appear in
   the table; **Export CSV** saves them.

## How it works

- The segment is converted to a binary labelmap in the volume's geometry,
  and arrays are cropped to the segment bounding box.
- First-order features use raw intensities (e.g. HU or SUV/Bq/mL).
- Texture features use fixed-bin-number discretisation of the ROI intensity
  range into *N* levels (1..N).
- GLCM and GLRLM use all 13 unique 3D directions (distance 1); GLCM is symmetric.
- GLSZM zones and the NGTDM neighbourhood use 26-connectivity; neighbours
  outside the segment are ignored.
- Surface area / mesh volume are computed with VTK marching cubes.
- Spacing is used for shape features only; for anisotropic voxels texture
  features are computed on the voxel grid.

## Testing

Developer mode → *Reload and Test* runs a self-test that checks the selected
segment (not the first one) is used.

## Dependencies

- 3D Slicer (>= 5.4) — NumPy, SciPy and VTK are bundled with Slicer.

# 3DSlicer-modules

A collection of scripted [3D Slicer](https://www.slicer.org/) modules for preclinical PET/CT image processing and analysis, developed at [@DobruckiLab](https://github.com/DobruckiLab).

All modules are pure-Python scripted modules — no compilation required. Install any of them by adding the module folder to **Edit → Application Settings → Modules → Additional module paths** and restarting Slicer.

---

## Modules

### [DynamicPETQuantifier](DynamicPETQuantifier/)

Quantifies **every frame** of a dynamic PET sequence against a pre-defined segmentation and exports all results to a single CSV or Excel file — replacing the manual frame-by-frame workflow.

- Iterates automatically over all frames in a Sequence Browser
- Computes voxel count, volume, min, max, mean, std dev, and median per segment per frame
- Supports uniform frame duration for a physical time axis
- Outputs a structured CSV or styled Excel file

**Module location in Slicer:** `Modules → Quantification → Dynamic PET Quantifier`

---

### [PADAnalysis](PADAnalysis/)

Splits a manually-drawn whole-hindlimb CT segmentation (left and/or right leg) into three anatomical zones — **Proximal**, **Distal**, and **Foot** — using cut-off planes positioned interactively on the CT (hip/upper-thigh boundary, knee joint, ankle joint), then quantifies PET tracer uptake in each zone.

- Interactive Markups plane widgets for the three per-leg cut-off boundaries, dragged directly onto the CT image
- Voxel-level classification against the planes (numpy), independent of plane orientation
- Generates six zone segments (`Left_Proximal`, `Left_Distal`, `Left_Foot`, `Right_Proximal`, `Right_Distal`, `Right_Foot`) in an output segmentation
- Quantifies PET uptake per zone via SegmentStatistics and exports a structured CSV

**Module location in Slicer:** `Modules → Quantification → PAD Analysis`

---

### [RadiomicsHeterogeneity](RadiomicsHeterogeneity/)

Extracts **54 IBSI-standardised 3D radiomics features** from a CT or PET volume using a segmentation, capturing tumour heterogeneity through first-order statistics, co-occurrence, run-length, grey-level zone, and neighbourhood grey-level difference matrices — all computed in 3D.

- **First-order** (19): mean, std dev, skewness, kurtosis, entropy, percentiles, and more
- **3D GLCM** (8): contrast, correlation, energy, homogeneity, entropy, and information measures
- **3D GLRLM** (10): run-length emphasis, short/long run emphasis, grey-level non-uniformity
- **3D GLZDM** (5): dependence emphasis, dependence percent, dependence entropy
- **NGTDM** (5): coarseness, contrast, complexity, strength, granularity
- **Morphometry** (7): voxel count, volume, bounding box, compactness, sphericity
- User-selectable grey-level quantisation (8, 16, 32, or 64 bins) for texture features
- Results table with one-click CSV export

**Module location in Slicer:** `Modules → Quantification → Radiomics Heterogeneity`

---

### [SlicerWorkflow](SlicerWorkflow/)

A pipeline orchestrator that lets you build, save, and run **multi-step image-analysis workflows** with pre-configured parameters — on a single volume or an entire folder of files.

Built-in steps include registration (BRAINSFit), N4 bias correction, Gaussian and median filtering, threshold and Otsu segmentation, volume/intensity statistics, and save volume. A Custom Python step allows arbitrary code injection.

- Drag-and-drop pipeline builder with typed parameter widgets
- Save/load/share pipelines as JSON files
- Batch-folder mode processes every matching file in a directory
- Extensible via a `STEP_CATALOGUE` dictionary

**Module location in Slicer:** `Modules → Utilities → Workflow`

---

### [SiemensImporter](SiemensImporter/)

Imports **Siemens Inveon** preclinical PET/CT datasets (`.img` / `.img.hdr`) and converts them to NIfTI format, with full support for multi-frame dynamic acquisitions.

- Parses Inveon ASCII headers and binary image data
- Handles PET (float32) and CT (int16 / Hounsfield Units) modalities
- Applies per-frame scale factors and optional scanner calibration
- Converts from LPS scanner coordinates to Slicer RAS automatically

**Module location in Slicer:** `Modules → Utilities → Siemens Inveon Importer`

---

### [SupertoroidalGlyphs](SupertoroidalGlyphs/)

Generates **supertoroidal glyph fields** from co-registered imaging data (e.g. CT + PET) and displays them as parameter-mapped 3D glyphs overlaid on a reference volume — visualising multi-parametric relationships in anatomical context.

- Map up to four co-registered volumes to glyph parameters: **color**, **N1**, **N2**, and **size**
- Configurable binning factor (NxNxN voxel blocks averaged into one glyph) for performance
- Each glyph's shape and appearance encodes the underlying multi-parametric data
- Export glyph table as CSV with spatial coordinates and all parameter values

**Module location in Slicer:** `Modules → Informatics → Supertoroidal Glyphs`

---

## Requirements

- 3D Slicer 5.0 or later
- Bundled Slicer Python packages: `SimpleITK`, `vtk`, `qt`, `ctk`, `numpy`
- `openpyxl` — required for Excel output in DynamicPETQuantifier (bundled in Slicer 5.x)
- `BRAINSFit` extension — required for the Registration step in SlicerWorkflow

---

## License

MIT License — free to use, modify, and distribute.

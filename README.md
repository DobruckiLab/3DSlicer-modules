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

### [SiemensImporter](SiemensImporter/)

Imports **Siemens Inveon** preclinical PET/CT datasets (`.img` / `.img.hdr`) and converts them to NIfTI format, with full support for multi-frame dynamic acquisitions.

- Parses Inveon ASCII headers and binary image data
- Handles PET (float32) and CT (int16 / Hounsfield Units) modalities
- Applies per-frame scale factors and optional scanner calibration
- Converts from LPS scanner coordinates to Slicer RAS automatically

**Module location in Slicer:** `Modules → Utilities → Siemens Inveon Importer`

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

## Requirements

- 3D Slicer 5.0 or later
- Bundled Slicer Python packages: `SimpleITK`, `vtk`, `qt`, `ctk`, `numpy`
- `openpyxl` — required for Excel output in DynamicPETQuantifier (bundled in Slicer 5.x)
- `BRAINSFit` extension — required for the Registration step in SlicerWorkflow

---

## License

MIT License — free to use, modify, and distribute.

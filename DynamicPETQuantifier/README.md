# Dynamic PET Quantifier — 3D Slicer Module

A scripted 3D Slicer module that quantifies **every frame** of a dynamic PET
sequence against a pre-defined segmentation and exports all results to a
**single CSV or Excel file** — replacing the manual frame-by-frame workflow.

---

## Features

| Feature | Details |
|---------|---------|
| Batch processing | Iterates automatically over all frames in a Sequence Browser |
| Flexible statistics | Voxel count, volume (mm³), min, max, mean, std dev, median — each optional |
| Time axis | Uses frame index, or computes mid-frame time from a uniform frame duration |
| Output formats | CSV (always available) or formatted Excel with header styling & frozen row |
| Non-destructive | Restores the original selected frame after running |

---

## Installation (no build required)

Because this is a **pure-Python scripted module** you can load it directly
without compiling anything.

1. Copy the file `DynamicPETQuantifier/DynamicPETQuantifier.py` anywhere on
   your computer (e.g. `~/SlicerModules/DynamicPETQuantifier.py`).
2. Open 3D Slicer → **Edit → Application Settings → Modules**.
3. Under *Additional module paths*, click **Add** and select the **folder**
   containing the `.py` file.
4. Click **OK** and restart Slicer when prompted.
5. The module appears under **Modules → Quantification → Dynamic PET Quantifier**.

---

## Workflow

### 1 — Prepare your scene

1. Load all dynamic PET frames (NIfTI or DICOM) into Slicer as individual
   scalar volumes.
2. Open the **Sequences** module and create a Sequence Browser:
   - Click *Add sequence* for each PET frame in time order.
   - The resulting *proxy volume* is the scalar volume node that updates as
     you move through frames.
3. Open the **Segment Editor** module, create a new segmentation, and define
   all organ/ROI segments on any representative frame.  The same segments
   will be applied to every frame automatically.

### 2 — Run Dynamic PET Quantifier

| Field | What to select |
|-------|----------------|
| **Sequence Browser** | The browser node created in step 2 above |
| **PET Proxy Volume** | The scalar volume proxy linked to that browser |
| **Segmentation** | The segmentation node with your organ segments |
| **Frame Duration** | (Optional) Uniform frame duration in seconds — leave 0 to use frame index |
| **Statistics** | Check/uncheck the metrics you need |
| **Output Format** | CSV or Excel |
| **Output File** | Path and filename for the result file |

Click **▶ Run Quantification**.  A progress bar tracks each frame.

### 3 — Output file structure

Each row represents one frame.  Columns are:

```
Frame | Time_s | <Seg1>_Mean | <Seg1>_StdDev | ... | <SegN>_Max | ...
```

- `Frame` — zero-based frame index
- `Time_s` — mid-frame time in seconds (or frame index when duration = 0)
- One column per (segment × statistic) combination

Example (3 segments, mean + stdev only):

| Frame | Time_s | Liver_Mean | Liver_StdDev | Kidney_Mean | Kidney_StdDev | Heart_Mean | Heart_StdDev |
|-------|--------|-----------|-------------|------------|--------------|-----------|-------------|
| 0     | 15.0   | 1.23      | 0.45        | 0.98       | 0.31         | 3.40      | 1.12        |
| 1     | 45.0   | 1.56      | 0.51        | 1.10       | 0.28         | 4.20      | 1.30        |
| …     | …      | …         | …           | …          | …            | …         | …           |

---

## Using the Logic from a Python Script / Jupyter

You can call the logic headlessly (e.g. from the Slicer Python console):

```python
import importlib, sys
sys.path.insert(0, "/path/to/folder/containing/DynamicPETQuantifier.py")
import DynamicPETQuantifier
importlib.reload(DynamicPETQuantifier)

logic = DynamicPETQuantifier.DynamicPETQuantifierLogic()
logic.run(
    browserNode      = slicer.util.getNode("SequenceBrowser"),
    proxyVolume      = slicer.util.getNode("PET_proxy"),
    segmentationNode = slicer.util.getNode("Segmentation"),
    outputPath       = "/tmp/pet_results.xlsx",
    frameDurationSec = 30.0,          # 30 s frames
    statsFlags       = {
        "voxel_count": False,
        "volume_mm3":  False,
        "minimum":     False,
        "maximum":     True,
        "mean":        True,
        "stdev":       True,
        "median":      False,
    },
    useExcel = True,
)
```

---

## Requirements

| Requirement | Notes |
|-------------|-------|
| 3D Slicer | 5.0 or later |
| Sequences extension | Included in Slicer 5.x by default |
| SegmentStatistics module | Included in Slicer 5.x by default |
| openpyxl | Bundled with Slicer 5.x — required only for Excel output |

---

## Troubleshooting

**"Could not find a Sequence node linked to the selected proxy volume"**  
→ Make sure you selected the correct Sequence Browser *and* the correct proxy
volume.  Open the Sequences module and verify that the volume appears in the
browser's proxy list.

**Excel output falls back to CSV**  
→ `openpyxl` is not available in your Slicer build.  The logic automatically
writes a `.csv` file instead.  Install openpyxl via the Slicer Python
interactor: `pip_install("openpyxl")`.

**Statistics columns are empty**  
→ Ensure the segmentation has a binary labelmap representation.  In Segment
Editor, right-click a segment → *Export to binary labelmap*.

---

## License

MIT License — free to use, modify, and distribute.

# SlicerWorkflow — Custom Workflow Orchestrator for 3D Slicer

A scripted module that lets you build, save, and run **multi-step image-analysis
pipelines** with pre-configured parameters on single volumes or entire folders.

---

## Features

| Feature | Detail |
|---------|--------|
| **Pipeline builder** | Add / remove / reorder steps by drag-and-drop |
| **Per-step parameters** | Typed widgets (spinbox, enum, path, node selector, code editor) |
| **Save / load workflows** | JSON files stored in your Slicer user-data folder |
| **Import / Export** | Share `.json` workflow files with collaborators |
| **Single-volume mode** | Process any node currently loaded in the scene |
| **Batch-folder mode** | Walk a directory and process every matching file |
| **Built-in steps** | Registration, N4 Bias, Gaussian/Median filtering, Threshold & Otsu segmentation, Volume & Intensity statistics, Save Volume, Custom Python |
| **Extensible** | Add your own steps by editing the `STEP_CATALOGUE` dict |

---

## Installation

1. Copy the `SlicerWorkflow/` folder anywhere on your machine.
2. Open 3D Slicer → **Edit → Application Settings → Modules**
3. Under **Additional module paths**, click **Add** and select the `SlicerWorkflow/` folder.
4. Click **OK** and **restart Slicer** when prompted.
5. The module now appears under **Modules → Utilities → Workflow** (or search "Workflow").

---

## Quick Start

### Build a simple denoise → segment → quantify pipeline

```
1. Open the Workflow module
2. Click ➕ (New Workflow), name it "MyPipeline"
3. Add steps:
     • Gaussian Smoothing  (sigma = 1.0)
     • Otsu Segmentation
     • Volume Statistics   (saveCSV = ✓, outputDir = /my/results)
4. Choose a volume node in the "Input" panel
5. Click ▶ Run Workflow
```

### Batch-process a folder

```
1. Select "Batch Folder" radio button
2. Browse to your folder (e.g. /data/study_group_A/)
3. Set file filter: *.nrrd *.nii.gz
4. Set output dir: /data/results/
5. Click ▶ Run Workflow
```

---

## Built-in Steps Reference

### Registration
Uses **BRAINSFit** CLI module.
| Param | Type | Default |
|-------|------|---------|
| transformType | Rigid / Affine / BSpline | Rigid |
| samplingPercentage | 0.01–1.0 | 0.25 |
| numberOfIterations | int | 200 |
| fixedVolumeNode | node selector | — |

### N4 Bias Correction
SimpleITK N4ITK MRI bias-field correction.
| Param | Default |
|-------|---------|
| numberOfIterations | 50 |
| convergenceThreshold | 0.001 |
| shrinkFactor | 4 |

### Gaussian Smoothing
`sitk.SmoothingRecursiveGaussian`  — `sigma` in mm.

### Median Filtering
`sitk.Median` — isotropic `radius` in voxels.

### Threshold Segmentation
Creates a segmentation node from an intensity range `[lowerThreshold, upperThreshold]`.

### Otsu Segmentation
Automatic two-class Otsu threshold via SimpleITK.

### Volume Statistics / Intensity Statistics
Uses **SegmentStatistics** logic. Exports CSV to `outputDir`.

### Save Volume
Saves current volume to `outputDir` in `nrrd`, `nii.gz`, or `mha` format.

### Custom Python Script
Arbitrary code with `volumeNode`, `segmentationNode`, `slicer`, and `logFn` injected.

---

## Adding Your Own Steps

Edit `STEP_CATALOGUE` in `WorkflowModule.py`:

```python
"My Custom Filter": {
    "description": "Applies my special algorithm",
    "params": {
        "strength": {"type": "float", "min": 0, "max": 10, "default": 2.5},
        "mode":     {"type": "enum", "values": ["fast","accurate"], "default": "fast"},
    },
},
```

Then add the execution logic in `_executeStep()`:

```python
elif stepName == "My Custom Filter":
    strength = float(params.get("strength", 2.5))
    # ... your SimpleITK / vtk / numpy code ...
    return outVolumeNode, segNode
```

---

## Workflow JSON Format

Workflows are stored in `<Slicer user dir>/SlicerWorkflows/<name>.json`:

```json
{
  "steps": [
    {
      "step": "Gaussian Smoothing",
      "params": { "sigma": 1.5 }
    },
    {
      "step": "Threshold Segmentation",
      "params": {
        "lowerThreshold": 200,
        "upperThreshold": 3000,
        "segmentName": "Bone"
      }
    },
    {
      "step": "Volume Statistics",
      "params": { "saveCSV": true, "outputDir": "/data/results" }
    }
  ]
}
```

---

## Extending with Claude Code

The recommended development workflow:

1. **Install Claude Code** in VS Code and open the `SlicerWorkflow/` folder.
2. Ask Claude to add new steps, fix bugs, or create batch post-processing logic.
3. In Slicer, enable **Developer Mode** (`Edit → Application Settings → Developer`)
   to get the **Reload** button — no restart needed during iteration.
4. Paste error messages back to Claude Code for rapid debugging.

---

## Requirements

- 3D Slicer ≥ 5.x (tested on 5.6)
- Bundled Slicer Python packages: `SimpleITK`, `vtk`, `qt`, `ctk`, `numpy`
- Optional: `BRAINSFit` extension (for Registration step)
- Optional: `SegmentStatistics` module (ships with Slicer core)

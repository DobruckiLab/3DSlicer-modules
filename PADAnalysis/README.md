# PAD Analysis — 3D Slicer Module

A scripted 3D Slicer module for **preclinical hindlimb PET/CT** studies (e.g.
peripheral arterial disease / hindlimb ischemia models). It splits a
manually-drawn whole-hindlimb CT segmentation into three anatomical zones —
**Proximal**, **Distal**, and **Foot** — using cut-off planes you position
interactively on the CT, then quantifies PET tracer uptake in each zone.

---

## Concept

```
                 upper thigh / hip cutoff plane
                              │
        (excluded — body)    │   Proximal zone
        ─────────────────────┼──────────────────
                              │
                         knee plane
                              │   Distal zone
                     ─────────┼──────────────────
                              │
                        ankle plane
                              │   Foot zone
                     ─────────┴──────────────────
```

Each leg (left and/or right) gets its own set of three planes. Everything
proximal to the hip/upper-thigh plane is **excluded** from all zones, so
trunk/body signal (e.g. bladder, muscle mass near the hip) doesn't
contaminate the Proximal hindlimb measurement.

---

## Installation (no build required)

1. Copy the `PADAnalysis/` folder anywhere on your computer.
2. Open 3D Slicer → **Edit → Application Settings → Modules**.
3. Under *Additional module paths*, click **Add** and select the
   `PADAnalysis/` folder.
4. Click **OK** and restart Slicer when prompted.
5. The module appears under **Modules → Quantification → PAD Analysis**.

---

## Workflow

### 1 — Prepare the whole-limb segmentation

Before opening PAD Analysis, use **Segment Editor** on your CT volume to
create a segmentation with one or both of:

- a segment covering the **whole left hindlimb** (hip to toes)
- a segment covering the **whole right hindlimb**

Thresholding (bone/soft tissue) followed by paint/scissors/islands cleanup
works well. Name the segments anything you like — you'll pick them from a
dropdown in the next step. This module does not modify these original
segments.

### 2 — Open PAD Analysis and set inputs

| Field | What to select |
|-------|----------------|
| **CT Volume** | The CT the whole-limb segmentation was drawn on |
| **PET Volume** | The PET volume to quantify uptake from |
| **Whole-Limb Segmentation** | The segmentation node from step 1 |
| **Left/Right Leg Segment** | Which segment in it is the left / right leg (choose "None" to skip a side) |

### 3 — Place the cut-off planes

For each leg you're processing, and for each of the three rows
(**Proximal Boundary**, **Knee Joint**, **Ankle / Foot Joint**):

1. Click **Create**. A plane appears on the CT with interactive handles, and
   the slice views jump to it.
2. Drag the center handle to reposition it, and the rotation handle to
   orient it, directly on the CT image, at:
   - **Proximal Boundary** — a consistent point in the upper thigh/hip
     region, above which tissue should be excluded entirely
   - **Knee Joint** — the knee
   - **Ankle / Foot Joint** — the ankle
3. Check that the plane's **normal arrow points toward the foot**. If it
   points toward the body instead, click that row's **Flip Normal** button
   to reverse it immediately.

Use **Jump To** to re-center the slice views on an existing plane, and
**Delete** to remove and start over. Planes are fully interactive — drag the
center handle to translate and the rotation handles to tilt the plane to
match the leg's actual orientation in that scan, which can vary from image
to image.

### 4 — Run

Choose/​create an **Output Segmentation** (this is where the six zone
segments — `Left_Proximal`, `Left_Distal`, `Left_Foot`, `Right_Proximal`,
`Right_Distal`, `Right_Foot` — will be written; existing zone segments with
the same name are replaced on re-run), pick which **PET Statistics** to
compute, choose an **Output CSV** path, and click **▶ Run PAD Analysis**.

### 5 — Output

- The output segmentation gains up to six new segments (three per
  processed leg), colored by zone.
- A CSV is written with one row per zone segment:

| Segment | Side | Zone | VoxelCount | Volume_mm3 | Min | Max | Mean | StdDev |
|---------|------|------|-----------|------------|-----|-----|------|--------|
| Left_Proximal | Left | Proximal | 18420 | 2211.4 | 0.01 | 4.32 | 1.85 | 0.61 |
| Left_Distal   | Left | Distal   | 9105  | 1092.6 | 0.02 | 5.10 | 2.40 | 0.88 |
| Left_Foot     | Left | Foot     | 3320  | 398.4  | 0.05 | 3.77 | 1.62 | 0.55 |
| …             | …    | …        | …     | …      | …    | …    | …    | …      |

Statistics are computed directly from the PET volume's own units (e.g.
Bq/cc if calibrated via SiemensImporter) — no SUV normalization is applied.

---

## Troubleshooting

**"No zone voxels were produced for any leg"**
→ The planes are likely mis-ordered or mis-oriented. Check that, for each
leg, Proximal → Knee → Ankle progress consistently toward the foot, and
that each plane's normal points toward the foot (use **Flip Normal** if not).

**A zone segment has 0 voxels**
→ That plane pair doesn't intersect the whole-limb segment as expected —
re-check its position, or that the Proximal/Knee/Ankle planes don't cross
each other.

**"The CT volume has an active transform"**
→ Harden the transform first: **Data module** → right-click the CT volume
→ **Harden Transform**.

**Statistics columns are empty**
→ Ensure the whole-limb segmentation has a binary labelmap representation
that overlaps the CT geometry (Segment Editor operates on this
automatically).

---

## Using the Logic from a Python Script / Jupyter

```python
import importlib, sys
sys.path.insert(0, "/path/to/folder/containing/PADAnalysis.py")
import PADAnalysis
importlib.reload(PADAnalysis)

logic = PADAnalysis.PADAnalysisLogic()

ctVolume  = slicer.util.getNode("CT")
petVolume = slicer.util.getNode("PET")
wholeLimb = slicer.util.getNode("WholeLimbSegmentation")
outputSeg = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "PADAnalysis_Segmentation")

# Planes can also be created interactively in the GUI first, then reused here:
proximalPlane = slicer.util.getNode("PAD_Left_Proximal")
kneePlane     = slicer.util.getNode("PAD_Left_Knee")
anklePlane    = slicer.util.getNode("PAD_Left_Ankle")

zoneCounts = logic.run(
    ctVolumeNode=ctVolume,
    petVolumeNode=petVolume,
    wholeLimbSegmentation=wholeLimb,
    legSegmentIDs={"L": "Segment_1"},
    planeNodes={"L": {"proximal": proximalPlane, "knee": kneePlane, "ankle": anklePlane}},
    outputSegmentationNode=outputSeg,
    outputCsvPath="/tmp/pad_results.csv",
)
print(zoneCounts)
```

---

## Requirements

| Requirement | Notes |
|-------------|-------|
| 3D Slicer | 5.2 or later (needs interactive Markups plane handles) |
| Markups extension | Included in Slicer core |
| SegmentStatistics module | Included in Slicer core |
| numpy | Bundled with Slicer |

---

## License

MIT License — free to use, modify, and distribute.

# PAD Analysis 2 — 3D Slicer Module

Semi-automatic hindlimb PET/CT analysis for preclinical **peripheral arterial
disease / hindlimb ischemia** studies. From the CT, the module builds a body
mask automatically and, from **3 clicks per leg** (hip, knee, ankle), splits
each hindlimb into **Proximal**, **Distal** and **Foot** zones. It then
quantifies PET uptake in each zone as **%ID/g**, **SUV** and
**ischemic / contralateral ratios**.

Compared with [PAD Analysis](../PADAnalysis/), you no longer have to draw the
whole limb by hand or position cut planes. The planes are derived from the
landmarks.

Assumes a **supine** animal with legs extended, a CT in **HU**, and PET
already aligned with the CT.

---

## Workflow

| Step | What happens | Output |
|------|--------------|--------|
| 1. Body mask | CT > threshold (default −300 HU) → opening (0.5 mm) → largest connected component → fill holes, on an isotropic working grid (default 0.25 mm) | `BodyMask` segment (hidden later); edit it in Segment Editor if the bed is still attached |
| 2. Landmarks | Click **Place**, then hip joint → knee joint → ankle joint on the CT, for each leg | `PAD2_Left_Landmarks`, `PAD2_Right_Landmarks` |
| 3. Zones | Body mask split per leg (see below) | `Left/Right_Proximal`, `_Distal`, `_Foot` segments |
| 4. Quantify | PET statistics per zone → %ID/g, SUV, ratios | CSV + `_ratios.csv`, tables `PAD2_Statistics`, `PAD2_Ratios` |

The zone segments can be corrected in Segment Editor before Step 4.

### How the zones are defined

H = hip, K = knee and A = ankle landmarks; *u* = body long axis (principal
axis of the body mask, pointing toward the feet).

- **Foot**: body voxels beyond the ankle plane (through A, normal along the
  tibia axis K→A), within *foot reach* of A (default = knee–ankle distance).
- **Distal**: body voxels between the knee plane and the ankle plane, within
  the *capture radius* of the K–A axis (default 35 % of the knee–ankle
  distance). The knee plane goes through K, perpendicular to the bisector of
  the femur and tibia directions.
- **Proximal**: body voxels on the knee side of a transverse plane placed
  *hip margin* mm (default 2 mm) above H, and above the knee plane. A sagittal
  **medial boundary** sits at a set fraction of the midline→hip distance
  (default 50 %). This gives the box-shaped thigh region of the reference
  figure and keeps bladder, genital and tail signal out.
- Each zone keeps only the connected piece nearest to its landmarks, so the
  tail or the other leg is not picked up. Zones never overlap.

### Quantification

```
net dose   = (syringe − residual) [× exp(−ln2 · Δt / T½) if decay-correcting to scan start]
%ID/g      = C [Bq/mL] / net dose [Bq] × 100            (1 g/mL tissue density)
SUV        = C [Bq/mL] × body weight [g] / net dose [Bq]
ratio      = ischemic / contralateral   (per zone, for Mean and Max)
```

- PET values are sampled with nearest-neighbour lookup of the zone masks at
  the PET voxel centres.
- **Exclude skin margin** shrinks each zone before measuring, to reduce
  partial-volume spill-out at the skin. The segments themselves are not
  changed.
- Supported PET units: Bq/mL, kBq/mL, MBq/mL, nCi/mL, µCi/mL.
- Supported dose units: MBq, kBq, mCi, µCi.
- Isotope half-lives are preset for F-18, C-11, N-13, O-15, Ga-68, Cu-64,
  Zr-89 and I-124. A custom half-life can be entered.
- **Append** mode adds rows to an existing CSV, with a `SubjectID` column, so
  one file can collect a whole cohort.

**CSV columns:** `SubjectID, Segment, Side, Zone, Condition, Volume_mm3,
PET_Voxels, Mean, StdDev, Median, Min, Max, PET_Units, %ID/g_Mean,
%ID/g_StdDev, %ID/g_Max, SUV_Mean, SUV_StdDev, SUV_Max, NetDose_MBq,
BodyWeight_g`

**Ratios CSV:** `SubjectID, Zone, IschemicSide, ControlSide, Ratio_Mean,
Ratio_Max`

---

## Installation

1. In 3D Slicer, go to **Edit → Application Settings → Modules → Additional
   module paths** and add the `PADAnalysis2/` folder.
2. Restart Slicer.

The module appears under **Modules → Quantification → PAD Analysis 2**.

---

## User Guide

This guide covers one animal from start to finish. Expect about 2–3 minutes
per animal once you are familiar with it.

### Before you start: what to prepare

| Item | Where it comes from |
|------|---------------------|
| CT volume in HU | scanner export (DICOM / NIfTI / NRRD) |
| PET volume aligned with the CT | same session; check the overlay (below) |
| PET image units | scanner / reconstruction settings (usually Bq/mL) |
| Syringe activity before and after injection, and their unit | dose-calibrator log |
| Body weight in g | animal record |
| Time from injection to PET scan start in min, and the isotope | study log |
| Which leg is ischemic (ligated) | study protocol |

### 0. Load and check the data

1. Load the CT and PET (**File → Add Data** or the **DICOM** module).
2. Show the CT as the *background* and the PET as the *foreground* in the
   slice views (layer dropdowns at the top of each view; set the opacity
   slider to ~50 %). Check that the PET signal sits on the anatomy.
3. If either volume is under a transform (for example from a manual
   registration), harden it: **Data** module → right-click the volume →
   **Harden transform**. The module refuses volumes with an active transform.
4. For landmark placement, a bone window on the CT helps: **Volumes** module →
   *Display* → *Window/Level* preset **CT-Bone** (or W ≈ 2000, L ≈ 500).

### 1. Inputs

Open **Modules → Quantification → PAD Analysis 2** and fill the **Inputs**
section:

- **CT Volume**: the CT.
- **PET Volume**: the PET.
- **Segmentation**: choose *Create new Segmentation* (it is named
  `PAD2_Segmentation`; rename it, for example `M01_PAD2`, if you like).
- **Subject ID**: an identifier such as `M01_Day7`. It is written to every
  CSV row and is what lets you pool many animals in one file.

### 2. Step 1 — Body Mask

1. Leave the defaults to start with:
   - Tissue threshold −300 HU
   - Opening radius 0.5 mm
   - Working resolution 0.25 mm
2. Click **Create Body Mask**. A red `BodyMask` segment appears.
3. **Check it** by scrolling through the axial and coronal views:
   - The whole animal, including feet and tail, should be covered.
   - The bed, nose cone and holder should not be.
4. If something is wrong:

   | Problem | Fix |
   |---------|-----|
   | Bed / holder included | Increase **Opening radius** to 1.0 mm and rerun, or open **Segment Editor**, select `BodyMask`, and remove the bed with **Scissors** (or **Islands → Keep largest island** after cutting). |
   | Thin parts (toes, tail) missing | Lower **Tissue threshold** (for example −400 HU) or set **Opening radius** to 0. |
   | Very slow on a large CT | Increase **Working resolution** to 0.3–0.5 mm. This is still much finer than the PET. |

   Edits made in Segment Editor are used by the next steps; you don't need to
   rerun Step 1 after editing.

### 3. Step 2 — Landmarks (3 clicks per leg)

Place each point in the **centre of the joint**. The coronal view is usually
the easiest, but check in a second view.

| Landmark | Where to click |
|----------|----------------|
| **Hip** | Centre of the femoral head, inside the acetabulum |
| **Knee** | Joint space between the femoral condyles and the tibial plateau (just behind the patella) |
| **Ankle** | Joint between the distal tibia/fibula and the talus |

1. Click **Place** next to *Left leg*. The mouse switches to placement mode.
2. Click **hip → knee → ankle**, in that order. The status line shows
   *1/3 … 2/3 …* and turns green with **✔ hip, knee, ankle placed**.
   Placement mode ends automatically after the third point.
3. Repeat for the *Right leg*.
4. **Adjust**: drag any point in a slice view. **Jump** centres the views on
   that leg's points.
   **Reset** clears a leg's three points so you can place them again.

> "Left" and "Right" mean the **animal's** left and right. In a supine animal
> shown in standard radiological orientation, the animal's left leg appears
> on the **right** side of the screen. The point labels (`L-Hip`, `R-Knee`, …)
> show which leg is which.

### 4. Step 3 — Generate Zones

1. Click **Generate Zones**. Six coloured segments appear:
   - Left leg (warm colours): `Left_Proximal`, `Left_Distal`, `Left_Foot`
   - Right leg (cool colours): `Right_Proximal`, `Right_Distal`, `Right_Foot`
   - The `BodyMask` is hidden.
2. Check the result in the **3D view**. Click **Show 3D** in the
   Segmentations module or in the segmentation's eye menu.
   - The thigh boxes should sit at the hips.
   - Each lower leg and foot should be coloured separately.
   - The tail and abdomen should not be coloured.
3. Fine-tune if needed, then click **Generate Zones** again (it replaces the
   previous zones):

   | Parameter | Default | Change it when … |
   |-----------|---------|------------------|
   | **Proximal top above hip** | 2 mm | Proximal box should start higher (+) or lower (−) |
   | **Medial boundary** | 50 % | Bladder / genital / tail-base signal leaks into Proximal → increase (60–70 %). Too little thigh muscle included → decrease |
   | **Limb capture radius** | Auto (35 % of knee–ankle length) | Distal misses part of the calf (increase) or picks up the tail (decrease) |
   | **Foot reach** | Auto (= knee–ankle length) | Toes cut off (increase) or the tail tip is included (decrease) |

4. Small local errors can be fixed directly in **Segment Editor** (Paint /
   Erase / Scissors on the zone segments). Quantification uses the segments
   exactly as they are.

### 5. Step 4 — PET Quantification

1. **PET image units**: units of the PET voxel values (check, for example,
   with the mouse-over value at the bottom of the slice view).
2. **Injected activity**:
   - **Syringe** is the activity before injection.
   - **Residual** is what was left in the syringe afterwards.
   - Select the unit (MBq, kBq, mCi or µCi).
   - Net dose = syringe − residual.
3. **Body weight** in grams.
4. **Decay**:
   - Keep **Decay-correct dose to scan start** checked if the PET is
     reconstructed with decay correction to scan start. This is the usual
     setting on most preclinical scanners.
   - Select the isotope (or *Custom* and type the half-life).
   - Enter **Injection → scan** time in minutes.
   - Uncheck it only if your PET is decay-corrected to injection time.
5. **Ischemic leg**: the ligated side. Ratios are ischemic ÷ contralateral.
6. **Exclude skin margin** (optional): 0.3–0.5 mm reduces spill-out at the
   skin edge for small limbs. Leave at 0 to measure the full zones.
7. **Output CSV**: choose a file. Keep **Append** checked to collect all
   animals of a study in the same file.
8. Click **Quantify PET**. The layout switches to show the `PAD2_Statistics`
   table, and two files are written:
   - `yourfile.csv` has one row per zone: volume, raw PET statistics, %ID/g,
     SUV, and the dose and weight used.
   - `yourfile_ratios.csv` has one row per zone pair: `Ratio_Mean` and
     `Ratio_Max` (ischemic / control).

   The ratios table is also available as `PAD2_Ratios` (select it in the
   table view).

### 6. Save and move to the next animal

1. **File → Save Data** saves the scene (`.mrb` is easiest). This keeps the
   CT, PET, landmarks and zones, so the analysis can be reviewed or redone
   later.
2. **File → Close Scene**, load the next animal, and repeat from step 0.
   Use a new **Subject ID** and the same **Output CSV**.
3. For longitudinal studies (for example Day 0 / 7 / 14), put the time point
   in the Subject ID (`M01_D7`) so it can be split out in Excel / R / Prism.

### Quick checklist

- [ ] PET and CT overlap correctly; no active transforms
- [ ] BodyMask covers the animal, not the bed
- [ ] 3 landmarks per leg, at joint centres, on the correct side
- [ ] Zones look right in 3D; tail and bladder excluded
- [ ] Units, dose, residual, weight, isotope and timing entered
- [ ] Correct ischemic leg selected
- [ ] Unique Subject ID; scene saved

---

## Scripted use

```python
import PADAnalysis2
logic = PADAnalysis2.PADAnalysis2Logic()
logic.createBodyMask(ctNode, segNode)
landmarks = {s: logic.landmarkPositions(slicer.util.getNode(logic.landmarkNodeName(s)))
             for s in ("L", "R")}
logic.generateZones(segNode, landmarks)
logic.quantify(segNode, petNode, petUnits="Bq/mL", syringeActivity=10.2,
               residualActivity=0.3, doseUnits="MBq", bodyWeightG=25.1,
               decayCorrect=True, halfLifeMin=109.77, injectionToScanMin=60,
               ischemicSide="L", subjectID="M01", outputCsvPath="/path/pad2.csv")
```

## Troubleshooting

- **Bed included in the body mask**: increase *Opening radius* (for example
  to 1 mm), or erase the bed from `BodyMask` in Segment Editor, then rerun
  Step 3.
- **Proximal zone picks up bladder signal**: increase *Medial boundary*
  (for example to 70 %).
- **Distal zone too thin or too wide**: set *Limb capture radius* manually.
- **Volume under a transform**: harden the transform first.

## Self-test

Open the module with **Developer mode** enabled and click **Reload and
Test**. The test builds a synthetic animal (torso, legs, tail and bed) and
checks the body mask, tail exclusion, zones, %ID/g, SUV, ratios and decay
correction.

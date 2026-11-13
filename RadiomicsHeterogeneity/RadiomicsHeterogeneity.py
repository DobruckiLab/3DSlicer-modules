"""Radiomics Heterogeneity module for 3D Slicer.

Computes volumetric (3D) radiomics heterogeneity features from a CT or PET
volume inside a selected segment of a segmentation node.

Feature definitions follow IBSI / PyRadiomics conventions:

* **First-order** statistics of the raw intensities inside the ROI.
* **Shape** — volume, surface area (marching cubes), sphericity, etc.
* **GLCM** — symmetric co-occurrence matrix, 13 3D directions, distance 1.
* **GLRLM** — run-length matrix, 13 3D directions.
* **GLSZM** — size-zone matrix, 26-connected 3D zones.
* **NGTDM** — neighbourhood grey-tone difference matrix, 26-neighbourhood.

Texture matrices are computed after fixed-bin-number (FBN) discretisation of
the ROI intensities into ``N`` grey levels (1..N).  For GLCM and GLRLM the
per-direction matrices are merged (summed) before feature calculation
(IBSI "3D merged" aggregation).
"""

from __future__ import annotations

import csv
import logging
import os

import ctk
import numpy as np
import qt
import slicer
import vtk
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
    ScriptedLoadableModuleWidget,
)

try:
    from scipy import ndimage as _ndimage
except ImportError:  # pragma: no cover - scipy ships with Slicer
    _ndimage = None


EPS = np.finfo(np.float64).eps

# 13 unique 3D directions (the other 13 are their negatives), in array
# index order (k, j, i) == (z, y, x).
DIRECTIONS_3D = [
    (0, 0, 1), (0, 1, 0), (1, 0, 0),
    (0, 1, 1), (0, 1, -1),
    (1, 0, 1), (1, 0, -1),
    (1, 1, 0), (1, -1, 0),
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1),
]


# ============================================================================
# Array helpers
# ============================================================================


def _shifted_pair_slices(shape, offset):
    """Slices (a, b) such that arr[b] is the neighbour of arr[a] at +offset."""
    a, b = [], []
    for n, d in zip(shape, offset):
        if d >= 0:
            a.append(slice(0, n - d))
            b.append(slice(d, n))
        else:
            a.append(slice(-d, n))
            b.append(slice(0, n + d))
    return tuple(a), tuple(b)


def _shift(arr, offset, fill=0):
    """out[x] = arr[x - offset] (zero/fill outside)."""
    out = np.full_like(arr, fill)
    src, dst = _shifted_pair_slices(arr.shape, offset)
    out[dst] = arr[src]
    return out


def crop_to_mask(intensity, mask, pad=1):
    """Crop both arrays to the mask bounding box plus ``pad`` voxels."""
    idx = np.nonzero(mask)
    lo = [max(int(i.min()) - pad, 0) for i in idx]
    hi = [min(int(i.max()) + pad + 1, n) for i, n in zip(idx, mask.shape)]
    sl = tuple(slice(a, b) for a, b in zip(lo, hi))
    return intensity[sl], mask[sl]


def discretize(intensity, mask, num_levels):
    """Fixed-bin-number discretisation of ROI intensities to 1..num_levels.

    Voxels outside the mask are set to 0.
    """
    roi = intensity[mask]
    vmin, vmax = float(roi.min()), float(roi.max())
    levels = np.zeros(intensity.shape, dtype=np.int32)
    if vmax == vmin:
        levels[mask] = 1
        return levels
    q = np.floor(num_levels * (roi - vmin) / (vmax - vmin)).astype(np.int32) + 1
    levels[mask] = np.clip(q, 1, num_levels)
    return levels


# ============================================================================
# First-order features
# ============================================================================


def compute_first_order(values: np.ndarray, num_bins: int = 64) -> dict:
    """First-order statistics of the raw ROI intensities."""
    v = values.astype(np.float64)
    n = v.size
    mean = v.mean()
    std = v.std()
    p5, p10, p25, median, p75, p90, p95 = np.percentile(v, [5, 10, 25, 50, 75, 90, 95])
    if std > 0:
        z = (v - mean) / std
        skew = float(np.mean(z ** 3))
        kurt = float(np.mean(z ** 4)) - 3.0  # excess kurtosis (IBSI)
    else:
        skew = kurt = 0.0

    # Robust MAD: within the 10-90th percentile range (PyRadiomics)
    rv = v[(v >= p10) & (v <= p90)]
    rmad = float(np.mean(np.abs(rv - rv.mean()))) if rv.size else 0.0

    hist, _ = np.histogram(v, bins=num_bins)
    p = hist[hist > 0] / n
    entropy = float(-np.sum(p * np.log2(p)))
    uniformity = float(np.sum(p ** 2))

    return {
        "firstorder_VoxelCount": int(n),
        "firstorder_Mean": float(mean),
        "firstorder_StdDev": float(std),
        "firstorder_Variance": float(v.var()),
        "firstorder_Median": float(median),
        "firstorder_Minimum": float(v.min()),
        "firstorder_Maximum": float(v.max()),
        "firstorder_Range": float(v.max() - v.min()),
        "firstorder_P5": float(p5),
        "firstorder_P10": float(p10),
        "firstorder_P25": float(p25),
        "firstorder_P75": float(p75),
        "firstorder_P90": float(p90),
        "firstorder_P95": float(p95),
        "firstorder_InterquartileRange": float(p75 - p25),
        "firstorder_MeanAbsoluteDeviation": float(np.mean(np.abs(v - mean))),
        "firstorder_RobustMeanAbsoluteDeviation": rmad,
        "firstorder_MedianAbsoluteDeviation": float(np.median(np.abs(v - median))),
        "firstorder_RootMeanSquared": float(np.sqrt(np.mean(v ** 2))),
        "firstorder_Skewness": skew,
        "firstorder_Kurtosis": kurt,
        "firstorder_Entropy": entropy,
        "firstorder_Uniformity": uniformity,
        "firstorder_CoefficientOfVariation": float(std / mean) if mean != 0 else 0.0,
        "firstorder_QuartileCoefficientOfDispersion": (
            float((p75 - p25) / (p75 + p25)) if (p75 + p25) != 0 else 0.0
        ),
    }


# ============================================================================
# Shape features
# ============================================================================


def compute_shape(mask: np.ndarray, spacing_kji) -> dict:
    """Shape features. ``spacing_kji`` is voxel spacing in array order."""
    count = int(mask.sum())
    voxel_volume = float(np.prod(spacing_kji))
    volume = count * voxel_volume

    idx = np.nonzero(mask)
    extent = [(int(i.max()) - int(i.min()) + 1) * s for i, s in zip(idx, spacing_kji)]
    bb_z, bb_y, bb_x = extent

    # Surface area / mesh volume via marching cubes on the padded mask
    padded = np.pad(mask.astype(np.uint8), 1)
    image = vtk.vtkImageData()
    image.SetDimensions(padded.shape[2], padded.shape[1], padded.shape[0])
    image.SetSpacing(spacing_kji[2], spacing_kji[1], spacing_kji[0])
    from vtk.util import numpy_support
    image.GetPointData().SetScalars(
        numpy_support.numpy_to_vtk(padded.ravel(), deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    )
    mc = vtk.vtkMarchingCubes()
    mc.SetInputData(image)
    mc.SetValue(0, 0.5)
    mc.ComputeNormalsOff()
    tri = vtk.vtkTriangleFilter()
    tri.SetInputConnection(mc.GetOutputPort())
    props = vtk.vtkMassProperties()
    props.SetInputConnection(tri.GetOutputPort())
    props.Update()
    area = float(props.GetSurfaceArea())
    mesh_volume = float(props.GetVolume())

    sphericity = (np.pi ** (1 / 3)) * (6 * mesh_volume) ** (2 / 3) / area if area > 0 else 0.0

    return {
        "shape_VoxelCount": count,
        "shape_VoxelVolume_mm3": volume,
        "shape_MeshVolume_mm3": mesh_volume,
        "shape_SurfaceArea_mm2": area,
        "shape_SurfaceVolumeRatio": area / mesh_volume if mesh_volume > 0 else 0.0,
        "shape_Sphericity": float(sphericity),
        "shape_BoundingBox_X_mm": float(bb_x),
        "shape_BoundingBox_Y_mm": float(bb_y),
        "shape_BoundingBox_Z_mm": float(bb_z),
    }


# ============================================================================
# GLCM
# ============================================================================


def compute_glcm(levels: np.ndarray, num_levels: int) -> np.ndarray:
    """Symmetric GLCM merged over the 13 3D directions (distance 1).

    ``levels`` uses 0 for background and 1..N inside the ROI.
    """
    glcm = np.zeros((num_levels, num_levels), dtype=np.float64)
    for d in DIRECTIONS_3D:
        a, b = _shifted_pair_slices(levels.shape, d)
        la, lb = levels[a], levels[b]
        valid = (la > 0) & (lb > 0)
        i0 = la[valid] - 1
        i1 = lb[valid] - 1
        glcm += np.bincount(
            i0 * num_levels + i1, minlength=num_levels * num_levels
        ).reshape(num_levels, num_levels)
    return glcm + glcm.T  # symmetric


GLCM_KEYS = [
    "glcm_JointEnergy", "glcm_JointEntropy", "glcm_Contrast", "glcm_Correlation",
    "glcm_Idm", "glcm_Id", "glcm_Idmn", "glcm_Idn", "glcm_InverseVariance",
    "glcm_DifferenceEntropy", "glcm_SumEntropy", "glcm_ClusterShade",
    "glcm_ClusterProminence", "glcm_Imc1", "glcm_Imc2", "glcm_Autocorrelation",
]


def compute_glcm_features(glcm: np.ndarray) -> dict:
    total = glcm.sum()
    if total == 0:
        return {k: 0.0 for k in GLCM_KEYS}
    p = glcm / total
    ng = p.shape[0]
    g = np.arange(1, ng + 1, dtype=np.float64)
    i, j = np.meshgrid(g, g, indexing="ij")

    px = p.sum(axis=1)
    py = p.sum(axis=0)
    ux = np.sum(g * px)
    uy = np.sum(g * py)
    sx = np.sqrt(np.sum((g - ux) ** 2 * px))
    sy = np.sqrt(np.sum((g - uy) ** 2 * py))

    diff = np.abs(i - j).astype(int)
    p_diff = np.bincount(diff.ravel(), weights=p.ravel(), minlength=ng)  # k = 0..ng-1
    ssum = (i + j).astype(int)
    p_sum = np.bincount(ssum.ravel(), weights=p.ravel(), minlength=2 * ng + 1)[2:]  # k=2..2ng
    k_diff = np.arange(ng)

    def H(arr):
        a = arr[arr > 0]
        return float(-np.sum(a * np.log2(a)))

    hxy = H(p)
    hx, hy = H(px), H(py)
    pxpy = np.outer(px, py)
    nz = pxpy > 0
    hxy1 = float(-np.sum(p[nz] * np.log2(pxpy[nz] + EPS)))
    hxy2 = float(-np.sum(pxpy[nz] * np.log2(pxpy[nz] + EPS)))

    f = {}
    f["glcm_JointEnergy"] = float(np.sum(p ** 2))
    f["glcm_JointEntropy"] = hxy
    f["glcm_Contrast"] = float(np.sum(p * (i - j) ** 2))
    f["glcm_Correlation"] = (
        float((np.sum(p * i * j) - ux * uy) / (sx * sy)) if sx * sy > 0 else 1.0
    )
    f["glcm_Idm"] = float(np.sum(p / (1 + (i - j) ** 2)))
    f["glcm_Id"] = float(np.sum(p / (1 + np.abs(i - j))))
    f["glcm_Idmn"] = float(np.sum(p / (1 + (i - j) ** 2 / ng ** 2)))
    f["glcm_Idn"] = float(np.sum(p / (1 + np.abs(i - j) / ng)))
    f["glcm_InverseVariance"] = float(np.sum(p_diff[1:] / k_diff[1:] ** 2))
    f["glcm_DifferenceEntropy"] = H(p_diff)
    f["glcm_SumEntropy"] = H(p_sum)
    f["glcm_ClusterShade"] = float(np.sum(p * (i + j - ux - uy) ** 3))
    f["glcm_ClusterProminence"] = float(np.sum(p * (i + j - ux - uy) ** 4))
    denom = max(hx, hy)
    f["glcm_Imc1"] = float((hxy - hxy1) / denom) if denom > 0 else 0.0
    f["glcm_Imc2"] = float(np.sqrt(max(0.0, 1 - np.exp(-2 * (hxy2 - hxy)))))
    f["glcm_Autocorrelation"] = float(np.sum(p * i * j))
    return f


# ============================================================================
# Run / zone matrix helpers
# ============================================================================


def compute_glrlm(levels: np.ndarray, num_levels: int) -> np.ndarray:
    """GLRLM merged over 13 3D directions. Shape (num_levels, max_run)."""
    max_run = max(levels.shape)
    rlm = np.zeros((num_levels, max_run), dtype=np.float64)
    roi = levels > 0
    for d in DIRECTIONS_3D:
        neg = tuple(-x for x in d)
        # cont[x]: voxel x continues the run of its predecessor x - d
        cont = roi & (_shift(levels, d) == levels)
        # nxt[x]: the successor x + d continues the run from x
        nxt = _shift(cont, neg, fill=False)
        # Position within run, propagated along the direction.
        pos = roi.astype(np.int32)
        active = cont.copy()
        while active.any():
            new = np.where(cont, _shift(pos, d) + 1, pos)
            active = new != pos
            pos = new
        ends = roi & ~nxt
        g = levels[ends] - 1
        r = pos[ends] - 1
        np.add.at(rlm, (g, r), 1)
    # Trim unused run lengths
    used = np.nonzero(rlm.sum(axis=0))[0]
    return rlm[:, : (used.max() + 1 if used.size else 1)]


def compute_glszm(levels: np.ndarray, num_levels: int) -> np.ndarray:
    """GLSZM with 26-connectivity. Shape (num_levels, max_zone_size)."""
    if _ndimage is None:
        raise RuntimeError("scipy is required for GLSZM computation")
    structure = np.ones((3, 3, 3), dtype=bool)
    rows = {}
    max_size = 1
    for g in range(1, num_levels + 1):
        m = levels == g
        if not m.any():
            continue
        lab, n = _ndimage.label(m, structure=structure)
        sizes = np.bincount(lab.ravel())[1:]
        rows[g] = sizes
        max_size = max(max_size, int(sizes.max()))
    szm = np.zeros((num_levels, max_size), dtype=np.float64)
    for g, sizes in rows.items():
        szm[g - 1] += np.bincount(sizes - 1, minlength=max_size)
    return szm


def _size_matrix_features(m: np.ndarray, n_voxels: float, prefix: str, short: str) -> dict:
    """Shared GLRLM/GLSZM feature set (rows: grey level, cols: run/zone size)."""
    total = m.sum()
    names = [
        f"{short}ShortEmphasis", f"{short}LongEmphasis",
        "GrayLevelNonUniformity", "GrayLevelNonUniformityNormalized",
        f"{short}NonUniformity", f"{short}NonUniformityNormalized",
        f"{short}Percentage", "GrayLevelVariance", f"{short}Variance",
        f"{short}Entropy", "LowGrayLevelEmphasis", "HighGrayLevelEmphasis",
        f"{short}ShortLowGrayLevelEmphasis", f"{short}ShortHighGrayLevelEmphasis",
        f"{short}LongLowGrayLevelEmphasis", f"{short}LongHighGrayLevelEmphasis",
    ]
    if total == 0:
        return {f"{prefix}_{n}": 0.0 for n in names}
    p = m / total
    g = np.arange(1, m.shape[0] + 1, dtype=np.float64)[:, None]
    s = np.arange(1, m.shape[1] + 1, dtype=np.float64)[None, :]
    pg = m.sum(axis=1)
    ps = m.sum(axis=0)
    mu_g = np.sum(p * g)
    mu_s = np.sum(p * s)
    pn = p[p > 0]
    vals = [
        np.sum(p / s ** 2), np.sum(p * s ** 2),
        np.sum(pg ** 2) / total, np.sum(pg ** 2) / total ** 2,
        np.sum(ps ** 2) / total, np.sum(ps ** 2) / total ** 2,
        total / n_voxels, np.sum(p * (g - mu_g) ** 2), np.sum(p * (s - mu_s) ** 2),
        -np.sum(pn * np.log2(pn)), np.sum(p / g ** 2), np.sum(p * g ** 2),
        np.sum(p / (g ** 2 * s ** 2)), np.sum(p * g ** 2 / s ** 2),
        np.sum(p * s ** 2 / g ** 2), np.sum(p * g ** 2 * s ** 2),
    ]
    return {f"{prefix}_{n}": float(v) for n, v in zip(names, vals)}


def compute_glrlm_features(rlm, n_voxels):
    # Merged over 13 directions: each voxel belongs to 13 runs.
    return _size_matrix_features(rlm, n_voxels * len(DIRECTIONS_3D), "glrlm", "Run")


def compute_glszm_features(szm, n_voxels):
    return _size_matrix_features(szm, n_voxels, "glszm", "Zone")


# ============================================================================
# NGTDM
# ============================================================================

NGTDM_KEYS = [
    "ngtdm_Coarseness", "ngtdm_Contrast", "ngtdm_Busyness",
    "ngtdm_Complexity", "ngtdm_Strength",
]


def compute_ngtdm_features(levels: np.ndarray, num_levels: int) -> dict:
    """NGTDM features, 26-neighbourhood, neighbours restricted to the ROI."""
    if _ndimage is None:
        raise RuntimeError("scipy is required for NGTDM computation")
    roi = levels > 0
    kernel = np.ones((3, 3, 3), dtype=np.float64)
    kernel[1, 1, 1] = 0
    nsum = _ndimage.convolve(levels.astype(np.float64), kernel, mode="constant")
    ncnt = _ndimage.convolve(roi.astype(np.float64), kernel, mode="constant")
    valid = roi & (ncnt > 0)
    if not valid.any():
        return {k: 0.0 for k in NGTDM_KEYS}
    lv = levels[valid]
    a = np.abs(lv - nsum[valid] / ncnt[valid])
    n = np.bincount(lv - 1, minlength=num_levels).astype(np.float64)
    s = np.bincount(lv - 1, weights=a, minlength=num_levels)
    nvp = n.sum()
    present = n > 0
    p = n[present] / nvp
    s = s[present]
    g = np.arange(1, num_levels + 1, dtype=np.float64)[present]
    ngp = int(present.sum())

    pi, pj = np.meshgrid(p, p, indexing="ij")
    gi, gj = np.meshgrid(g, g, indexing="ij")
    si, sj = np.meshgrid(s, s, indexing="ij")
    ps_sum = float(np.sum(p * s))
    s_tot = float(s.sum())

    f = {}
    f["ngtdm_Coarseness"] = 1.0 / ps_sum if ps_sum > 0 else 1e6
    f["ngtdm_Contrast"] = (
        float(np.sum(pi * pj * (gi - gj) ** 2) / (ngp * (ngp - 1)) * s_tot / nvp)
        if ngp > 1 else 0.0
    )
    bden = float(np.sum(np.abs(gi * pi - gj * pj)))
    f["ngtdm_Busyness"] = ps_sum / bden if bden > 0 else 0.0
    f["ngtdm_Complexity"] = float(
        np.sum(np.abs(gi - gj) * (pi * si + pj * sj) / (pi + pj)) / nvp
    )
    f["ngtdm_Strength"] = (
        float(np.sum((pi + pj) * (gi - gj) ** 2) / s_tot) if s_tot > 0 else 0.0
    )
    return f


# ============================================================================
# Logic
# ============================================================================


class RadiomicsHeterogeneityLogic(ScriptedLoadableModuleLogic):
    """Feature extraction for one volume / segment pair."""

    STEPS = ["First-order", "Shape", "GLCM", "GLRLM", "GLSZM", "NGTDM"]

    @staticmethod
    def extractArrays(volumeNode, segmentationNode, segmentId):
        """Return (intensity, mask, spacing_kji) in the volume's geometry."""
        if not segmentId:
            raise ValueError("No segment selected.")
        if segmentationNode.GetSegmentation().GetSegment(segmentId) is None:
            raise ValueError(f"Segment '{segmentId}' not found in segmentation.")
        intensity = slicer.util.arrayFromVolume(volumeNode).astype(np.float64)
        mask = slicer.util.arrayFromSegmentBinaryLabelmap(
            segmentationNode, segmentId, volumeNode
        ).astype(bool)
        if mask.shape != intensity.shape:
            raise RuntimeError(
                f"Mask shape {mask.shape} does not match volume shape {intensity.shape}."
            )
        spacing_kji = tuple(reversed(volumeNode.GetSpacing()))
        return intensity, mask, spacing_kji

    def computeFeaturesFromArrays(self, intensity, mask, spacing_kji, numLevels=64,
                                  progressCallback=None):
        """Compute all features from numpy arrays (testable outside Slicer UI)."""
        def progress(step):
            if progressCallback:
                progressCallback(step, self.STEPS[step] if step < len(self.STEPS) else "Done")

        if not mask.any():
            raise ValueError("The selected segment contains no voxels in this volume.")

        intensity, mask = crop_to_mask(intensity, mask)
        n_voxels = int(mask.sum())
        features = {}

        progress(0)
        features.update(compute_first_order(intensity[mask]))
        progress(1)
        features.update(compute_shape(mask, spacing_kji))

        levels = discretize(intensity, mask, numLevels)
        progress(2)
        features.update(compute_glcm_features(compute_glcm(levels, numLevels)))
        progress(3)
        features.update(compute_glrlm_features(compute_glrlm(levels, numLevels), n_voxels))
        progress(4)
        features.update(compute_glszm_features(compute_glszm(levels, numLevels), n_voxels))
        progress(5)
        features.update(compute_ngtdm_features(levels, numLevels))
        progress(6)
        return features

    def computeFeatures(self, volumeNode, segmentationNode, segmentId, numLevels=64,
                        progressCallback=None):
        intensity, mask, spacing = self.extractArrays(volumeNode, segmentationNode, segmentId)
        return self.computeFeaturesFromArrays(
            intensity, mask, spacing, numLevels, progressCallback
        )


# ============================================================================
# Module
# ============================================================================


class RadiomicsHeterogeneity(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Radiomics Heterogeneity"
        self.parent.categories = ["Quantification"]
        self.parent.dependencies = ["Segmentations"]
        self.parent.contributors = [
            "Wawrzyniec Dobrucki <br> Experimental Molecular Imaging Laboratory "
            "<br> University of Illinois at Urbana-Champaign"
        ]
        self.parent.helpText = """
<p>Computes 3D radiomics heterogeneity features of a CT or PET volume inside a
selected segment.</p>
<ol>
  <li>Select the <b>Volume</b>.</li>
  <li>Select the <b>Segmentation</b> and the <b>Segment</b> to analyse.</li>
  <li>Choose the number of <b>Gray levels</b> used to discretise intensities
      for texture features (fixed bin number over the ROI range).</li>
  <li>Click <b>Compute Features</b>. Use <b>Export CSV</b> to save results.</li>
</ol>
<p>Feature classes: first-order, shape, GLCM and GLRLM (13 directions, merged),
GLSZM (26-connected zones) and NGTDM (26-neighbourhood). Definitions follow
IBSI / PyRadiomics; values may differ slightly from PyRadiomics because no
resampling is performed and directional matrices are merged.</p>
"""
        self.parent.acknowledgementText = (
            "Developed for volumetric radiomics analysis of preclinical PET/CT."
        )


# ============================================================================
# Widget
# ============================================================================


class RadiomicsHeterogeneityWidget(ScriptedLoadableModuleWidget):
    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        self.logic = None
        self._results: dict = {}
        self._resultsInfo: dict = {}

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = RadiomicsHeterogeneityLogic()

        # ── Parameters ─────────────────────────────────────────
        paramsBox = ctk.ctkCollapsibleButton()
        paramsBox.text = "Parameters"
        self.layout.addWidget(paramsBox)
        form = qt.QFormLayout(paramsBox)

        # MRML-aware selectors keep their selection across scene changes.
        self.volumeSelector = slicer.qMRMLNodeComboBox()
        self.volumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.volumeSelector.showChildNodeTypes = False  # exclude label maps
        self.volumeSelector.addEnabled = False
        self.volumeSelector.removeEnabled = False
        self.volumeSelector.noneEnabled = False
        self.volumeSelector.setMRMLScene(slicer.mrmlScene)
        self.volumeSelector.setToolTip("CT or PET volume to analyse.")
        form.addRow("Volume:", self.volumeSelector)

        self.segmentSelector = slicer.qMRMLSegmentSelectorWidget()
        self.segmentSelector.setMRMLScene(slicer.mrmlScene)
        self.segmentSelector.setToolTip("Segmentation and segment defining the ROI.")
        form.addRow("Segment:", self.segmentSelector)

        self.grayLevelsComboBox = qt.QComboBox()
        self.grayLevelsComboBox.addItems(["8", "16", "32", "64", "128"])
        self.grayLevelsComboBox.setCurrentText("64")
        self.grayLevelsComboBox.setToolTip(
            "Number of grey levels (fixed bin number over the ROI intensity range) "
            "used for GLCM / GLRLM / GLSZM / NGTDM."
        )
        form.addRow("Gray levels:", self.grayLevelsComboBox)

        self.computeButton = qt.QPushButton("Compute Features")
        form.addRow(self.computeButton)

        self.progressBar = qt.QProgressBar()
        self.progressBar.setRange(0, len(RadiomicsHeterogeneityLogic.STEPS))
        self.progressBar.hide()
        form.addRow(self.progressBar)

        self.statusLabel = qt.QLabel("")
        self.statusLabel.wordWrap = True
        form.addRow(self.statusLabel)

        # ── Results ────────────────────────────────────────────
        resultsBox = ctk.ctkCollapsibleButton()
        resultsBox.text = "Results"
        self.layout.addWidget(resultsBox)
        rform = qt.QVBoxLayout(resultsBox)

        self.resultsTable = qt.QTableWidget()
        self.resultsTable.setColumnCount(2)
        self.resultsTable.setHorizontalHeaderLabels(["Feature", "Value"])
        self.resultsTable.verticalHeader().setVisible(False)
        self.resultsTable.horizontalHeader().setStretchLastSection(True)
        self.resultsTable.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)
        self.resultsTable.setMinimumHeight(300)
        rform.addWidget(self.resultsTable)

        self.exportButton = qt.QPushButton("Export CSV")
        self.exportButton.enabled = False
        rform.addWidget(self.exportButton)

        self.layout.addStretch(1)

        # ── Connections ────────────────────────────────────────
        self.computeButton.connect("clicked()", self.onCompute)
        self.exportButton.connect("clicked()", self.onExport)
        self.volumeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateButtonStates)
        self.segmentSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateButtonStates)
        self.segmentSelector.connect("currentSegmentChanged(QString)", self.updateButtonStates)
        self.updateButtonStates()

    def cleanup(self):
        pass

    def updateButtonStates(self, *_):
        ok = (
            self.volumeSelector.currentNode() is not None
            and self.segmentSelector.currentNode() is not None
            and bool(self.segmentSelector.currentSegmentID())
        )
        self.computeButton.enabled = ok

    # ── Actions ───────────────────────────────────────────────

    def onCompute(self):
        # Read the user's selection once, up-front.
        volumeNode = self.volumeSelector.currentNode()
        segNode = self.segmentSelector.currentNode()
        segmentId = self.segmentSelector.currentSegmentID()
        numLevels = int(self.grayLevelsComboBox.currentText)
        if volumeNode is None or segNode is None or not segmentId:
            self.statusLabel.text = "Select a volume, segmentation and segment."
            return
        segmentName = segNode.GetSegmentation().GetSegment(segmentId).GetName()

        def progress(step, name):
            self.progressBar.value = step
            self.statusLabel.text = f"Computing {name}..." if step < 6 else "Done."
            slicer.app.processEvents()

        self.computeButton.enabled = False
        self.progressBar.value = 0
        self.progressBar.show()
        try:
            with slicer.util.tryWithErrorDisplay("Feature computation failed.", waitCursor=True):
                features = self.logic.computeFeatures(
                    volumeNode, segNode, segmentId, numLevels, progress
                )
                self._resultsInfo = {
                    "Volume": volumeNode.GetName(),
                    "Segmentation": segNode.GetName(),
                    "Segment": segmentName,
                    "GrayLevels": numLevels,
                }
                self._displayResults(features)
                self.statusLabel.text = (
                    f"{len(features)} features computed for segment '{segmentName}' "
                    f"in volume '{volumeNode.GetName()}'."
                )
        finally:
            self.progressBar.hide()
            self.updateButtonStates()

    def _displayResults(self, features):
        self._results = features
        rows = list(self._resultsInfo.items()) + list(features.items())
        self.resultsTable.setRowCount(len(rows))
        for r, (key, val) in enumerate(rows):
            text = f"{val:.6g}" if isinstance(val, float) else str(val)
            self.resultsTable.setItem(r, 0, qt.QTableWidgetItem(key))
            self.resultsTable.setItem(r, 1, qt.QTableWidgetItem(text))
        self.resultsTable.resizeColumnToContents(0)
        self.exportButton.enabled = True

    def onExport(self):
        if not self._results:
            return
        default_dir = os.path.expanduser("~")
        name = "radiomics_{}_{}.csv".format(
            self._resultsInfo.get("Volume", ""), self._resultsInfo.get("Segment", "")
        )
        name = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)
        path = qt.QFileDialog.getSaveFileName(
            slicer.util.mainWindow(), "Export Features",
            os.path.join(default_dir, name), "CSV Files (*.csv)",
        )
        if isinstance(path, (tuple, list)):  # PySide-style return
            path = path[0]
        if not path:
            return
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Feature", "Value"])
            for key, val in list(self._resultsInfo.items()) + list(self._results.items()):
                w.writerow([key, f"{val:.10g}" if isinstance(val, float) else val])
        self.statusLabel.text = f"Saved to {path}"


# ============================================================================
# Tests
# ============================================================================


class RadiomicsHeterogeneityTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_SelectedSegmentIsUsed()

    def test_SelectedSegmentIsUsed(self):
        """Two segments of different size must give different voxel counts."""
        self.delayDisplay("Creating synthetic volume")
        arr = np.random.default_rng(0).normal(100, 10, (40, 40, 40)).astype(np.float32)
        vol = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", "TestVol")
        slicer.util.updateVolumeFromArray(vol, arr)
        vol.SetSpacing(1.0, 1.0, 2.0)

        seg = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode", "TestSeg")
        seg.SetReferenceImageGeometryParameterFromVolumeNode(vol)
        ids = []
        for name, (lo, hi) in [("small", (10, 15)), ("large", (5, 30))]:
            m = np.zeros(arr.shape, dtype=np.uint8)
            m[lo:hi, lo:hi, lo:hi] = 1
            sid = seg.GetSegmentation().AddEmptySegment(name, name)
            slicer.util.updateSegmentBinaryLabelmapFromArray(m, seg, sid, vol)
            ids.append(sid)

        logic = RadiomicsHeterogeneityLogic()
        f_small = logic.computeFeatures(vol, seg, ids[0], 32)
        f_large = logic.computeFeatures(vol, seg, ids[1], 32)
        self.assertEqual(f_small["shape_VoxelCount"], 5 ** 3)
        self.assertEqual(f_large["shape_VoxelCount"], 25 ** 3)
        self.assertAlmostEqual(f_large["shape_VoxelVolume_mm3"], 25 ** 3 * 2.0)
        self.delayDisplay("Test passed")

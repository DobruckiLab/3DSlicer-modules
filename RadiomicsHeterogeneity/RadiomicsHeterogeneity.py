"""Radiomics Heterogeneity module for 3D Slicer.

Computes volumetric (3D) radiomics heterogeneity features from a CT or PET
volume using a user-provided segmentation (segment object, not labelmap).

Feature categories
------------------
* **First-order statistics** — mean, std, skewness, kurtosis, percentiles,
  entropy, CoV, robust CoV, etc.
* **3D GLCM** — contrast, correlation, energy, homogeneity, entropy,
  inverse-difference-moment variance, IMC1/IMC2 (averaged over all 13
  unique 3D neighbour directions).
* **3D GLRLM** — run length emphasis, short/long run emphasis, gray-level
  and run-length non-uniformities, etc.
* **3D GLZDM** — short/long zone emphasis, zone percentage, zone variance,
  zone entropy.
* **NGTDM** — coarseness, contrast, complexity, strength, granularity.
* **Morphological** — voxel count, volume (mm³), bounding-box dimensions,
  compactness, sphericity.
"""

from __future__ import annotations

import logging

import ctk
import numpy as np
import qt
import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleWidget,
)

from slicer import vtkMRMLScalarVolumeNode
from slicer import vtkMRMLSegmentationNode


# ============================================================================
# Core feature computation logic (inlined from RadiomicsHeterogeneityLogic.py)
# ============================================================================

try:
    from scipy import stats as scipy_stats

    def _scipy_stats_skew(a):
        return float(scipy_stats.skew(a, bias=False))

    def _scipy_stats_kurtosis(a, fisher=True):
        return float(scipy_stats.kurtosis(a, fisher=fisher))

except ImportError:  # pragma: no cover — fallback when scipy is absent

    def _scipy_stats_skew(a):
        n = len(a)
        m = np.mean(a)
        s = np.std(a, ddof=0)
        if s == 0:
            return 0.0
        return float(np.mean(((a - m) / s) ** 3) * n * (n - 1)) / (
            (n - 2) * (n - 3) / (n * (n + 1))
        )

    def _scipy_stats_kurtosis(a, fisher=True):
        n = len(a)
        m = np.mean(a)
        s = np.std(a, ddof=0)
        if s == 0:
            return 0.0
        k = float(np.mean(((a - m) / s) ** 4))
        correction = (3.0 * (n - 1) ** 2) / ((n - 2) * (n - 3)) if n > 3 else 0
        return k - 3 + correction if fisher else k


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def _clip_array(arr: np.ndarray) -> np.ndarray:
    """Cast to int32 to avoid memory issues with very large arrays."""
    return arr.astype(np.int32, copy=False)


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------


def extract_mask_and_intensity(volume_node, segment_node):
    """Extract a 3D intensity array masked by the given segment.

    Parameters
    ----------
    volume_node :
        A ``vtkMRMLScalarVolumeNode`` or ``vtkMRMLVolumeNode``.
    segment_node :
        A ``vtkMRMLSegmentationNode`` whose first segment is used.

    Returns
    -------
    tuple of (mask_3d, intensity_3d, shape, spacing)
        ``mask_3d`` is a boolean array; ``intensity_3d`` is the volume data
        restricted to the same geometry.
    """
    volume_array = slicer.util.arrayFromVolume(volume_node)
    shape = volume_array.shape  # (X, Y, Z)
    spacing = volume_node.GetSpacing()

    # Get binary labelmap array for the first segment
    segmentation = segment_node.GetSegmentation()
    segment_id = segmentation.GetNthSegmentID(0)
    mask_array = slicer.util.arrayFromSegmentBinaryLabelmap(
        segment_node, segment_id, referenceVolumeNode=volume_node
    )

    return mask_array.astype(bool), volume_array.astype(np.float64), shape, spacing




# ---------------------------------------------------------------------------
# First-order features
# ---------------------------------------------------------------------------

FIRST_ORDER_KEYS = [
    "firstorder_Count",
    "firstorder_Mean", "firstorder_StdDev", "firstorder_Variance",
    "firstorder_Median", "firstorder_Minimum", "firstorder_Maximum",
    "firstorder_Range", "firstorder_InterquartileRange",
    "firstorder_P5", "firstorder_P25", "firstorder_P75", "firstorder_P95",
    "firstorder_MeanAbsoluteDeviation",
    "firstorder_Skewness", "firstorder_Kurtosis", "firstorder_Entropy",
    "firstorder_CoV", "firstorder_RobustSpreadEstimator",
]


def compute_first_order(values: np.ndarray) -> dict[str, float]:
    """Compute first-order intensity statistics."""
    if values.size == 0:
        return {k: 0.0 for k in FIRST_ORDER_KEYS}

    n = len(values)
    features: dict[str, float] = {}
    features["firstorder_Count"] = float(n)
    features["firstorder_Mean"] = float(np.mean(values))
    features["firstorder_StdDev"] = float(np.std(values, ddof=0))
    features["firstorder_Variance"] = float(np.var(values, ddof=0))
    features["firstorder_Median"] = float(np.median(values))
    features["firstorder_Minimum"] = float(np.min(values))
    features["firstorder_Maximum"] = float(np.max(values))
    features["firstorder_Range"] = (
        features["firstorder_Maximum"] - features["firstorder_Minimum"]
    )
    features["firstorder_InterquartileRange"] = float(
        np.percentile(values, 75) - np.percentile(values, 25)
    )
    features["firstorder_P5"] = float(np.percentile(values, 5))
    features["firstorder_P25"] = float(np.percentile(values, 25))
    features["firstorder_P75"] = float(np.percentile(values, 75))
    features["firstorder_P95"] = float(np.percentile(values, 95))
    features["firstorder_MeanAbsoluteDeviation"] = float(
        np.median(np.abs(values - np.median(values)))
    )
    features["firstorder_Skewness"] = float(_scipy_stats_skew(values))
    features["firstorder_Kurtosis"] = float(
        _scipy_stats_kurtosis(values, fisher=True)
    )
    features["firstorder_Entropy"] = float(compute_shannon_entropy(values))
    features["firstorder_CoV"] = (
        features["firstorder_StdDev"] / features["firstorder_Mean"]
        if features["firstorder_Mean"] != 0
        else 0.0
    )
    features["firstorder_RobustSpreadEstimator"] = (
        features["firstorder_InterquartileRange"] / features["firstorder_Median"]
        if features["firstorder_Median"] != 0
        else 0.0
    )
    return features


# ---------------------------------------------------------------------------
# Morphological features
# ---------------------------------------------------------------------------


def compute_morphological(shape: tuple, spacing, mask: np.ndarray) -> dict[str, float]:
    """Compute shape-related metrics from the binary mask."""
    voxel_volume = float(np.prod(spacing))
    count = int(np.sum(mask))
    volume_mm3 = count * voxel_volume

    # Bounding box
    nz, ny, nx = np.where(mask)
    if count == 0:
        return {
            "morphometry_VoxelCount": 0,
            "morphometry_Volume_mm3": 0.0,
            "morphometry_BoundingBox_X_mm": 0.0,
            "morphometry_BoundingBox_Y_mm": 0.0,
            "morphometry_BoundingBox_Z_mm": 0.0,
            "morphometry_Compactness": 0.0,
            "morphometry_Sphericity": 0.0,
        }

    bb_x = (nx.max() - nx.min() + 1) * spacing[0]
    bb_y = (ny.max() - ny.min() + 1) * spacing[1]
    bb_z = (nz.max() - nz.min() + 1) * spacing[2]

    # Compactness = surface_area_of_bbox / (6 * volume)  (1 = perfect cube)
    bb_surface = 2 * (bb_x * bb_y + bb_y * bb_z + bb_x * bb_z)
    compactness = bb_surface / (6.0 * volume_mm3) if volume_mm3 > 0 else 0.0

    # Sphericity = (pi^(1/3) * (6*V)^(2/3)) / A_bbox
    # where A_bbox = bb_surface
    sphericity = (
        (np.pi ** (1 / 3)) * (6 * volume_mm3) ** (2 / 3) / bb_surface
        if bb_surface > 0
        else 0.0
    )

    return {
        "morphometry_VoxelCount": count,
        "morphometry_Volume_mm3": volume_mm3,
        "morphometry_BoundingBox_X_mm": bb_x,
        "morphometry_BoundingBox_Y_mm": bb_y,
        "morphometry_BoundingBox_Z_mm": bb_z,
        "morphometry_Compactness": compactness,
        "morphometry_Sphericity": sphericity,
    }


# ---------------------------------------------------------------------------
# GLCM (3D)
# ---------------------------------------------------------------------------

GLCM_KEYS = [
    "glcm_Contrast", "glcm_Correlation", "glcm_Energy",
    "glcm_Homogeneity", "glcm_Entropy", "glcm_IdmVariance",
    "glcm_InformationMeasureCorrelation1", "glcm_InformationMeasureCorrelation2",
]


def _compute_3d_glcms(intensity, mask, num_levels: int) -> list[np.ndarray]:
    """Compute 3D GLCM for each unique direction.

    Returns a list of GLCM matrices (one per direction), each of shape
    (num_levels, num_levels).  Uses only positive shifts since GLCM pairs
    are undirected (we count both (i,j) and (j,i) implicitly by using
    positive offsets only).
    """
    vmin, vmax = intensity[mask].min(), intensity[mask].max()
    if vmax == vmin:
        return [np.ones((1, 1))]

    step = (vmax - vmin) / num_levels
    levels = np.clip(
        np.floor((intensity - vmin) / step).astype(np.int32), 0, num_levels - 1
    )

    # Use only the 7 unique positive-shift directions (face + edge + corner)
    pos_dirs = [
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 1, 0), (1, 0, 1), (0, 1, 1),
        (1, 1, 1),
    ]

    glcms = []
    for dx, dy, dz in pos_dirs:
        # Current voxel region: must leave room for the shifted neighbour
        x1, x2 = dx, intensity.shape[0]
        y1, y2 = dy, intensity.shape[1]
        z1, z2 = dz, intensity.shape[2]

        sub_mask = (
            mask[x1:x2, y1:y2, z1:z2]
            & mask[x1 - dx:x2 - dx, y1 - dy:y2 - dy, z1 - dz:z2 - dz]
        )
        if not np.any(sub_mask):
            glcms.append(np.zeros((num_levels, num_levels)))
            continue

        i0 = levels[x1:x2, y1:y2, z1:z2][sub_mask]
        i1 = levels[x1 - dx:x2 - dx, y1 - dy:y2 - dy, z1 - dz:z2 - dz][sub_mask]

        cm = np.bincount(i0 * num_levels + i1, minlength=num_levels * num_levels)
        glcms.append(cm.reshape(num_levels, num_levels))

    return glcms


def compute_glcM_features(glcms: list[np.ndarray]) -> dict[str, float]:
    """Compute texture features averaged across all 3D directions."""
    if not glcms or all(np.sum(cm) == 0 for cm in glcms):
        return {k: 0.0 for k in GLCM_KEYS}

    # Normalize each GLCM to a probability matrix and average
    n_dirs = len(glcms)
    p_avg = np.zeros_like(glcms[0], dtype=np.float64)
    for cm in glcms:
        s = cm.sum()
        if s > 0:
            p_avg += cm.astype(np.float64) / s
    p_avg /= n_dirs

    n = p_avg.shape[0]
    if n == 0 or p_avg.sum() == 0:
        return {k: 0.0 for k in GLCM_KEYS}

    # Precompute helpers
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    p_ij = p_avg
    p_i = p_ij.sum(axis=1)
    p_j = p_ij.sum(axis=0)

    features: dict[str, float] = {}

    # ASM / Energy
    features["glcm_Energy"] = float(np.sum(p_ij ** 2))

    # Contrast
    features["glcm_Contrast"] = float(np.sum(p_ij * ((i - j) ** 2)))

    # Correlation
    mu_x = float(np.sum(i * p_i))
    mu_y = float(np.sum(j * p_j))
    sigma_x = float(np.sqrt(np.sum(((i - mu_x) ** 2) * p_i)))
    sigma_y = float(np.sqrt(np.sum(((j - mu_y) ** 2) * p_j)))
    if sigma_x * sigma_y > 0:
        features["glcm_Correlation"] = float(
            np.sum(p_ij * (i - mu_x) * (j - mu_y)) / (sigma_x * sigma_y)
        )
    else:
        features["glcm_Correlation"] = 0.0

    # Homogeneity (Inverse Difference Moment)
    features["glcm_Homogeneity"] = float(np.sum(p_ij / (1 + (i - j) ** 2)))

    # Entropy
    nonzero = p_ij > 0
    features["glcm_Entropy"] = float(-np.sum(p_ij[nonzero] * np.log(p_ij[nonzero])))

    # Idm Variance (variance of IDM)
    idm = p_ij / (1 + (i - j) ** 2)
    features["glcm_IdmVariance"] = float(np.var(idm))

    # IMC1 (Informational Measure of Correlation #1)
    hx = float(-np.sum(p_i * np.log(p_i + 1e-12)))
    hy = float(-np.sum(p_j * np.log(p_j + 1e-12)))
    hxy = features["glcm_Entropy"]
    features["glcm_InformationMeasureCorrelation1"] = (
        (hxy - max(hx, hy)) / max(abs(hxy), 1e-12) if hxy > 0 else 0.0
    )

    # IMC2 (Informational Measure of Correlation #2)
    hxy_max = max(hx, hy)
    features["glcm_InformationMeasureCorrelation2"] = (
        (1 - np.exp(-2 * (hxy_max - hxy))) ** 0.5 if hxy_max > hxy else 0.0
    )

    return features


# ---------------------------------------------------------------------------
# GLRLM (3D)
# ---------------------------------------------------------------------------

GLRLM_KEYS = [
    "glrlm_RunLengthEmphasis", "glrlm_ShortRunEmphasis", "glrlm_LongRunEmphasis",
    "glrlm_GrayLevelNonUniformity", "glrlm_RunLengthNonUniformity",
    "glrlm_RunLengthNonUniformityNormalized", "glrlm_ShortRunLargeGrayEmphasis",
    "glrlm_LongRunLargeGrayEmphasis", "glrlm_RunLengthGrayLevelNonUniformity",
    "glrlm_GrayLevelVariance",
]


def _compute_3d_glrlms(intensity, mask, num_levels: int) -> list[np.ndarray]:
    """Compute 3D GLRLM matrices for each unique positive direction.

    Uses the 7 unique positive shifts (face + edge + corner).
    Returns a list of GLRLM matrices (one per direction), shape
    (num_levels, max_run_length).
    """
    vmin, vmax = intensity[mask].min(), intensity[mask].max()
    if vmax == vmin:
        return [np.zeros((1, 1))]

    step = (vmax - vmin) / num_levels
    levels = np.clip(
        np.floor((intensity - vmin) / step).astype(np.int32), 0, num_levels - 1
    )

    # 7 unique positive-shift directions (same set used for GLCM)
    pos_dirs = [
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 1, 0), (1, 0, 1), (0, 1, 1),
        (1, 1, 1),
    ]

    glrlms = []
    max_possible_run = max(intensity.shape)

    for dx, dy, dz in pos_dirs:
        rlm = np.zeros((num_levels, max_possible_run), dtype=np.int64)

        # Shifted mask: prev_mask[x] = mask[x - shift], so a run starts
        # where mask is True but the predecessor (x - shift) is False or
        # out of bounds.
        prev_mask = np.zeros_like(mask)
        if dx > 0:
            prev_mask[dx:, :, :] = mask[:-dx, :, :]
        if dy > 0:
            prev_mask[:, dy:, :] = mask[:, :-dy, :]
        if dz > 0:
            prev_mask[:, :, dz:] = mask[:, :, :-dz]

        # Run starts: in mask, but predecessor is NOT in mask
        run_starts = mask & ~prev_mask

        nz, ny, nx = np.where(run_starts)
        for idx in range(len(nz)):
            x0, y0, z0 = nz[idx], ny[idx], nx[idx]
            gray = levels[x0, y0, z0]

            # Trace the run forward
            run_len = 0
            x, y, z = x0, y0, z0
            while (
                0 <= x < levels.shape[0]
                and 0 <= y < levels.shape[1]
                and 0 <= z < levels.shape[2]
                and mask[x, y, z]
                and levels[x, y, z] == gray
            ):
                run_len += 1
                x += dx
                y += dy
                z += dz
                if run_len >= max_possible_run:
                    break

            if run_len > 0:
                rlm[gray, run_len - 1] += 1

        glrlms.append(rlm)

    return glrlms


def compute_glrlM_features(glrlms: list[np.ndarray]) -> dict[str, float]:
    """Compute GLRLM texture features averaged across directions."""
    if not glrlms or all(np.sum(rl) == 0 for rl in glrlms):
        return {k: 0.0 for k in GLRLM_KEYS}

    n_dirs = len(glrlms)
    # Sum all RLMs and normalize
    rlm_sum = np.zeros_like(glrlms[0], dtype=np.float64)
    for rl in glrlms:
        rlm_sum += rl.astype(np.float64)
    total = rlm_sum.sum()
    if total == 0:
        return {k: 0.0 for k in GLRLM_KEYS}

    p = rlm_sum / total
    n_grays = p.shape[0]  # number of gray levels
    n_runs = p.shape[1]   # max run length
    gray_vals = np.arange(n_grays)
    run_lens = np.arange(1, n_runs + 1)

    features: dict[str, float] = {}

    # RLE (Run Length Emphasis) = sum of l^2 * p(l, g)
    features["glrlm_RunLengthEmphasis"] = float(np.sum(
        p * (run_lens ** 2)[None, :]
    ))

    # SRE (Short Run Emphasis)
    features["glrlm_ShortRunEmphasis"] = float(
        np.sum(p / (run_lens[None, :] ** 2)) / n_grays
    )

    # LRE (Long Run Emphasis)
    features["glrlm_LongRunEmphasis"] = float(
        np.sum((run_lens[None, :] ** 2) * p) / n_grays
    )

    # GLN (Gray Level Non-Uniformity)
    p_g = p.sum(axis=1)
    features["glrlm_GrayLevelNonUniformity"] = float(np.sqrt(np.sum(p_g ** 2)))

    # RLN (Run Length Non-Uniformity)
    p_l = p.sum(axis=0)
    features["glrlm_RunLengthNonUniformity"] = float(np.sqrt(np.sum(p_l ** 2)))

    # RLNNG (Run Length Non-Uniformity Normalized Gray)
    features["glrlm_RunLengthNonUniformityNormalized"] = (
        features["glrlm_RunLengthNonUniformity"] / np.sqrt(n_grays)
        if n_grays > 0
        else 0.0
    )

    # SRLGE (Short Run Large Gray Emphasis)
    features["glrlm_ShortRunLargeGrayEmphasis"] = float(np.sum(
        p * gray_vals[:, None] ** 2 / (run_lens[None, :] ** 2)
    ))

    # LRLGE (Long Run Large Gray Emphasis)
    features["glrlm_LongRunLargeGrayEmphasis"] = float(np.sum(
        p * gray_vals[:, None] ** 2 * (run_lens[None, :] ** 2)
    ))

    # RLGNonU (Run Length Gray-Level Non-Uniformity)
    features["glrlm_RunLengthGrayLevelNonUniformity"] = float(np.sum(p ** 2))

    # GLNonU (Gray-Level Non-Uniformity)
    features["glrlm_GrayLevelVariance"] = float(np.sum(p_g ** 2))

    return features


# ---------------------------------------------------------------------------
# GLZDM (3D)
# ---------------------------------------------------------------------------

GLZDM_KEYS = [
    "gldm_SmallDependenceEmphasis", "gldm_LargeDependenceEmphasis",
    "gldm_DependencePercent", "gldm_DependenceVariance", "gldm_DependenceEntropy",
]


def _compute_3d_glzdm(intensity, mask, num_levels: int) -> np.ndarray:
    """Compute the 3D Gray Level Zone Dependency Matrix.

    Uses the 7 unique positive-shift directions.
    A 'zone' is a connected component of equal gray level voxels in a given
    direction. GLZDM[z, g] = number of zones of length z at gray level g.
    """
    vmin, vmax = intensity[mask].min(), intensity[mask].max()
    if vmax == vmin:
        return np.zeros((1, 1))

    step = (vmax - vmin) / num_levels
    levels = np.clip(
        np.floor((intensity - vmin) / step).astype(np.int32), 0, num_levels - 1
    )

    max_zone = max(intensity.shape)
    glzdm = np.zeros((max_zone, num_levels), dtype=np.int64)

    # 7 unique positive-shift directions
    pos_dirs = [
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (1, 1, 0), (1, 0, 1), (0, 1, 1),
        (1, 1, 1),
    ]

    for dx, dy, dz in pos_dirs:
        # Build the "previous voxel in direction" mask
        prev_mask = np.zeros_like(mask)
        if dx > 0:
            prev_mask[dx:, :, :] = mask[:-dx, :, :]
        if dy > 0:
            prev_mask[:, dy:, :] = mask[:, :-dy, :]
        if dz > 0:
            prev_mask[:, :, dz:] = mask[:, :, :-dz]

        zone_starts = mask & ~prev_mask
        nz, ny, nx = np.where(zone_starts)

        for idx in range(len(nz)):
            x0, y0, z0 = nz[idx], ny[idx], nx[idx]
            gray = levels[x0, y0, z0]
            zone_len = 0
            x, y, z = x0, y0, z0
            while (
                0 <= x < levels.shape[0]
                and 0 <= y < levels.shape[1]
                and 0 <= z < levels.shape[2]
                and mask[x, y, z]
                and levels[x, y, z] == gray
            ):
                zone_len += 1
                x += dx
                y += dy
                z += dz
                if zone_len >= max_zone:
                    break
            if zone_len > 0:
                glzdm[zone_len - 1, gray] += 1

    return glzdm


def compute_glzdm_features(glzdm: np.ndarray) -> dict[str, float]:
    """Compute GLZDM features."""
    if glzdm.size == 0 or glzdm.sum() == 0:
        return {k: 0.0 for k in GLZDM_KEYS}

    total = glzdm.sum()
    p = glzdm.astype(np.float64) / total
    zone_lens = np.arange(1, p.shape[0] + 1)
    gray_vals = np.arange(p.shape[1])

    features: dict[str, float] = {}

    # SZE (Small Dependence Emphasis)
    features["gldm_SmallDependenceEmphasis"] = float(
        np.sum(p / (zone_lens[:, None] ** 2))
    )

    # LZE (Large Dependence Emphasis)
    features["gldm_LargeDependenceEmphasis"] = float(
        np.sum(p * (zone_lens[:, None] ** 2))
    )

    # ZP (Zone Percentage) = 1 - sum(p^2)
    features["gldm_DependencePercent"] = float(1.0 - np.sum(p ** 2))

    # ZV (Zone Variance)
    mean_z = float(np.sum(p * zone_lens[:, None]))
    features["gldm_DependenceVariance"] = float(
        np.sum(p * (zone_lens[:, None] - mean_z) ** 2)
    )

    # ZE (Zone Entropy)
    nonzero = p > 0
    features["gldm_DependenceEntropy"] = float(
        -np.sum(p[nonzero] * np.log(p[nonzero]))
    )

    return features


# ---------------------------------------------------------------------------
# NGTDM
# ---------------------------------------------------------------------------

NGTDM_KEYS = [
    "ngtdm_Coarseness", "ngtdm_Contrast", "ngtdm_Complexity",
    "ngtdm_Strength", "ngtdm_Granularity",
]


def compute_ngtdM_features(intensity, mask, num_levels: int) -> dict[str, float]:
    """Compute NGTDM features.

    For each masked voxel, the 'neighboring gray tone difference' is the
    average gray level of its (face-adjacent) neighbors inside the mask,
    minus the voxel's own gray level. The NGTDM is a histogram of these
    differences.
    """
    vmin, vmax = intensity[mask].min(), intensity[mask].max()
    if vmax == vmin:
        return {k: 0.0 for k in NGTDM_KEYS}

    # Compute face-neighbor sums and counts inside the mask
    local_sum = np.zeros_like(intensity, dtype=np.float64)
    local_count = np.zeros_like(mask, dtype=np.int32)
    for dx, dy, dz in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
        if dx > 0:
            lo, hi = dx, intensity.shape[0]
        elif dx < 0:
            lo, hi = 0, intensity.shape[0] + dx
        else:
            lo, hi = 0, intensity.shape[0]

        if dy > 0:
            lo2, hi2 = dy, intensity.shape[1]
        elif dy < 0:
            lo2, hi2 = 0, intensity.shape[1] + dy
        else:
            lo2, hi2 = 0, intensity.shape[1]

        if dz > 0:
            lo3, hi3 = dz, intensity.shape[2]
        elif dz < 0:
            lo3, hi3 = 0, intensity.shape[2] + dz
        else:
            lo3, hi3 = 0, intensity.shape[2]

        local_sum[lo:hi, lo2:hi2, lo3:hi3] += intensity[
            lo - dx:hi - dx, lo2 - dy:hi2 - dy, lo3 - dz:hi3 - dz
        ]
        local_count[lo:hi, lo2:hi2, lo3:hi3] += 1

    # Neighborhood gray tone
    with np.errstate(divide="ignore", invalid="ignore"):
        ngt = np.where(
            local_count > 0,
            local_sum / local_count - intensity,
            0.0,
        )

    # Quantize the difference
    ngt_flat = ngt[mask]
    if ngt_flat.size == 0:
        return {k: 0.0 for k in NGTDM_KEYS}

    # Create histogram of NGT values (quantized)
    n_bins = 100
    ngt_min, ngt_max = ngt_flat.min(), ngt_flat.max()
    if ngt_max == ngt_min:
        return {k: 0.0 for k in NGTDM_KEYS}

    # Bin the NGT values
    ngt_quantized = np.clip(
        np.floor(
            (ngt_flat - ngt_min) / (ngt_max - ngt_min) * (n_bins - 1)
        ).astype(np.int32),
        0, n_bins - 1,
    )
    hist, _ = np.histogram(ngt_quantized, bins=n_bins)
    p_ngt = hist / hist.sum() if hist.sum() > 0 else hist

    features: dict[str, float] = {}
    bin_centers = (np.arange(n_bins) + 0.5) * (ngt_max - ngt_min) / n_bins

    # Coarseness = 1 / sum(p * |bin_center|)
    abs_centers = np.abs(bin_centers)
    denom = np.sum(p_ngt * abs_centers)
    features["ngtdm_Coarseness"] = 1.0 / denom if denom > 0 else 0.0

    # Contrast = sum(p * |bin_center|^2)
    features["ngtdm_Contrast"] = float(np.sum(p_ngt * bin_centers ** 2))

    # Complexity
    mean_ngt = float(np.sum(p_ngt * bin_centers))
    var_ngt = float(np.sum(p_ngt * (bin_centers - mean_ngt) ** 2))
    features["ngtdm_Complexity"] = float(
        np.sum(p_ngt * np.abs(bin_centers) * np.sqrt(var_ngt + 1e-12))
    )

    # Strength
    features["ngtdm_Strength"] = float(
        np.sum(p_ngt * np.sqrt(abs_centers + 1e-12))
    )

    # Granularity
    features["ngtdm_Granularity"] = (
        float(hist.max() / hist.sum()) if hist.sum() > 0 else 0.0
    )

    return features


# ---------------------------------------------------------------------------
# Shannon entropy helper
# ---------------------------------------------------------------------------


def compute_shannon_entropy(values: np.ndarray, num_bins: int = 256) -> float:
    """Shannon entropy of the intensity distribution."""
    if values.size == 0:
        return 0.0
    hist, _ = np.histogram(values, bins=num_bins)
    p = hist / hist.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def compute_all_features(
    volume_node,
    segment_node,
    num_levels: int = 256,
) -> dict:
    """Compute the full feature set for a volume/segment pair.

    Returns a flat dictionary of feature name -> value suitable for a table.
    """
    import vtk

    mask, intensity, shape, spacing = extract_mask_and_intensity(
        volume_node, segment_node
    )

    if not np.any(mask):
        return {}

    # Subset to masked voxels
    values = intensity[mask]

    # 1. First-order
    features = compute_first_order(values)

    # 2. Morphological
    features.update(compute_morphological(shape, spacing, mask))

    # 3. 3D GLCM
    glcms = _compute_3d_glcms(intensity, mask, num_levels)
    features.update(compute_glcM_features(glcms))

    # 4. 3D GLRLM
    glrlms = _compute_3d_glrlms(intensity, mask, num_levels)
    features.update(compute_glrlM_features(glrlms))

    # 5. 3D GLZDM
    glzdm = _compute_3d_glzdm(intensity, mask, num_levels)
    features.update(compute_glzdm_features(glzdm))

    # 6. NGTDM
    features.update(compute_ngtdM_features(intensity, mask, num_levels))

    return features


# ============================================================================
# Logic wrapper (provides a class interface compatible with the previous
# RadiomicsHeterogeneityLogic.py)
# ============================================================================


class RadiomicsHeterogeneityLogic:
    """Thin wrapper that delegates to :func:`compute_all_features`."""

    def compute_all_features(self, volume_node, segment_node, num_levels=256):
        return compute_all_features(
            volume_node, segment_node, num_levels=num_levels
        )


# ============================================================================
# Slicer module class
# ============================================================================


log = logging.getLogger(__name__)


class RadiomicsHeterogeneity(ScriptedLoadableModule):
    """radiomics heterogeneity analysis."""

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "Radiomics Heterogeneity"
        parent.categories = ["Quantification"]
        parent.dependencies = []
        parent.contributors = ["Your Name"]
        parent.helpText = self._help_text()
        parent.acknowledgementText = (
            "This module was developed for volumetric radiomics analysis."
        )

    # ── helpText ─────────────────────────────────────────────────

    @staticmethod
    def _help_text() -> str:
        return """
        <h2>Radiomics Heterogeneity</h2>
        <p>
        Computes 3D radiomics heterogeneity features from a CT or PET volume
        using a segmentation as the region-of-interest.
        </p>

        <h3>How to use</h3>
        <ol>
          <li>Load a CT or PET volume into the scene.</li>
          <li>Create or load a segmentation and make sure at least one segment
              is defined (use the <em>Segment Editor</em> module to draw one).</li>
          <li>In this module's panel, select the <b>Volume</b> and the
              <b>Segmentation</b> from the drop-down lists.</li>
          <li>Optionally change the <b>Gray Levels</b> (number of bins used to
              quantize the intensity for texture features — choose from
              8, 16, 32, or 64).</li>
          <li>Click <b>Compute Features</b>.</li>
          <li>Results appear in the <em>Results</em> table. Use
              <b>Export CSV</b> to save them to disk.</li>
        </ol>

        <h3>Feature categories</h3>
        <ul>
          <li><b>First-order statistics</b> — shape of the intensity histogram
              (mean, std, skewness, kurtosis, percentiles, entropy, CoV, etc.).</li>
          <li><b>3D GLCM</b> — Gray Level Co-occurrence Matrix features
              averaged over all 13 unique 3D neighbour directions.</li>
          <li><b>3D GLRLM</b> — Gray Level Run Length Matrix features.</li>
          <li><b>3D GLZDM</b> — Gray Level Zone Dependency Matrix features.</li>
          <li><b>NGTDM</b> — Neighboring Gray Tone Difference Matrix features.</li>
          <li><b>Morphological</b> — voxel count, volume, bounding-box size,
              compactness, sphericity.</li>
        </ul>

        <h3>Notes</h3>
        <ul>
          <li>All features are computed over the <b>3D</b> volume defined by
              the selected segment — not slice-by-slice.</li>
          <li>Texture features (GLCM, GLRLM, GLZDM) are averaged over all
              unique 3D directions for rotation-invariant results.</li>
          <li>Intensities are quantized to <i>N</i> gray levels before
              texture computation. Lower values speed up computation and
              reduce noise sensitivity.</li>
          <li>For large volumes the computation may take a few seconds.</li>
        </ul>
        """


# ============================================================================
# Widget class (inlined from RadiomicsHeterogeneityGUI.py)
# Slicer discovers this by the <ModuleName>Widget naming convention.
# ============================================================================


class RadiomicsHeterogeneityWidget(ScriptedLoadableModuleWidget):
    """Widget for the Radiomics Heterogeneity module."""

    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        self.logic = None
        self._volumeNodes: list = []
        self._segmentationNodes: list = []
        self._results: dict[str, float] = {}

    # ── setup ─────────────────────────────────────────────────────

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        # Register the logic
        self.logic = RadiomicsHeterogeneityLogic()

        # ── Layout ──────────────────────────────────────────────
        self.layout = self.parent.layout()
        self.layout.setContentsMargins(4, 4, 4, 4)

        # Collapsible parameter panel
        self.paramsCollapsibleButton = ctk.ctkCollapsibleButton()
        self.paramsCollapsibleButton.text = "Parameters"
        self.paramsCollapsibleButton.setMaximumWidth(350)
        self.layout.addWidget(self.paramsCollapsibleButton)
        self.paramsFormLayout = qt.QFormLayout(self.paramsCollapsibleButton)

        # Volume selector
        self.volumeSelector = qt.QComboBox()
        self.volumeSelector.connect(
            "currentIndexChanged(int)", self._onVolumeChanged
        )
        self.paramsFormLayout.addRow("Volume:", self.volumeSelector)

        # Segmentation selector
        self.segmentationSelector = qt.QComboBox()
        self.segmentationSelector.connect(
            "currentIndexChanged(int)", self._onSegmentationChanged
        )
        self.paramsFormLayout.addRow("Segmentation:", self.segmentationSelector)

        # Segment selector
        self.segmentSelector = qt.QComboBox()
        self.segmentSelector.connect(
            "currentIndexChanged(int)", self._onSegmentChanged
        )
        self.paramsFormLayout.addRow("Segment:", self.segmentSelector)

        # Gray levels bin count (dropdown)
        self.grayLevelsComboBox = qt.QComboBox()
        self.grayLevelsComboBox.addItems(["8", "16", "32", "64"])
        self.grayLevelsComboBox.setCurrentIndex(3)  # default: 64
        self.paramsFormLayout.addRow("Gray Levels:", self.grayLevelsComboBox)

        # Compute button
        self.computeButton = qt.QPushButton("Compute Features")
        self.computeButton.enabled = True
        self.computeButton.connect("clicked()", self._onCompute)
        self.paramsFormLayout.addRow(self.computeButton)

        # Progress bar
        self.progressBar = qt.QProgressBar()
        self.progressBar.hide()
        self.progressBar.minimum = 0
        self.progressBar.maximum = 6
        self.progressBar.value = 0
        self.layout.addWidget(self.progressBar)

        # ── Results section ─────────────────────────────────────
        self.resultsCollapsibleButton = ctk.ctkCollapsibleButton()
        self.resultsCollapsibleButton.text = "Results"
        self.resultsCollapsibleButton.setMaximumWidth(350)
        self.layout.addWidget(self.resultsCollapsibleButton)
        self.resultsFormLayout = qt.QFormLayout(self.resultsCollapsibleButton)

        self.resultsTable = qt.QTableWidget()
        self.resultsTable.maximumHeight = 400
        self.resultsTable.horizontalHeader().stretchLastSection = True
        self.resultsFormLayout.addRow(self.resultsTable)

        # Export button
        self.exportButton = qt.QPushButton("Export CSV")
        self.exportButton.enabled = False
        self.exportButton.connect("clicked()", self._onExport)
        self.resultsFormLayout.addRow(self.exportButton)

        # ── Status ──────────────────────────────────────────────
        self.statusLabel = qt.QLabel("")
        self.statusLabel.maximumWidth = 350
        self.layout.addWidget(self.statusLabel)

        # Populate selectors from current scene
        self._populateSelectors()

        # Watch for changes in the scene
        scene = slicer.mrmlScene
        scene.AddObserver(
            slicer.vtkMRMLScene.NodeAddedEvent, self._onSceneChanged
        )
        scene.AddObserver(
            slicer.vtkMRMLScene.NodeRemovedEvent, self._onSceneChanged
        )

    def _populateSelectors(self):
        """Fill volume and segmentation combo boxes from the current scene."""
        vol_nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        seg_nodes = slicer.util.getNodesByClass("vtkMRMLSegmentationNode")

        self.populateVolumeSelector(vol_nodes)
        self.populateSegmentationSelector(seg_nodes)

    def _onSceneChanged(self, _caller=None, _event=None):
        """Re-populate selectors whenever a node is added or removed."""
        self._populateSelectors()

    # ── enter / exit ──────────────────────────────────────────────

    def enter(self):
        pass

    def exit(self):
        pass

    # ── Population helpers ────────────────────────────────────────

    def populateVolumeSelector(self, nodes):
        """Fill the volume combo box from a list of ``vtkMRMLScalarVolumeNode``."""
        self._volumeNodes = nodes
        self.volumeSelector.clear()
        for node in nodes:
            self.volumeSelector.addItem(node.GetName())
        if len(nodes) > 0:
            self.volumeSelector.setCurrentIndex(0)

    def populateSegmentationSelector(self, nodes):
        """Fill the segmentation combo box from a list of ``vtkMRMLSegmentationNode``."""
        self._segmentationNodes = nodes
        self.segmentationSelector.clear()
        for node in nodes:
            self.segmentationSelector.addItem(node.GetName())
        if len(nodes) > 0:
            self.segmentationSelector.setCurrentIndex(0)
            self._populateSegmentSelector()

    def _populateSegmentSelector(self):
        """Fill the segment combo box from the selected segmentation's segments."""
        if not self._segmentationNodes:
            self.segmentSelector.clear()
            return

        seg_node = self._segmentationNodes[self.segmentationSelector.currentIndex]
        segmentation = seg_node.GetSegmentation()
        num_segments = segmentation.GetNumberOfSegments()

        self.segmentSelector.clear()
        for i in range(num_segments):
            segment = segmentation.GetNthSegment(i)
            self.segmentSelector.addItem(segment.GetName() if segment.GetName() else f"Segment {i+1}")

        if num_segments > 0:
            self.segmentSelector.setCurrentIndex(0)

    # ── Callbacks ─────────────────────────────────────────────────

    def _onVolumeChanged(self):
        self.statusLabel.text = "Volume changed."

    def _onSegmentationChanged(self):
        self.statusLabel.text = "Segmentation changed."
        self._populateSegmentSelector()

    def _onSegmentChanged(self):
        self.statusLabel.text = "Segment changed."

    def _onCompute(self):
        self.statusLabel.text = "Computing features..."
        self.progressBar.show()
        self.progressBar.value = 0
        self.computeButton.enabled = False

        volume_node = self._volumeNodes[self.volumeSelector.currentIndex]
        seg_node = self._segmentationNodes[self.segmentationSelector.currentIndex]
        num_levels = int(self.grayLevelsComboBox.currentText)

        # Get the selected segment ID
        segmentation = seg_node.GetSegmentation()
        segment_id = segmentation.GetNthSegmentID(
            self.segmentSelector.currentIndex
        )

        try:
            # Step 1: First-order
            self.progressBar.value = 1
            mask, intensity, shape, spacing = extract_mask_and_intensity(
                volume_node, seg_node
            )

            if not np.any(mask):
                self.statusLabel.text = "No voxels found in segment."
                return

            values = intensity[mask]

            features = compute_first_order(values)

            # Step 2: Morphological
            self.progressBar.value = 2
            features.update(compute_morphological(shape, spacing, mask))

            # Step 3: 3D GLCM
            self.progressBar.value = 3
            glcms = _compute_3d_glcms(intensity, mask, num_levels)
            features.update(compute_glcM_features(glcms))

            # Step 4: 3D GLRLM
            glrlms = _compute_3d_glrlms(intensity, mask, num_levels)
            features.update(compute_glrlM_features(glrlms))

            # Step 5: 3D GLZDM
            self.progressBar.value = 4
            glzdm = _compute_3d_glzdm(intensity, mask, num_levels)
            features.update(compute_glzdm_features(glzdm))

            # Step 6: NGTDM
            self.progressBar.value = 5
            features.update(compute_ngtdM_features(intensity, mask, num_levels))

            self.progressBar.value = 6
            self._display_results(features)
            self.statusLabel.text = (
                f"Done — {len(features)} features computed."
            )
        except Exception as exc:
            self.statusLabel.text = f"Error: {exc}"
            import traceback
            traceback.print_exc()
        finally:
            self.progressBar.hide()
            self.computeButton.enabled = True

    def _display_results(self, features: dict[str, float]):
        self._results = features
        self.exportButton.enabled = True
        self.resultsTable.clear()
        self.resultsTable.setColumnCount(2)
        self.resultsTable.verticalHeader().setVisible(False)
        self.resultsTable.setHorizontalHeaderLabels(["Feature", "Value"])

        row = 0
        for key, val in features.items():
            self.resultsTable.insertRow(row)
            name_item = qt.QTableWidgetItem(key)
            name_item.setFlags(name_item.flags() & ~qt.Qt.ItemIsEditable)
            self.resultsTable.setItem(row, 0, name_item)

            if isinstance(val, float):
                text = f"{val:.6g}"
            else:
                text = str(val)
            value_item = qt.QTableWidgetItem(text)
            value_item.setFlags(value_item.flags() & ~qt.Qt.ItemIsEditable)
            self.resultsTable.setItem(row, 1, value_item)
            row += 1

        self.resultsTable.resizeColumnsToContents()

    def _onExport(self):
        if not self._results:
            return

        # Ask user where to save, defaulting to Desktop with a filename
        import os
        default_dir = os.path.expanduser("~/Desktop")
        if not os.path.isdir(default_dir):
            default_dir = os.path.expanduser("~")

        path = qt.QFileDialog.getSaveFileName(
            None,
            "Export Features",
            os.path.join(default_dir, "radiomics_heterogeneity.csv"),
            "CSV Files (*.csv)",
        )
        if not path or not path[0]:
            return

        filepath = path[0]

        with open(path, "w") as f:
            f.write("Feature,Value\n")
            for key, val in self._results.items():
                if isinstance(val, float):
                    f.write(f"{key},{val:.10g}\n")
                else:
                    f.write(f"{key},{val}\n")

        self.statusLabel.text = f"Saved to {path}"

"""
PAD Analysis - 3D Slicer Scripted Module
=========================================
Splits a manually-drawn whole-hindlimb CT segmentation (left and/or right leg)
into three anatomical zones - Proximal, Distal, and Foot - using interactively
positioned cut-off planes (upper thigh/hip boundary, knee joint, ankle joint),
then quantifies PET tracer uptake in each resulting zone via SegmentStatistics.

Typical workflow:
  1. In Segment Editor, create a segmentation with a "Left Hindlimb" and/or
     "Right Hindlimb" segment covering the whole leg on the CT image.
  2. In this module, create the three cut-off planes per leg and drag them
     into position directly on the CT (hip/upper-thigh boundary, knee, ankle).
  3. Click Run: each leg's whole-limb segment is split into
     <Side>_Proximal / <Side>_Distal / <Side>_Foot segments, and PET uptake
     statistics are computed for each and exported to CSV.

Author: Generated for 3D Slicer (Python scripted module)
Slicer version: 5.x+
"""

import os
import csv
import logging

import numpy as np
import qt
import ctk
import vtk
import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleWidget,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
)
from slicer.util import VTKObservationMixin


# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

SIDES = [("L", "Left"), ("R", "Right")]
PLANE_ROLES = [
    ("proximal", "Proximal Boundary (hip / upper thigh)"),
    ("knee", "Knee Joint"),
    ("ankle", "Ankle / Foot Joint"),
]
# Order defines the zone-counting logic in PADAnalysisLogic.run() - do not reorder
# without updating the classification loop.
ZONE_NAMES = ["Proximal", "Distal", "Foot"]

PLANE_COLORS = {
    "proximal": (1.00, 0.60, 0.00),   # orange
    "knee":     (0.10, 0.90, 0.10),   # green
    "ankle":    (0.20, 0.60, 1.00),   # blue
}
ZONE_COLORS = {
    "Proximal": (0.20, 0.60, 1.00),   # blue
    "Distal":   (0.20, 0.80, 0.20),   # green
    "Foot":     (1.00, 0.55, 0.00),   # orange
}


# ---------------------------------------------------------------------------
# Module descriptor
# ---------------------------------------------------------------------------

class PADAnalysis(ScriptedLoadableModule):
    """Module metadata shown in the Slicer module finder."""

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "PAD Analysis"
        self.parent.categories = ["Quantification"]
        self.parent.dependencies = ["Markups", "SegmentStatistics"]
        self.parent.contributors = ["Wawrzyniec Dobrucki <br> Experimental Molecular Imaging Laboratory <br> University of Illinois at Urbana-Champaign"]
        self.parent.helpText = (
            "Splits a whole-hindlimb CT segmentation into <b>Proximal</b>, "
            "<b>Distal</b>, and <b>Foot</b> zones using three manually "
            "positioned cut-off planes per leg, then quantifies PET tracer "
            "uptake in each zone."
            "<br><br><b>Before running this module:</b> use Segment Editor to "
            "create a segmentation with a <i>Left Hindlimb</i> and/or "
            "<i>Right Hindlimb</i> segment that covers the whole leg on the "
            "CT image (thresholding + paint/scissors works well)."
            "<br><br><b>Workflow:</b>"
            "<ol>"
            "<li>Select the CT volume, PET volume, and whole-limb segmentation.</li>"
            "<li>Pick which segment corresponds to the left / right leg.</li>"
            "<li>For each leg, click <b>Create</b> next to each of the three "
            "boundary planes, then drag the plane's handles directly on the "
            "CT image to position it at the hip/upper-thigh boundary, the "
            "knee joint, and the ankle joint.</li>"
            "<li>Make sure each plane's normal (the arrow handle) points "
            "toward the foot. Use <b>Flip Normal</b> if it points toward "
            "the body instead.</li>"
            "<li>Click <b>Run</b> to generate the six zone segments and "
            "export a CSV of PET uptake statistics per zone.</li>"
            "</ol>"
            "<br><b>Tissue proximal to the hip/upper-thigh plane is excluded "
            "entirely</b> (not assigned to any zone) so that body-trunk "
            "signal does not contaminate the Proximal hindlimb measurement."
        )
        self.parent.acknowledgementText = (
            "<p>Developed by:</p>"
            "<p><b>Wawrzyniec Dobrucki</b><br>"
            "Experimental Molecular Imaging Laboratory (EMIL)<br>"
            "University of Illinois at Urbana-Champaign</p>"
        )


# ---------------------------------------------------------------------------
# Widget (GUI)
# ---------------------------------------------------------------------------

class PADAnalysisWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """User interface for the module."""

    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None

        # side -> role -> vtkMRMLMarkupsPlaneNode | None
        self.planeNodes = {s: {r: None for r, _ in PLANE_ROLES} for s, _ in SIDES}
        # side -> role -> {"create": btn, "jump": btn, "delete": btn}
        self.planeButtons = {s: {} for s, _ in SIDES}
        # side -> QComboBox (leg segment picker)
        self.legSegmentCombos = {}

    # ------------------------------------------------------------------
    # setup
    # ------------------------------------------------------------------

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = PADAnalysisLogic()

        # ── Input section ─────────────────────────────────────────────
        inputCollapsible = ctk.ctkCollapsibleButton()
        inputCollapsible.text = "Inputs"
        self.layout.addWidget(inputCollapsible)
        inputForm = qt.QFormLayout(inputCollapsible)

        self.ctVolumeSelector = self._makeVolumeSelector(
            "Select the CT volume used to draw the whole-limb segmentation "
            "and to define the geometry the cut planes are classified against."
        )
        inputForm.addRow("CT Volume:", self.ctVolumeSelector)

        self.petVolumeSelector = self._makeVolumeSelector(
            "Select the PET volume to quantify uptake from."
        )
        inputForm.addRow("PET Volume:", self.petVolumeSelector)

        self.wholeLimbSegSelector = slicer.qMRMLNodeComboBox()
        self.wholeLimbSegSelector.nodeTypes = ["vtkMRMLSegmentationNode"]
        self.wholeLimbSegSelector.selectNodeUponCreation = False
        self.wholeLimbSegSelector.addEnabled = False
        self.wholeLimbSegSelector.removeEnabled = False
        self.wholeLimbSegSelector.noneEnabled = True
        self.wholeLimbSegSelector.showHidden = False
        self.wholeLimbSegSelector.setMRMLScene(slicer.mrmlScene)
        self.wholeLimbSegSelector.setToolTip(
            "Segmentation containing the whole (undivided) left and/or right "
            "hindlimb segment(s), created beforehand in Segment Editor."
        )
        inputForm.addRow("Whole-Limb Segmentation:", self.wholeLimbSegSelector)

        for sideKey, sideLabel in SIDES:
            combo = qt.QComboBox()
            combo.addItem("-- None --", None)
            combo.setToolTip(
                f"Segment in the Whole-Limb Segmentation that corresponds to "
                f"the {sideLabel.lower()} hindlimb."
            )
            self.legSegmentCombos[sideKey] = combo
            inputForm.addRow(f"{sideLabel} Leg Segment:", combo)

        # ── Cut planes section ───────────────────────────────────────
        planesCollapsible = ctk.ctkCollapsibleButton()
        planesCollapsible.text = "Cut-Off Planes"
        self.layout.addWidget(planesCollapsible)
        planesLayout = qt.QVBoxLayout(planesCollapsible)

        note = qt.QLabel(
            "Create each plane, then drag its handles directly on the CT "
            "image. The plane's normal arrow should point toward the foot."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color:#555;font-style:italic;")
        planesLayout.addWidget(note)

        for sideKey, sideLabel in SIDES:
            group = qt.QGroupBox(f"{sideLabel} Leg")
            groupLayout = qt.QVBoxLayout(group)
            for roleKey, roleLabel in PLANE_ROLES:
                groupLayout.addLayout(
                    self._buildPlaneRow(sideKey, roleKey, roleLabel)
                )
            planesLayout.addWidget(group)

        # ── Output section ───────────────────────────────────────────
        outputCollapsible = ctk.ctkCollapsibleButton()
        outputCollapsible.text = "Output"
        self.layout.addWidget(outputCollapsible)
        outputForm = qt.QFormLayout(outputCollapsible)

        self.outputSegSelector = slicer.qMRMLNodeComboBox()
        self.outputSegSelector.nodeTypes = ["vtkMRMLSegmentationNode"]
        self.outputSegSelector.selectNodeUponCreation = True
        self.outputSegSelector.addEnabled = True
        self.outputSegSelector.removeEnabled = True
        self.outputSegSelector.renameEnabled = True
        self.outputSegSelector.noneEnabled = False
        self.outputSegSelector.baseName = "PADAnalysis_Segmentation"
        self.outputSegSelector.setMRMLScene(slicer.mrmlScene)
        self.outputSegSelector.setToolTip(
            "Segmentation node that will receive the six zone segments "
            "(existing zone segments with the same name are replaced)."
        )
        outputForm.addRow("Output Segmentation:", self.outputSegSelector)

        statsRow = qt.QHBoxLayout()
        self.cbVoxelCount = self._makeCheckBox("Voxel Count", True)
        self.cbVolume     = self._makeCheckBox("Volume (mm³)", True)
        self.cbMin        = self._makeCheckBox("Min", True)
        self.cbMax        = self._makeCheckBox("Max", True)
        self.cbMean       = self._makeCheckBox("Mean", True)
        self.cbStd        = self._makeCheckBox("StdDev", True)
        self.cbMedian     = self._makeCheckBox("Median", False)
        for cb in (self.cbVoxelCount, self.cbVolume, self.cbMin,
                   self.cbMax, self.cbMean, self.cbStd, self.cbMedian):
            statsRow.addWidget(cb)
        outputForm.addRow("PET Statistics:", statsRow)

        fileRow = qt.QHBoxLayout()
        self.outputCsvEdit = qt.QLineEdit()
        self.outputCsvEdit.setPlaceholderText("Choose output CSV file path …")
        browseBtn = qt.QPushButton("Browse…")
        browseBtn.clicked.connect(self.onBrowseOutput)
        fileRow.addWidget(self.outputCsvEdit)
        fileRow.addWidget(browseBtn)
        outputForm.addRow("Output CSV:", fileRow)

        # ── Run button ────────────────────────────────────────────────
        self.runButton = qt.QPushButton("▶  Run PAD Analysis")
        self.runButton.toolTip = "Split whole-limb segments into zones and quantify PET uptake."
        self.runButton.enabled = False
        self.runButton.setStyleSheet(
            "QPushButton { background-color: #2ecc71; color: white; "
            "font-weight: bold; padding: 8px; border-radius: 4px; }"
            "QPushButton:disabled { background-color: #95a5a6; }"
            "QPushButton:hover { background-color: #27ae60; }"
        )
        self.layout.addWidget(self.runButton)

        self.statusLabel = qt.QLabel("")
        self.statusLabel.setWordWrap(True)
        self.statusLabel.setStyleSheet("color: #555; font-style: italic;")
        self.layout.addWidget(self.statusLabel)

        self.layout.addStretch()

        # Connections
        self.ctVolumeSelector.currentNodeChanged.connect(self._onInputChanged)
        self.petVolumeSelector.currentNodeChanged.connect(self._onInputChanged)
        self.wholeLimbSegSelector.currentNodeChanged.connect(self._onWholeLimbSegChanged)
        for sideKey, _ in SIDES:
            self.legSegmentCombos[sideKey].currentIndexChanged.connect(self._onInputChanged)
        self.outputCsvEdit.textChanged.connect(self._onInputChanged)
        self.runButton.clicked.connect(self.onRunClicked)

        self._onWholeLimbSegChanged(None)

    # ------------------------------------------------------------------
    # UI builder helpers
    # ------------------------------------------------------------------

    def _makeVolumeSelector(self, tooltip):
        selector = slicer.qMRMLNodeComboBox()
        selector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        selector.selectNodeUponCreation = True
        selector.addEnabled = False
        selector.removeEnabled = False
        selector.noneEnabled = True
        selector.showHidden = False
        selector.setMRMLScene(slicer.mrmlScene)
        selector.setToolTip(tooltip)
        return selector

    def _makeCheckBox(self, label, default):
        cb = qt.QCheckBox(label)
        cb.checked = default
        return cb

    def _buildPlaneRow(self, sideKey, roleKey, roleLabel):
        row = qt.QHBoxLayout()
        label = qt.QLabel(roleLabel)
        label.setMinimumWidth(220)
        row.addWidget(label)

        createBtn = qt.QPushButton("Create")
        jumpBtn = qt.QPushButton("Jump To")
        deleteBtn = qt.QPushButton("Delete")
        flipBtn = qt.QPushButton("Flip Normal")
        jumpBtn.enabled = False
        deleteBtn.enabled = False
        flipBtn.enabled = False
        flipBtn.setToolTip(
            "Reverse the plane's normal arrow. The normal should point "
            "toward the foot (distal direction) for correct classification."
        )

        createBtn.clicked.connect(
            lambda checked=False, s=sideKey, r=roleKey: self._onCreatePlane(s, r)
        )
        jumpBtn.clicked.connect(
            lambda checked=False, s=sideKey, r=roleKey: self._onJumpToPlane(s, r)
        )
        deleteBtn.clicked.connect(
            lambda checked=False, s=sideKey, r=roleKey: self._onDeletePlane(s, r)
        )
        flipBtn.clicked.connect(
            lambda checked=False, s=sideKey, r=roleKey: self._onFlipPlane(s, r)
        )

        row.addWidget(createBtn)
        row.addWidget(jumpBtn)
        row.addWidget(deleteBtn)
        row.addWidget(flipBtn)

        self.planeButtons[sideKey][roleKey] = {
            "create": createBtn, "jump": jumpBtn, "delete": deleteBtn, "flip": flipBtn,
        }
        return row

    # ------------------------------------------------------------------
    # Slots - inputs
    # ------------------------------------------------------------------

    def _onWholeLimbSegChanged(self, node):
        segNode = self.wholeLimbSegSelector.currentNode()
        for sideKey, _ in SIDES:
            combo = self.legSegmentCombos[sideKey]
            previous = combo.currentText
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("-- None --", None)
            if segNode is not None:
                segmentation = segNode.GetSegmentation()
                for i in range(segmentation.GetNumberOfSegments()):
                    segId = segmentation.GetNthSegmentID(i)
                    segName = segmentation.GetSegment(segId).GetName()
                    combo.addItem(segName, segId)
            idx = combo.findText(previous)
            combo.setCurrentIndex(idx if idx >= 0 else 0)
            combo.blockSignals(False)
        self._onInputChanged()

    def _onInputChanged(self, *args):
        anySideReady = False
        for sideKey, _ in SIDES:
            hasSegment = self.legSegmentCombos[sideKey].currentData is not None
            hasAllPlanes = all(
                self.planeNodes[sideKey][r] is not None for r, _ in PLANE_ROLES
            )
            if hasSegment and hasAllPlanes:
                anySideReady = True

        ready = (
            self.ctVolumeSelector.currentNode() is not None
            and self.petVolumeSelector.currentNode() is not None
            and self.wholeLimbSegSelector.currentNode() is not None
            and anySideReady
            and self.outputCsvEdit.text.strip() != ""
        )
        self.runButton.enabled = ready

    def onBrowseOutput(self):
        path = qt.QFileDialog.getSaveFileName(
            self.parent, "Save Output CSV", "", "CSV (*.csv)"
        )
        if path:
            if not path.lower().endswith(".csv"):
                path += ".csv"
            self.outputCsvEdit.text = path

    # ------------------------------------------------------------------
    # Slots - planes
    # ------------------------------------------------------------------

    def _onCreatePlane(self, sideKey, roleKey):
        ctVolume = self.ctVolumeSelector.currentNode()
        if ctVolume is None:
            slicer.util.errorDisplay("Select a CT volume first.")
            return

        sideLabel = dict(SIDES)[sideKey]
        roleLabel = dict(PLANE_ROLES)[roleKey]
        nodeName = f"PAD_{sideLabel}_{roleLabel.split(' ')[0]}"

        existing = self.planeNodes[sideKey][roleKey]
        if existing is not None:
            slicer.mrmlScene.RemoveNode(existing)

        planeNode = self.logic.createDefaultPlane(
            ctVolume, sideKey, roleKey, nodeName, PLANE_COLORS[roleKey]
        )
        self.planeNodes[sideKey][roleKey] = planeNode

        buttons = self.planeButtons[sideKey][roleKey]
        buttons["jump"].enabled = True
        buttons["delete"].enabled = True
        buttons["flip"].enabled = True
        buttons["create"].setText("Recreate")

        self._onJumpToPlane(sideKey, roleKey)
        self._onInputChanged()

    def _onJumpToPlane(self, sideKey, roleKey):
        planeNode = self.planeNodes[sideKey][roleKey]
        if planeNode is None or planeNode.GetNumberOfControlPoints() == 0:
            return
        markupsLogic = slicer.modules.markups.logic()
        markupsLogic.JumpSlicesToNthPointInMarkup(planeNode.GetID(), 0, True)

    def _onFlipPlane(self, sideKey, roleKey):
        planeNode = self.planeNodes[sideKey][roleKey]
        if planeNode is None:
            return
        normal = [0.0, 0.0, 0.0]
        planeNode.GetNormalWorld(normal)
        planeNode.SetNormalWorld([-n for n in normal])

    def _onDeletePlane(self, sideKey, roleKey):
        planeNode = self.planeNodes[sideKey][roleKey]
        if planeNode is not None:
            slicer.mrmlScene.RemoveNode(planeNode)
        self.planeNodes[sideKey][roleKey] = None

        buttons = self.planeButtons[sideKey][roleKey]
        buttons["jump"].enabled = False
        buttons["delete"].enabled = False
        buttons["flip"].enabled = False
        buttons["create"].setText("Create")

        self._onInputChanged()

    # ------------------------------------------------------------------
    # Slots - run
    # ------------------------------------------------------------------

    def onRunClicked(self):
        ctVolume = self.ctVolumeSelector.currentNode()
        petVolume = self.petVolumeSelector.currentNode()
        wholeLimbSeg = self.wholeLimbSegSelector.currentNode()
        outputSeg = self.outputSegSelector.currentNode()
        outputCsv = self.outputCsvEdit.text.strip()
        if not outputCsv.lower().endswith(".csv"):
            outputCsv += ".csv"
            self.outputCsvEdit.text = outputCsv

        legSegmentIDs = {}
        planeNodes = {}
        for sideKey, _ in SIDES:
            segId = self.legSegmentCombos[sideKey].currentData
            hasAllPlanes = all(
                self.planeNodes[sideKey][r] is not None for r, _ in PLANE_ROLES
            )
            if segId is not None and hasAllPlanes:
                legSegmentIDs[sideKey] = segId
                planeNodes[sideKey] = dict(self.planeNodes[sideKey])

        statsFlags = {
            "voxel_count": self.cbVoxelCount.checked,
            "volume_mm3":  self.cbVolume.checked,
            "minimum":     self.cbMin.checked,
            "maximum":     self.cbMax.checked,
            "mean":        self.cbMean.checked,
            "stdev":       self.cbStd.checked,
            "median":      self.cbMedian.checked,
        }

        self.runButton.enabled = False
        self.statusLabel.text = "Running …"
        self.statusLabel.setStyleSheet("color:#555;font-style:italic;")
        slicer.app.processEvents()

        try:
            zoneCounts = self.logic.run(
                ctVolumeNode=ctVolume,
                petVolumeNode=petVolume,
                wholeLimbSegmentation=wholeLimbSeg,
                legSegmentIDs=legSegmentIDs,
                planeNodes=planeNodes,
                outputSegmentationNode=outputSeg,
                outputCsvPath=outputCsv,
                statsFlags=statsFlags,
                progressCallback=self.updateStatus,
            )
            summary = ", ".join(
                f"{k}={v}" for k, v in zoneCounts.items()
            )
            self.statusLabel.text = f"✔ Done. Zone voxel counts: {summary}. CSV saved to: {outputCsv}"
            self.statusLabel.setStyleSheet("color:#27ae60;font-weight:bold;")
        except Exception as e:
            slicer.util.errorDisplay(f"PAD Analysis failed:\n{e}")
            self.statusLabel.text = f"✘ Error: {e}"
            self.statusLabel.setStyleSheet("color:#e74c3c;font-weight:bold;")
            logging.exception("PADAnalysis error")
        finally:
            self._onInputChanged()

    def updateStatus(self, message):
        self.statusLabel.text = message
        slicer.app.processEvents()

    def cleanup(self):
        pass


# ---------------------------------------------------------------------------
# Logic
# ---------------------------------------------------------------------------

class PADAnalysisLogic(ScriptedLoadableModuleLogic):
    """
    Core computation logic - no GUI dependencies.
    Can be used from Python scripts / Jupyter without the widget.
    """

    # ------------------------------------------------------------------
    # Plane creation
    # ------------------------------------------------------------------

    def createDefaultPlane(self, ctVolumeNode, sideKey, roleKey, nodeName, color):
        """Create a vtkMRMLMarkupsPlaneNode with a sensible default position
        derived from the CT volume's bounds, ready for the user to drag into
        place."""
        bounds = [0.0] * 6
        ctVolumeNode.GetRASBounds(bounds)
        cx = (bounds[0] + bounds[1]) / 2.0
        cy = (bounds[2] + bounds[3]) / 2.0
        cz = (bounds[4] + bounds[5]) / 2.0
        sizeX = bounds[1] - bounds[0]
        sizeY = bounds[3] - bounds[2]
        sizeZ = bounds[5] - bounds[4]

        # Spread the three default planes along Z so they don't overlap,
        # and offset left/right sides along X for visual separation.
        zOffsetByRole = {"proximal": sizeZ * 0.25, "knee": 0.0, "ankle": -sizeZ * 0.25}
        xOffsetBySide = {"L": -sizeX * 0.15, "R": sizeX * 0.15}
        originX = cx + xOffsetBySide.get(sideKey, 0.0)
        originZ = cz + zOffsetByRole.get(roleKey, 0.0)

        planeNode = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLMarkupsPlaneNode", nodeName
        )
        planeNode.SetPlaneType(slicer.vtkMRMLMarkupsPlaneNode.PlaneTypePointNormal)
        planeNode.CreateDefaultDisplayNodes()
        planeNode.SetOriginWorld([originX, cy, originZ])
        planeNode.SetNormalWorld([0.0, 0.0, 1.0])
        planeSize = max(min(sizeX, sizeY) * 0.7, 1.0)
        planeNode.SetSize(planeSize, planeSize)
        planeNode.SetAttribute("PADAnalysis.side", sideKey)
        planeNode.SetAttribute("PADAnalysis.role", roleKey)

        displayNode = planeNode.GetDisplayNode()
        if displayNode is not None:
            displayNode.SetHandlesInteractive(True)
            displayNode.SetColor(*color)
            displayNode.SetSelectedColor(*color)
            displayNode.SetOpacity(0.4)
            if hasattr(displayNode, "SetTranslationHandleVisibility"):
                displayNode.SetTranslationHandleVisibility(True)
            if hasattr(displayNode, "SetRotationHandleVisibility"):
                displayNode.SetRotationHandleVisibility(True)
            if hasattr(displayNode, "SetScaleHandleVisibility"):
                displayNode.SetScaleHandleVisibility(True)
            if hasattr(displayNode, "SetPropertiesLabelVisibility"):
                displayNode.SetPropertiesLabelVisibility(True)

        return planeNode

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(
        self,
        ctVolumeNode,
        petVolumeNode,
        wholeLimbSegmentation,
        legSegmentIDs,
        planeNodes,
        outputSegmentationNode,
        outputCsvPath,
        statsFlags=None,
        progressCallback=None,
    ):
        """
        Parameters
        ----------
        ctVolumeNode            : vtkMRMLScalarVolumeNode
        petVolumeNode           : vtkMRMLScalarVolumeNode
        wholeLimbSegmentation   : vtkMRMLSegmentationNode
        legSegmentIDs           : dict  side ("L"/"R") -> segmentID
        planeNodes              : dict  side -> role ("proximal"/"knee"/"ankle") -> vtkMRMLMarkupsPlaneNode
                                   (normal of each plane must point toward the foot)
        outputSegmentationNode  : vtkMRMLSegmentationNode  (zone segments written here)
        outputCsvPath           : str
        statsFlags              : dict[str, bool] | None
        progressCallback        : callable(message) | None

        Returns
        -------
        dict  "<Side>_<Zone>" -> voxel count, for the zones that were created.
        """
        if not legSegmentIDs:
            raise RuntimeError(
                "No leg is ready to process. For at least one side, select a "
                "leg segment and create all three cut-off planes."
            )
        if ctVolumeNode.GetParentTransformNode() is not None:
            raise RuntimeError(
                "The CT volume has an active transform. Please harden the "
                "transform on the CT volume (Data module -> right-click -> "
                "Harden Transform) before running PAD Analysis."
            )

        if statsFlags is None:
            statsFlags = {
                "voxel_count": True, "volume_mm3": True,
                "minimum": True, "maximum": True,
                "mean": True, "stdev": True, "median": False,
            }

        def progress(msg):
            logging.info("PADAnalysis: %s", msg)
            if progressCallback:
                progressCallback(msg)

        outputSegmentationNode.CreateDefaultDisplayNodes()
        outputSegmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(ctVolumeNode)

        ijkToRas = vtk.vtkMatrix4x4()
        ctVolumeNode.GetIJKToRASMatrix(ijkToRas)
        m = np.array(
            [[ijkToRas.GetElement(r, c) for c in range(4)] for r in range(4)]
        )

        sideFullName = dict(SIDES)
        zoneVoxelCounts = {}

        for sideKey, segmentID in legSegmentIDs.items():
            progress(f"Classifying {sideFullName[sideKey]} leg voxels …")

            maskArray = slicer.util.arrayFromSegmentBinaryLabelmap(
                wholeLimbSegmentation, segmentID, ctVolumeNode
            )
            if maskArray is None or not np.any(maskArray):
                progress(f"  {sideFullName[sideKey]} leg segment is empty - skipping.")
                continue

            coords = np.argwhere(maskArray)  # N x 3, columns are (k, j, i)
            kk, jj, ii = coords[:, 0], coords[:, 1], coords[:, 2]
            ijk1 = np.column_stack([ii, jj, kk, np.ones(len(ii))])
            ras = (ijk1 @ m.T)[:, :3]

            zoneCount = np.zeros(len(coords), dtype=np.int8)
            for roleKey, _ in PLANE_ROLES:
                planeNode = planeNodes[sideKey][roleKey]
                origin = [0.0, 0.0, 0.0]
                normal = [0.0, 0.0, 0.0]
                planeNode.GetOriginWorld(origin)
                planeNode.GetNormalWorld(normal)
                origin = np.array(origin)
                normal = np.array(normal)
                d = (ras - origin) @ normal
                zoneCount += (d >= 0).astype(np.int8)

            for zoneIdx, zoneName in enumerate(ZONE_NAMES, start=1):
                segmentName = f"{sideFullName[sideKey]}_{zoneName}"
                self._removeSegmentByName(outputSegmentationNode, segmentName)

                sel = coords[zoneCount == zoneIdx]
                voxelCount = len(sel)
                zoneVoxelCounts[segmentName] = voxelCount
                if voxelCount == 0:
                    progress(f"  {segmentName}: 0 voxels (planes may be out of order or mis-oriented).")
                    continue

                outMask = np.zeros_like(maskArray, dtype=np.uint8)
                outMask[sel[:, 0], sel[:, 1], sel[:, 2]] = 1

                newSegmentID = outputSegmentationNode.GetSegmentation().AddEmptySegment(
                    "", segmentName, ZONE_COLORS[zoneName]
                )
                slicer.util.updateSegmentBinaryLabelmapFromArray(
                    outMask, outputSegmentationNode, newSegmentID, ctVolumeNode
                )
                progress(f"  {segmentName}: {voxelCount} voxels")

        if not any(v > 0 for v in zoneVoxelCounts.values()):
            raise RuntimeError(
                "No zone voxels were produced for any leg. Check that the "
                "planes are positioned inside the whole-limb segment and "
                "ordered proximal -> knee -> ankle, and that each plane's "
                "normal points toward the foot (use Flip Normal if needed)."
            )

        progress("Computing PET statistics …")
        self._quantifyAndExport(
            outputSegmentationNode, petVolumeNode, statsFlags, outputCsvPath
        )
        progress("Done.")

        return zoneVoxelCounts

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _removeSegmentByName(self, segmentationNode, name):
        segmentation = segmentationNode.GetSegmentation()
        for i in range(segmentation.GetNumberOfSegments() - 1, -1, -1):
            segId = segmentation.GetNthSegmentID(i)
            if segmentation.GetSegment(segId).GetName() == name:
                segmentation.RemoveSegment(segId)

    def _quantifyAndExport(self, segmentationNode, petVolumeNode, statsFlags, outputCsvPath):
        import SegmentStatistics
        import tempfile

        ssLogic = SegmentStatistics.SegmentStatisticsLogic()
        pNode = ssLogic.getParameterNode()
        pNode.SetParameter("Segmentation", segmentationNode.GetID())
        pNode.SetParameter("ScalarVolume", petVolumeNode.GetID())
        pNode.SetParameter("LabelmapSegmentStatisticsPlugin.enabled", "True")
        pNode.SetParameter("ScalarVolumeSegmentStatisticsPlugin.enabled", "True")
        for meas in ["voxel_count", "volume_mm3", "volume_cm3"]:
            pNode.SetParameter(f"LabelmapSegmentStatisticsPlugin.{meas}.enabled", "True")
        for meas in ["voxel_count", "volume_mm3", "min", "max", "mean", "median", "stdev"]:
            pNode.SetParameter(f"ScalarVolumeSegmentStatisticsPlugin.{meas}.enabled", "True")

        ssLogic.computeStatistics()
        rawCsvPath = tempfile.mktemp(suffix=".csv")
        ssLogic.exportToCSVFile(rawCsvPath)

        with open(rawCsvPath, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames or []
            rows = list(reader)
        try:
            os.remove(rawCsvPath)
        except OSError:
            pass

        STAT_CSV_SUBSTRINGS = {
            "voxel_count": ["voxel count", "voxelcount", "voxel_count", "number of voxels"],
            "volume_mm3":  ["volume mm3", "volume_mm3", "volume [mm3]", "volume[mm3]"],
            "minimum":     ["minimum", "min"],
            "maximum":     ["maximum", "max"],
            "mean":        ["mean"],
            "stdev":       ["standard deviation", "stdev", "std"],
            "median":      ["median"],
        }
        STAT_LABELS = {
            "voxel_count": "VoxelCount", "volume_mm3": "Volume_mm3",
            "minimum": "Min", "maximum": "Max",
            "mean": "Mean", "stdev": "StdDev", "median": "Median",
        }
        enabledStatKeys = [k for k, v in statsFlags.items() if v]

        resolvedCsvCols = {}
        for statKey in enabledStatKeys:
            substrings = STAT_CSV_SUBSTRINGS.get(statKey, [])
            for col in headers:
                colLower = col.strip().lower()
                if any(sub in colLower for sub in substrings):
                    resolvedCsvCols[statKey] = col
                    break

        segNameCol = None
        for col in headers:
            if col.strip().lower() in ("segment", "name", "segment name"):
                segNameCol = col
                break
        if segNameCol is None and headers:
            segNameCol = headers[0]

        outHeaders = ["Segment", "Side", "Zone"] + [
            STAT_LABELS[k] for k in enabledStatKeys if k in resolvedCsvCols
        ]
        outRows = []
        for row in rows:
            segName = row.get(segNameCol, "").strip()
            side, _, zone = segName.partition("_")
            outRow = [segName, side, zone]
            for statKey in enabledStatKeys:
                if statKey in resolvedCsvCols:
                    rawVal = row.get(resolvedCsvCols[statKey], "")
                    try:
                        outRow.append(float(rawVal))
                    except (ValueError, TypeError):
                        outRow.append(rawVal)
            outRows.append(outRow)

        os.makedirs(os.path.dirname(os.path.abspath(outputCsvPath)), exist_ok=True)
        with open(outputCsvPath, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(outHeaders)
            writer.writerows(outRows)


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

class PADAnalysisTest(ScriptedLoadableModuleTest):
    """Self-test using a synthetic volume/segmentation - no real data needed."""

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_ZoneClassification()

    def test_ZoneClassification(self):
        self.delayDisplay("Building synthetic CT/PET volumes and whole-limb segment …")
        import tempfile

        # 10 x 10 x 60 volume, 1 mm isotropic spacing, identity orientation.
        nx, ny, nz = 10, 10, 60
        ctArray = np.zeros((nz, ny, nx), dtype=np.int16)
        ctVolumeNode = slicer.util.addVolumeFromArray(ctArray, name="Test_CT")
        petArray = np.ones((nz, ny, nx), dtype=np.float32)
        petVolumeNode = slicer.util.addVolumeFromArray(petArray, name="Test_PET")

        segNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")
        segNode.CreateDefaultDisplayNodes()
        segNode.SetReferenceImageGeometryParameterFromVolumeNode(ctVolumeNode)
        segId = segNode.GetSegmentation().AddEmptySegment("", "Left Hindlimb", (1, 0, 0))
        wholeLimbMask = np.ones((nz, ny, nx), dtype=np.uint8)
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            wholeLimbMask, segNode, segId, ctVolumeNode
        )

        logic = PADAnalysisLogic()

        def makePlane(name, zPos, color):
            p = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsPlaneNode", name)
            p.SetPlaneType(slicer.vtkMRMLMarkupsPlaneNode.PlaneTypePointNormal)
            p.CreateDefaultDisplayNodes()
            p.SetOriginWorld([4.5, 4.5, zPos])
            p.SetNormalWorld([0.0, 0.0, 1.0])
            return p

        proximalPlane = makePlane("proxPlane", 10.0, PLANE_COLORS["proximal"])
        kneePlane     = makePlane("kneePlane", 30.0, PLANE_COLORS["knee"])
        anklePlane    = makePlane("anklePlane", 50.0, PLANE_COLORS["ankle"])

        outputSeg = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")

        self.delayDisplay("Running PADAnalysisLogic.run() …")
        csvPath = tempfile.mktemp(suffix=".csv")
        result = logic.run(
            ctVolumeNode=ctVolumeNode,
            petVolumeNode=petVolumeNode,
            wholeLimbSegmentation=segNode,
            legSegmentIDs={"L": segId},
            planeNodes={"L": {"proximal": proximalPlane, "knee": kneePlane, "ankle": anklePlane}},
            outputSegmentationNode=outputSeg,
            outputCsvPath=csvPath,
        )

        # z in [0,10) excluded, [10,30) Proximal (2000 vox), [30,50) Distal (2000),
        # [50,60) Foot (1000).
        assert result["Left_Proximal"] == 2000, result
        assert result["Left_Distal"] == 2000, result
        assert result["Left_Foot"] == 1000, result
        assert os.path.isfile(csvPath), "CSV not written"

        with open(csvPath, newline="", encoding="utf-8") as f:
            csvRows = list(csv.DictReader(f))
        assert len(csvRows) == 3, f"expected 3 CSV rows, got {len(csvRows)}"

        os.remove(csvPath)
        self.delayDisplay("Test passed ✓")

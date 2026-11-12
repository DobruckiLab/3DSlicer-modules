"""
PAD Analysis 2 - 3D Slicer Scripted Module
==========================================
Semi-automatic hindlimb segmentation of preclinical PET/CT (supine animal,
legs extended) followed by PET uptake quantification (%ID/g, SUV, and
ischemic / contralateral ratios).

Workflow
--------
  1. Body mask   - CT threshold + morphological opening + largest connected
                   component (removes bed / noise), stored as a "BodyMask"
                   segment that can be inspected / edited in Segment Editor.
  2. Landmarks   - 3 clicks per leg on the CT: hip joint, knee joint, ankle.
  3. Zones       - the body mask is split per leg into
                     <Side>_Proximal  (hip level -> knee, inside the body
                                       outline, lateral to a medial cut)
                     <Side>_Distal    (knee -> ankle)
                     <Side>_Foot      (beyond the ankle)
                   using landmark-derived cut planes, a capture radius
                   around the bone axis and connected-component checks
                   (keeps the tail / other leg out).
  4. Quantify    - PET statistics per zone, converted to %ID/g and SUV,
                   ischemic / control ratios per zone; written to CSV and
                   shown as Slicer tables.

Author: Wawrzyniec Dobrucki (EMIL, UIUC)
Slicer version: 5.x+
"""

import os
import csv
import math
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
LANDMARKS = ["Hip", "Knee", "Ankle"]          # click order
ZONE_NAMES = ["Proximal", "Distal", "Foot"]
BODY_SEGMENT_NAME = "BodyMask"
WORKING_VOLUME_REF = "PADAnalysis2WorkingVolume"

# Colours roughly matching the reference figure (warm = left, cool = right)
ZONE_COLORS = {
    ("L", "Proximal"): (1.00, 0.60, 0.10),
    ("L", "Distal"):   (0.85, 0.90, 0.10),
    ("L", "Foot"):     (0.20, 0.75, 0.20),
    ("R", "Proximal"): (0.30, 0.85, 0.80),
    ("R", "Distal"):   (0.20, 0.65, 1.00),
    ("R", "Foot"):     (0.10, 0.20, 0.90),
}
LANDMARK_COLORS = {"L": (1.0, 0.5, 0.0), "R": (0.0, 0.6, 1.0)}

# Unit conversions
PET_UNITS_TO_BQML = {
    "Bq/mL": 1.0,
    "kBq/mL": 1.0e3,
    "MBq/mL": 1.0e6,
    "nCi/mL": 37.0,
    "µCi/mL": 3.7e4,
}
DOSE_UNITS_TO_BQ = {
    "MBq": 1.0e6,
    "kBq": 1.0e3,
    "mCi": 3.7e7,
    "µCi": 3.7e4,
}
# Physical half-lives in minutes
ISOTOPE_HALF_LIFE_MIN = {
    "F-18": 109.77,
    "C-11": 20.364,
    "N-13": 9.965,
    "O-15": 2.037,
    "Ga-68": 67.71,
    "Cu-64": 762.06,
    "Zr-89": 4704.6,
    "I-124": 6013.4,
    "Custom": None,
}

STATS_HEADER = [
    "SubjectID", "Segment", "Side", "Zone", "Condition",
    "Volume_mm3", "PET_Voxels",
    "Mean", "StdDev", "Median", "Min", "Max", "PET_Units",
    "%ID/g_Mean", "%ID/g_StdDev", "%ID/g_Max",
    "SUV_Mean", "SUV_StdDev", "SUV_Max",
    "NetDose_MBq", "BodyWeight_g",
]
RATIO_HEADER = [
    "SubjectID", "Zone", "IschemicSide", "ControlSide",
    "Ratio_Mean", "Ratio_Max",
]


# ---------------------------------------------------------------------------
# Module descriptor
# ---------------------------------------------------------------------------

class PADAnalysis2(ScriptedLoadableModule):
    """Module metadata shown in the Slicer module finder."""

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "PAD Analysis 2"
        self.parent.categories = ["Quantification"]
        self.parent.dependencies = ["Markups", "Segmentations"]
        self.parent.contributors = [
            "Wawrzyniec Dobrucki <br> Experimental Molecular Imaging Laboratory "
            "<br> University of Illinois at Urbana-Champaign"
        ]
        self.parent.helpText = (
            "Semi-automatic hindlimb PET/CT analysis for PAD / hindlimb "
            "ischemia models (supine animal, legs extended)."
            "<ol>"
            "<li><b>Body mask</b>: thresholds the CT and keeps the animal "
            "(bed removed). Check / edit the <i>BodyMask</i> segment in "
            "Segment Editor if needed.</li>"
            "<li><b>Landmarks</b>: for each leg click the <i>hip joint</i>, "
            "<i>knee joint</i> and <i>ankle joint</i> (in this order) on "
            "the CT. Points can be dragged afterwards.</li>"
            "<li><b>Zones</b>: generates <i>Left/Right_Proximal</i>, "
            "<i>_Distal</i> and <i>_Foot</i> segments. They can be refined "
            "in Segment Editor before quantification.</li>"
            "<li><b>Quantify</b>: PET mean/SD/median/min/max per zone, "
            "converted to %ID/g and SUV, plus ischemic / control ratios. "
            "Results go to CSV and to Slicer tables.</li>"
            "</ol>"
            "%ID/g and SUV assume a tissue density of 1 g/mL."
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

class PADAnalysis2Widget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """User interface for the module."""

    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self.landmarkNodes = {s: None for s, _ in SIDES}
        self.landmarkStatusLabels = {}

    # ------------------------------------------------------------------
    # setup
    # ------------------------------------------------------------------

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = PADAnalysis2Logic()

        # ── Inputs ────────────────────────────────────────────────────
        inputBox = ctk.ctkCollapsibleButton()
        inputBox.text = "Inputs"
        self.layout.addWidget(inputBox)
        inputForm = qt.QFormLayout(inputBox)

        self.ctSelector = self._makeVolumeSelector("CT volume (calibrated in HU).")
        inputForm.addRow("CT Volume:", self.ctSelector)
        self.petSelector = self._makeVolumeSelector(
            "PET volume (already aligned with the CT)."
        )
        inputForm.addRow("PET Volume:", self.petSelector)

        self.segSelector = slicer.qMRMLNodeComboBox()
        self.segSelector.nodeTypes = ["vtkMRMLSegmentationNode"]
        self.segSelector.selectNodeUponCreation = True
        self.segSelector.addEnabled = True
        self.segSelector.removeEnabled = True
        self.segSelector.renameEnabled = True
        self.segSelector.noneEnabled = False
        self.segSelector.baseName = "PAD2_Segmentation"
        self.segSelector.setMRMLScene(slicer.mrmlScene)
        self.segSelector.setToolTip(
            "Segmentation that receives the BodyMask and the six zone segments."
        )
        inputForm.addRow("Segmentation:", self.segSelector)

        self.subjectEdit = qt.QLineEdit()
        self.subjectEdit.setPlaceholderText("e.g. Mouse01_Day7")
        inputForm.addRow("Subject ID:", self.subjectEdit)

        # ── Step 1: body mask ─────────────────────────────────────────
        bodyBox = ctk.ctkCollapsibleButton()
        bodyBox.text = "Step 1 - Body Mask (CT)"
        self.layout.addWidget(bodyBox)
        bodyForm = qt.QFormLayout(bodyBox)

        self.thresholdSpin = self._makeDoubleSpin(-1000, 3000, -300, 0, " HU", 10)
        self.thresholdSpin.setToolTip(
            "Voxels above this CT value are considered tissue."
        )
        bodyForm.addRow("Tissue threshold:", self.thresholdSpin)

        self.openingSpin = self._makeDoubleSpin(0, 5, 0.5, 2, " mm", 0.1)
        self.openingSpin.setToolTip(
            "Morphological opening radius used to detach the bed / thin "
            "structures from the animal (0 = off)."
        )
        bodyForm.addRow("Opening radius:", self.openingSpin)

        self.workSpacingSpin = self._makeDoubleSpin(0, 2, 0.25, 3, " mm", 0.05)
        self.workSpacingSpin.setToolTip(
            "Isotropic working resolution for segmentation. Use a value "
            "coarser than the CT to speed things up (PET resolution is much "
            "coarser anyway). 0 = native CT grid."
        )
        bodyForm.addRow("Working resolution:", self.workSpacingSpin)

        self.bodyButton = self._makeActionButton("Create Body Mask")
        bodyForm.addRow(self.bodyButton)

        # ── Step 2: landmarks ─────────────────────────────────────────
        lmBox = ctk.ctkCollapsibleButton()
        lmBox.text = "Step 2 - Landmarks (3 clicks per leg)"
        self.layout.addWidget(lmBox)
        lmLayout = qt.QVBoxLayout(lmBox)
        note = qt.QLabel(
            "Click <b>Place</b>, then click on the CT in this order: "
            "<b>hip joint → knee joint → ankle joint</b>. "
            "Points can be dragged afterwards."
        )
        note.setWordWrap(True)
        lmLayout.addWidget(note)
        for sideKey, sideLabel in SIDES:
            row = qt.QHBoxLayout()
            row.addWidget(qt.QLabel(f"{sideLabel} leg:"))
            placeBtn = qt.QPushButton("Place")
            resetBtn = qt.QPushButton("Reset")
            jumpBtn = qt.QPushButton("Jump")
            status = qt.QLabel("")
            status.setMinimumWidth(200)
            placeBtn.clicked.connect(lambda checked=False, s=sideKey: self.onPlaceLandmarks(s))
            resetBtn.clicked.connect(lambda checked=False, s=sideKey: self.onResetLandmarks(s))
            jumpBtn.clicked.connect(lambda checked=False, s=sideKey: self.onJumpLandmarks(s))
            row.addWidget(placeBtn)
            row.addWidget(resetBtn)
            row.addWidget(jumpBtn)
            row.addWidget(status, 1)
            self.landmarkStatusLabels[sideKey] = status
            lmLayout.addLayout(row)

        # ── Step 3: zones ─────────────────────────────────────────────
        zoneBox = ctk.ctkCollapsibleButton()
        zoneBox.text = "Step 3 - Generate Zones"
        self.layout.addWidget(zoneBox)
        zoneForm = qt.QFormLayout(zoneBox)

        self.hipMarginSpin = self._makeDoubleSpin(-10, 20, 2.0, 1, " mm", 0.5)
        self.hipMarginSpin.setToolTip(
            "Upper boundary of the Proximal zone, measured from the hip "
            "joint toward the head along the body axis."
        )
        zoneForm.addRow("Proximal top above hip:", self.hipMarginSpin)

        self.medialSpin = self._makeDoubleSpin(0, 100, 50, 0, " %", 5)
        self.medialSpin.setToolTip(
            "Medial boundary of the Proximal zone, as % of the distance from "
            "the body midline to the hip joint (0 % = midline, 100 % = at "
            "the hip joint). Increase to keep bladder / genital signal out."
        )
        zoneForm.addRow("Medial boundary:", self.medialSpin)

        self.captureSpin = self._makeDoubleSpin(0, 50, 0, 1, " mm", 0.5)
        self.captureSpin.setSpecialValueText("Auto")
        self.captureSpin.setToolTip(
            "Max. distance from the knee–ankle axis for Distal-zone "
            "voxels. Auto = 35 % of the knee–ankle distance."
        )
        zoneForm.addRow("Limb capture radius:", self.captureSpin)

        self.footReachSpin = self._makeDoubleSpin(0, 100, 0, 1, " mm", 0.5)
        self.footReachSpin.setSpecialValueText("Auto")
        self.footReachSpin.setToolTip(
            "Max. distance from the ankle landmark for Foot-zone voxels. "
            "Auto = knee–ankle distance."
        )
        zoneForm.addRow("Foot reach:", self.footReachSpin)

        self.zonesButton = self._makeActionButton("Generate Zones")
        zoneForm.addRow(self.zonesButton)

        # ── Step 4: quantification ────────────────────────────────────
        qBox = ctk.ctkCollapsibleButton()
        qBox.text = "Step 4 - PET Quantification"
        self.layout.addWidget(qBox)
        qForm = qt.QFormLayout(qBox)

        self.petUnitsCombo = qt.QComboBox()
        self.petUnitsCombo.addItems(list(PET_UNITS_TO_BQML.keys()))
        self.petUnitsCombo.setToolTip("Units of the PET image voxel values.")
        qForm.addRow("PET image units:", self.petUnitsCombo)

        doseRow = qt.QHBoxLayout()
        self.doseSpin = self._makeDoubleSpin(0, 1e6, 0, 3, "", 0.1)
        self.residualSpin = self._makeDoubleSpin(0, 1e6, 0, 3, "", 0.1)
        self.doseUnitsCombo = qt.QComboBox()
        self.doseUnitsCombo.addItems(list(DOSE_UNITS_TO_BQ.keys()))
        doseRow.addWidget(qt.QLabel("Syringe:"))
        doseRow.addWidget(self.doseSpin)
        doseRow.addWidget(qt.QLabel("Residual:"))
        doseRow.addWidget(self.residualSpin)
        doseRow.addWidget(self.doseUnitsCombo)
        self.doseSpin.setToolTip("Activity in the syringe before injection.")
        self.residualSpin.setToolTip("Activity left in the syringe after injection.")
        qForm.addRow("Injected activity:", doseRow)

        self.weightSpin = self._makeDoubleSpin(0, 5000, 0, 2, " g", 0.5)
        qForm.addRow("Body weight:", self.weightSpin)

        decayRow = qt.QHBoxLayout()
        self.decayCheck = qt.QCheckBox("Decay-correct dose to scan start")
        self.decayCheck.setToolTip(
            "Enable if the PET image is decay-corrected to scan start (the "
            "usual case). The net dose is then decayed from injection time "
            "to scan start. Disable if the PET is decay-corrected to "
            "injection time."
        )
        self.decayCheck.checked = True
        self.isotopeCombo = qt.QComboBox()
        self.isotopeCombo.addItems(list(ISOTOPE_HALF_LIFE_MIN.keys()))
        self.halfLifeSpin = self._makeDoubleSpin(0.01, 1e6, 109.77, 2, " min", 1)
        self.intervalSpin = self._makeDoubleSpin(0, 1e6, 60, 1, " min", 1)
        self.intervalSpin.setToolTip("Time from injection to PET scan start.")
        decayRow.addWidget(self.decayCheck)
        decayRow.addWidget(self.isotopeCombo)
        decayRow.addWidget(self.halfLifeSpin)
        qForm.addRow("Decay:", decayRow)
        qForm.addRow("Injection → scan:", self.intervalSpin)

        self.ischemicCombo = qt.QComboBox()
        self.ischemicCombo.addItem("Left", "L")
        self.ischemicCombo.addItem("Right", "R")
        self.ischemicCombo.setToolTip(
            "Ischemic (ligated) leg. Ratios are ischemic / contralateral."
        )
        qForm.addRow("Ischemic leg:", self.ischemicCombo)

        self.erodeSpin = self._makeDoubleSpin(0, 5, 0, 2, " mm", 0.1)
        self.erodeSpin.setToolTip(
            "Shrink each zone by this margin before measuring PET, to reduce "
            "partial-volume spill-out at the skin (segments are not changed)."
        )
        qForm.addRow("Exclude skin margin:", self.erodeSpin)

        fileRow = qt.QHBoxLayout()
        self.csvEdit = qt.QLineEdit()
        self.csvEdit.setPlaceholderText("Output CSV path …")
        browseBtn = qt.QPushButton("Browse…")
        browseBtn.clicked.connect(self.onBrowseCsv)
        fileRow.addWidget(self.csvEdit)
        fileRow.addWidget(browseBtn)
        qForm.addRow("Output CSV:", fileRow)
        self.appendCheck = qt.QCheckBox(
            "Append to existing file (one file for many animals)"
        )
        self.appendCheck.checked = True
        qForm.addRow("", self.appendCheck)

        self.quantButton = self._makeActionButton("Quantify PET")
        qForm.addRow(self.quantButton)

        # ── Status ────────────────────────────────────────────────────
        self.statusLabel = qt.QLabel("")
        self.statusLabel.setWordWrap(True)
        self.statusLabel.setStyleSheet("color:#555;font-style:italic;")
        self.layout.addWidget(self.statusLabel)
        self.layout.addStretch()

        # Connections
        self.bodyButton.clicked.connect(self.onCreateBodyMask)
        self.zonesButton.clicked.connect(self.onGenerateZones)
        self.quantButton.clicked.connect(self.onQuantify)
        self.isotopeCombo.currentIndexChanged.connect(self.onIsotopeChanged)
        self.decayCheck.toggled.connect(self.onDecayToggled)

        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneClosed)
        self._refreshLandmarkNodes()

    # ------------------------------------------------------------------
    # UI helpers
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

    def _makeDoubleSpin(self, lo, hi, value, decimals, suffix, step):
        s = qt.QDoubleSpinBox()
        s.setRange(lo, hi)
        s.setDecimals(decimals)
        s.setSingleStep(step)
        s.setValue(value)
        if suffix:
            s.setSuffix(suffix)
        return s

    def _makeActionButton(self, text):
        b = qt.QPushButton(text)
        b.setStyleSheet(
            "QPushButton { background-color: #2ecc71; color: white; "
            "font-weight: bold; padding: 6px; border-radius: 4px; }"
            "QPushButton:hover { background-color: #27ae60; }"
        )
        return b

    def _setStatus(self, message, kind="info"):
        style = {
            "info": "color:#555;font-style:italic;",
            "ok": "color:#27ae60;font-weight:bold;",
            "error": "color:#e74c3c;font-weight:bold;",
        }[kind]
        self.statusLabel.setStyleSheet(style)
        self.statusLabel.text = message
        slicer.app.processEvents()

    def _run(self, title, func):
        """Run func() with a wait cursor and uniform error reporting."""
        self._setStatus(f"{title} …")
        try:
            with slicer.util.tryWithErrorDisplay(f"{title} failed.", waitCursor=True):
                message = func()
                self._setStatus(f"✔ {message}", "ok")
                return
        except Exception as e:  # tryWithErrorDisplay re-raises
            logging.exception("PADAnalysis2: %s failed", title)
            self._setStatus(f"✘ {title} failed: {e}", "error")

    # ------------------------------------------------------------------
    # Landmarks
    # ------------------------------------------------------------------

    def _refreshLandmarkNodes(self):
        for sideKey, _ in SIDES:
            node = slicer.mrmlScene.GetFirstNodeByName(self.logic.landmarkNodeName(sideKey))
            self._setLandmarkNode(sideKey, node)

    def _setLandmarkNode(self, sideKey, node):
        old = self.landmarkNodes[sideKey]
        if old is not None and old is not node:
            self.removeObserver(old, vtk.vtkCommand.ModifiedEvent, self._onLandmarksModified)
            for ev in self._markupsEvents():
                self.removeObserver(old, ev, self._onLandmarksModified)
        self.landmarkNodes[sideKey] = node
        if node is not None and old is not node:
            for ev in self._markupsEvents():
                self.addObserver(node, ev, self._onLandmarksModified)
        self._updateLandmarkStatus()

    def _markupsEvents(self):
        cls = slicer.vtkMRMLMarkupsNode
        return [cls.PointPositionDefinedEvent, cls.PointPositionUndefinedEvent,
                cls.PointModifiedEvent, cls.PointRemovedEvent]

    def _onLandmarksModified(self, caller=None, event=None):
        self._updateLandmarkStatus()
        # leave place mode once the active leg has all 3 points
        if caller is not None and self.logic.landmarksComplete(caller):
            selectionNode = slicer.app.applicationLogic().GetSelectionNode()
            if selectionNode.GetActivePlaceNodeID() == caller.GetID():
                interactionNode = slicer.app.applicationLogic().GetInteractionNode()
                if interactionNode.GetCurrentInteractionMode() == interactionNode.Place:
                    interactionNode.SwitchToViewTransformMode()

    def _updateLandmarkStatus(self):
        for sideKey, _ in SIDES:
            node = self.landmarkNodes[sideKey]
            label = self.landmarkStatusLabels.get(sideKey)
            if label is None:
                continue
            placed = self.logic.placedLandmarks(node) if node else []
            n = len(placed)
            if n == 3:
                label.text = "✔ hip, knee, ankle placed"
                label.setStyleSheet("color:#27ae60;")
            else:
                nextName = LANDMARKS[n].lower() if n < 3 else ""
                label.text = f"{n}/3 placed — next: {nextName}"
                label.setStyleSheet("color:#e67e22;")

    def onPlaceLandmarks(self, sideKey):
        node = self.landmarkNodes[sideKey]
        if node is None or slicer.mrmlScene.GetNodeByID(node.GetID()) is None:
            node = self.logic.createLandmarkNode(sideKey)
            self._setLandmarkNode(sideKey, node)
        if self.logic.landmarksComplete(node):
            self._setStatus(
                f"{dict(SIDES)[sideKey]} landmarks already placed - drag them "
                "to adjust, or Reset to place again."
            )
            return
        selectionNode = slicer.app.applicationLogic().GetSelectionNode()
        selectionNode.SetReferenceActivePlaceNodeClassName("vtkMRMLMarkupsFiducialNode")
        selectionNode.SetActivePlaceNodeID(node.GetID())
        interactionNode = slicer.app.applicationLogic().GetInteractionNode()
        interactionNode.SetPlaceModePersistence(1)
        interactionNode.SetCurrentInteractionMode(interactionNode.Place)
        self._setStatus(f"Placing {dict(SIDES)[sideKey].lower()} leg landmarks: hip → knee → ankle")

    def onResetLandmarks(self, sideKey):
        node = self.landmarkNodes[sideKey]
        if node is None:
            return
        for i in range(node.GetNumberOfControlPoints()):
            node.UnsetNthControlPointPosition(i)
        self._updateLandmarkStatus()

    def onJumpLandmarks(self, sideKey):
        node = self.landmarkNodes[sideKey]
        if node is None:
            return
        placed = self.logic.placedLandmarks(node)
        if placed:
            slicer.modules.markups.logic().JumpSlicesToNthPointInMarkup(
                node.GetID(), placed[0], True
            )

    def onSceneClosed(self, caller=None, event=None):
        for sideKey, _ in SIDES:
            self._setLandmarkNode(sideKey, None)

    # ------------------------------------------------------------------
    # Actions
    # ------------------------------------------------------------------

    def onCreateBodyMask(self):
        ct = self.ctSelector.currentNode()
        seg = self.segSelector.currentNode()
        if ct is None or seg is None:
            slicer.util.errorDisplay("Select a CT volume and a segmentation first.")
            return

        def work():
            n = self.logic.createBodyMask(
                ct, seg,
                thresholdHU=self.thresholdSpin.value,
                openingRadiusMm=self.openingSpin.value,
                workingSpacingMm=self.workSpacingSpin.value,
                progressCallback=lambda m: self._setStatus(m),
            )
            return (f"Body mask created ({n:,} voxels). Check it in Segment "
                    "Editor if needed, then place the landmarks.")
        self._run("Body mask", work)

    def onGenerateZones(self):
        seg = self.segSelector.currentNode()
        if seg is None:
            slicer.util.errorDisplay("Select the segmentation first.")
            return
        landmarks = {}
        for sideKey, _ in SIDES:
            node = self.landmarkNodes[sideKey]
            if node is not None and self.logic.landmarksComplete(node):
                landmarks[sideKey] = self.logic.landmarkPositions(node)

        def work():
            counts = self.logic.generateZones(
                seg, landmarks,
                hipMarginMm=self.hipMarginSpin.value,
                medialFraction=self.medialSpin.value / 100.0,
                captureRadiusMm=self.captureSpin.value,
                footReachMm=self.footReachSpin.value,
                progressCallback=lambda m: self._setStatus(m),
            )
            return "Zones generated: " + ", ".join(f"{k}={v:,}" for k, v in counts.items())
        self._run("Zone generation", work)

    def onQuantify(self):
        seg = self.segSelector.currentNode()
        pet = self.petSelector.currentNode()
        csvPath = self.csvEdit.text.strip()
        if seg is None or pet is None:
            slicer.util.errorDisplay("Select the PET volume and the segmentation.")
            return
        if not csvPath:
            slicer.util.errorDisplay("Choose an output CSV file.")
            return
        if not csvPath.lower().endswith(".csv"):
            csvPath += ".csv"
            self.csvEdit.text = csvPath

        def work():
            result = self.logic.quantify(
                seg, pet,
                petUnits=self.petUnitsCombo.currentText,
                syringeActivity=self.doseSpin.value,
                residualActivity=self.residualSpin.value,
                doseUnits=self.doseUnitsCombo.currentText,
                bodyWeightG=self.weightSpin.value,
                decayCorrect=self.decayCheck.checked,
                halfLifeMin=self.halfLifeSpin.value,
                injectionToScanMin=self.intervalSpin.value,
                ischemicSide=self.ischemicCombo.currentData,
                skinMarginMm=self.erodeSpin.value,
                subjectID=self.subjectEdit.text.strip(),
                outputCsvPath=csvPath,
                appendToCsv=self.appendCheck.checked,
                progressCallback=lambda m: self._setStatus(m),
            )
            self.logic.showTable(result["statsTable"])
            return (f"Quantified {len(result['stats'])} zones; "
                    f"CSV: {csvPath} (+ _ratios.csv)")
        self._run("Quantification", work)

    def onBrowseCsv(self):
        path = qt.QFileDialog.getSaveFileName(
            self.parent, "Output CSV", "", "CSV (*.csv)"
        )
        if path:
            if not path.lower().endswith(".csv"):
                path += ".csv"
            self.csvEdit.text = path

    def onIsotopeChanged(self, index):
        hl = ISOTOPE_HALF_LIFE_MIN.get(self.isotopeCombo.currentText)
        if hl is not None:
            self.halfLifeSpin.setValue(hl)
        self.halfLifeSpin.setEnabled(self.decayCheck.checked)

    def onDecayToggled(self, checked):
        for w in (self.isotopeCombo, self.halfLifeSpin, self.intervalSpin):
            w.setEnabled(checked)

    def cleanup(self):
        self.removeObservers()


# ---------------------------------------------------------------------------
# Logic
# ---------------------------------------------------------------------------

def _matrixToNumpy(vtkMat):
    return np.array([[vtkMat.GetElement(r, c) for c in range(4)] for r in range(4)])


def _ijkToRas(volumeNode):
    m = vtk.vtkMatrix4x4()
    volumeNode.GetIJKToRASMatrix(m)
    return _matrixToNumpy(m)


def _unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-9:
        raise ValueError("Degenerate direction (landmarks coincide?)")
    return v / n


class _Grid:
    """Geometry helper for a (k, j, i) numpy array on a volume grid,
    optionally cropped. All index arguments/returns are local (cropped)."""

    def __init__(self, M, shape, offsetIJK=(0, 0, 0)):
        M = np.array(M, dtype=float)
        M[:3, 3] = M[:3, 3] + M[:3, :3] @ np.asarray(offsetIJK, dtype=float)
        self.M = M
        self.Minv = np.linalg.inv(M)
        self.shape = tuple(shape)                 # (nk, nj, ni)
        self.spacing = np.linalg.norm(M[:3, :3], axis=0)   # (si, sj, sk)

    def ras(self, kji):
        """N x 3 (k, j, i) indices -> N x 3 RAS."""
        ijk = np.asarray(kji, dtype=float)[:, ::-1]
        return ijk @ self.M[:3, :3].T + self.M[:3, 3]

    def ijk(self, ras):
        ras = np.atleast_2d(ras)
        return ras @ self.Minv[:3, :3].T + self.Minv[:3, 3]

    def planeField(self, origin, normal, box=None):
        """Signed distance (mm) to plane for every voxel of box
        (((k0,k1),(j0,j1),(i0,i1)); default whole grid)."""
        if box is None:
            box = tuple((0, n) for n in self.shape)
        (k0, k1), (j0, j1), (i0, i1) = box
        normal = _unit(normal)
        a = (normal @ self.M[:3, :3]).astype(np.float32)
        c = np.float32(normal @ (self.M[:3, 3] - np.asarray(origin, dtype=float)))
        i = np.arange(i0, i1, dtype=np.float32)[None, None, :] * a[0]
        j = np.arange(j0, j1, dtype=np.float32)[None, :, None] * a[1]
        k = np.arange(k0, k1, dtype=np.float32)[:, None, None] * a[2]
        return i + j + k + c

    def boxAroundRas(self, lo, hi):
        """Index box (clipped to grid) covering an axis-aligned RAS box."""
        corners = np.array([[x, y, z] for x in (lo[0], hi[0])
                            for y in (lo[1], hi[1]) for z in (lo[2], hi[2])])
        ijk = self.ijk(corners)
        mn = np.floor(ijk.min(axis=0)).astype(int) - 1
        mx = np.ceil(ijk.max(axis=0)).astype(int) + 2
        box = []
        for axis, n in zip((2, 1, 0), self.shape):     # k, j, i
            box.append((int(np.clip(mn[axis], 0, n)), int(np.clip(mx[axis], 0, n))))
        return tuple(box)


def _keepComponentNearest(mask, grid, seedRas):
    """Keep only the face-connected component of mask closest to seedRas."""
    import SimpleITK as sitk
    if not mask.any():
        return mask
    img = sitk.GetImageFromArray(mask.astype(np.uint8))
    labels = sitk.GetArrayFromImage(sitk.ConnectedComponent(img, False))
    idx = np.argwhere(mask)
    seedKji = grid.ijk(seedRas)[0][::-1]
    spacingKji = grid.spacing[::-1]
    d2 = (((idx - seedKji) * spacingKji) ** 2).sum(axis=1)
    nearest = idx[np.argmin(d2)]
    return labels == labels[tuple(nearest)]


def _distanceToSegment(P, A, B):
    AB = B - A
    denom = float(AB @ AB)
    t = np.clip(((P - A) @ AB) / denom, 0.0, 1.0) if denom > 0 else np.zeros(len(P))
    return np.linalg.norm(P - (A + t[:, None] * AB), axis=1)


class PADAnalysis2Logic(ScriptedLoadableModuleLogic):
    """Computation logic - usable from the Python console without the GUI."""

    # ------------------------------------------------------------------
    # Landmarks
    # ------------------------------------------------------------------

    @staticmethod
    def landmarkNodeName(sideKey):
        return f"PAD2_{dict(SIDES)[sideKey]}_Landmarks"

    def createLandmarkNode(self, sideKey):
        """Point list with 3 predefined (unplaced) points: hip, knee, ankle."""
        node = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLMarkupsFiducialNode", self.landmarkNodeName(sideKey)
        )
        node.CreateDefaultDisplayNodes()
        for name in LANDMARKS:
            idx = node.AddControlPoint(vtk.vtkVector3d(0, 0, 0), f"{sideKey}-{name}")
            node.UnsetNthControlPointPosition(idx)
        if hasattr(node, "SetMaximumNumberOfControlPoints"):
            node.SetMaximumNumberOfControlPoints(len(LANDMARKS))
        if hasattr(node, "SetFixedNumberOfControlPoints"):
            node.SetFixedNumberOfControlPoints(True)
        dn = node.GetDisplayNode()
        if dn is not None:
            dn.SetSelectedColor(*LANDMARK_COLORS[sideKey])
            dn.SetColor(*LANDMARK_COLORS[sideKey])
        node.SetAttribute("PADAnalysis2.side", sideKey)
        return node

    @staticmethod
    def placedLandmarks(node):
        defined = slicer.vtkMRMLMarkupsNode.PositionDefined
        return [i for i in range(min(node.GetNumberOfControlPoints(), len(LANDMARKS)))
                if node.GetNthControlPointPositionStatus(i) == defined]

    def landmarksComplete(self, node):
        return node is not None and len(self.placedLandmarks(node)) == len(LANDMARKS)

    @staticmethod
    def landmarkPositions(node):
        """-> dict {"Hip": ras, "Knee": ras, "Ankle": ras} (world coordinates)."""
        out = {}
        for i, name in enumerate(LANDMARKS):
            p = [0.0, 0.0, 0.0]
            node.GetNthControlPointPositionWorld(i, p)
            out[name] = np.array(p)
        return out

    # ------------------------------------------------------------------
    # Step 1: body mask
    # ------------------------------------------------------------------

    def createBodyMask(self, ctVolumeNode, segmentationNode, thresholdHU=-300.0,
                       openingRadiusMm=0.5, workingSpacingMm=0.25,
                       progressCallback=None):
        import SimpleITK as sitk
        import sitkUtils

        progress = self._progress(progressCallback)
        self._checkNoTransform(ctVolumeNode, "CT")

        # -- working grid --------------------------------------------------
        ctSpacing = np.array(ctVolumeNode.GetSpacing())
        workNode = ctVolumeNode
        if workingSpacingMm > 0 and np.any(ctSpacing < workingSpacingMm * 0.95):
            progress(f"Resampling CT to {workingSpacingMm:g} mm working grid …")
            img = sitkUtils.PullVolumeFromSlicer(ctVolumeNode)
            img = sitk.Cast(img, sitk.sitkFloat32)
            newSpacing = [float(workingSpacingMm)] * 3
            newSize = [int(math.ceil(sz * sp / workingSpacingMm))
                       for sz, sp in zip(img.GetSize(), img.GetSpacing())]
            res = sitk.Resample(img, newSize, sitk.Transform(), sitk.sitkLinear,
                                img.GetOrigin(), newSpacing, img.GetDirection(),
                                -1024.0, sitk.sitkFloat32)
            name = f"PAD2_WorkingCT_{ctVolumeNode.GetName()}"
            old = slicer.mrmlScene.GetFirstNodeByName(name)
            if old is not None:
                slicer.mrmlScene.RemoveNode(old)
            workNode = sitkUtils.PushVolumeToSlicer(res, name=name)
            workNode.SetHideFromEditors(True)

        # -- threshold + cleanup ------------------------------------------
        progress("Thresholding CT and removing bed / noise …")
        ctArr = slicer.util.arrayFromVolume(workNode)
        mask = (ctArr > thresholdHU).astype(np.uint8)
        if not mask.any():
            raise RuntimeError("No voxels above the threshold - check CT units (HU).")
        minSpacing = float(np.min(workNode.GetSpacing()))
        img = sitk.GetImageFromArray(mask)
        r = int(round(openingRadiusMm / minSpacing))
        if r > 0:
            img = sitk.BinaryMorphologicalOpening(img, [r] * 3)
        cc = sitk.RelabelComponent(sitk.ConnectedComponent(img), sortByObjectSize=True)
        img = sitk.Cast(cc == 1, sitk.sitkUInt8)
        img = sitk.BinaryFillhole(img)
        body = sitk.GetArrayFromImage(img).astype(np.uint8)

        # -- write segment -------------------------------------------------
        segmentationNode.CreateDefaultDisplayNodes()
        segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(workNode)
        segmentationNode.SetNodeReferenceID(WORKING_VOLUME_REF, workNode.GetID())
        self._removeSegmentByName(segmentationNode, BODY_SEGMENT_NAME)
        segId = segmentationNode.GetSegmentation().AddEmptySegment(
            "", BODY_SEGMENT_NAME, (0.85, 0.20, 0.20)
        )
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            body, segmentationNode, segId, workNode
        )
        progress("Body mask done.")
        return int(body.sum())

    # ------------------------------------------------------------------
    # Step 3: zones
    # ------------------------------------------------------------------

    def generateZones(self, segmentationNode, landmarks, hipMarginMm=2.0,
                      medialFraction=0.5, captureRadiusMm=0.0, footReachMm=0.0,
                      progressCallback=None):
        """
        landmarks : dict side ("L"/"R") -> {"Hip": ras, "Knee": ras, "Ankle": ras}
        Returns dict "<Side>_<Zone>" -> voxel count.
        """
        progress = self._progress(progressCallback)
        if not landmarks:
            raise RuntimeError("Place all three landmarks (hip, knee, ankle) "
                               "for at least one leg.")
        workNode = self.workingVolume(segmentationNode)
        bodyId = self._segmentIdByName(segmentationNode, BODY_SEGMENT_NAME)
        if workNode is None or bodyId is None:
            raise RuntimeError("No BodyMask found - run Step 1 first.")

        bodyFull = slicer.util.arrayFromSegmentBinaryLabelmap(
            segmentationNode, bodyId, workNode
        ).astype(bool)
        if not bodyFull.any():
            raise RuntimeError("BodyMask is empty.")

        # Crop to the body bounding box
        nz = np.argwhere(bodyFull)
        lo, hi = nz.min(axis=0), nz.max(axis=0) + 1
        crop = tuple(slice(a, b) for a, b in zip(lo, hi))
        body = bodyFull[crop]
        grid = _Grid(_ijkToRas(workNode), body.shape, offsetIJK=lo[::-1])

        lm = {s: {k: np.asarray(v, dtype=float) for k, v in d.items()}
              for s, d in landmarks.items()}

        # Body long axis (PCA of body voxels), pointing toward the feet
        progress("Estimating body axis …")
        idx = np.argwhere(body)
        if len(idx) > 200000:
            idx = idx[np.random.default_rng(0).choice(len(idx), 200000, replace=False)]
        pts = grid.ras(idx)
        evals, evecs = np.linalg.eigh(np.cov((pts - pts.mean(axis=0)).T))
        u = evecs[:, np.argmax(evals)]
        footward = np.mean([d["Knee"] - d["Hip"] + d["Ankle"] - d["Knee"]
                            for d in lm.values()], axis=0)
        if u @ footward < 0:
            u = -u

        # Body midline point at hip level
        if len(lm) == 2:
            mid = (lm["L"]["Hip"] + lm["R"]["Hip"]) / 2.0
        else:
            H = next(iter(lm.values()))["Hip"]
            slab = body & (np.abs(grid.planeField(H, u)) <= max(1.0, grid.spacing.max()))
            sIdx = np.argwhere(slab)
            mid = grid.ras(sIdx).mean(axis=0) if len(sIdx) else pts.mean(axis=0)

        label = np.zeros(body.shape, dtype=np.uint8)
        codes = {}
        code = 0
        for sideKey, d in lm.items():
            for zone in ZONE_NAMES:
                code += 1
                codes[(sideKey, zone)] = code

        for sideKey, d in lm.items():
            sideName = dict(SIDES)[sideKey]
            H, K, A = d["Hip"], d["Knee"], d["Ankle"]
            f = _unit(K - H)
            t = _unit(A - K)
            tibiaLen = float(np.linalg.norm(A - K))
            nKnee = _unit(f + t) if np.linalg.norm(f + t) > 1e-6 else f
            R = captureRadiusMm if captureRadiusMm > 0 else 0.35 * tibiaLen
            footReach = footReachMm if footReachMm > 0 else tibiaLen

            # -- Foot: beyond the ankle plane, near the ankle -------------
            progress(f"{sideName} foot …")
            box = grid.boxAroundRas(A - footReach, A + footReach)
            cand = self._halfSpaces(grid, body, box, [(A, t)])
            cand = self._withinDistance(grid, cand, box, A, A, footReach)
            self._assign(grid, label, cand, box, A, codes[(sideKey, "Foot")])

            # -- Distal: between knee and ankle planes, near the tibia ----
            progress(f"{sideName} distal …")
            box = grid.boxAroundRas(np.minimum(K, A) - R, np.maximum(K, A) + R)
            cand = self._halfSpaces(grid, body, box, [(K, nKnee), (A, -t)])
            cand = self._withinDistance(grid, cand, box, K, A, R)
            self._assign(grid, label, cand, box, (K + A) / 2.0, codes[(sideKey, "Distal")])

            # -- Proximal: hip level -> knee, lateral to the medial cut ---
            progress(f"{sideName} proximal …")
            lat = (H - mid) - ((H - mid) @ u) * u
            if np.linalg.norm(lat) < 1e-3:
                raise RuntimeError(f"{sideName} hip landmark lies on the body "
                                   "midline - check the landmarks.")
            halfWidth = float(np.linalg.norm(lat))
            lat = lat / halfWidth
            medialPoint = mid + medialFraction * halfWidth * lat
            hipPoint = H - hipMarginMm * u
            box = tuple((0, n) for n in body.shape)
            cand = self._halfSpaces(grid, body, box,
                                    [(hipPoint, u), (K, -nKnee), (medialPoint, lat)])
            self._assign(grid, label, cand, box, (H + K) / 2.0, codes[(sideKey, "Proximal")])

        # -- write segments -----------------------------------------------
        progress("Writing segments …")
        counts = {}
        for (sideKey, zone), c in codes.items():
            name = f"{dict(SIDES)[sideKey]}_{zone}"
            self._removeSegmentByName(segmentationNode, name)
            zoneMask = label == c
            counts[name] = int(zoneMask.sum())
            if counts[name] == 0:
                logging.warning("PADAnalysis2: %s is empty", name)
                continue
            full = np.zeros(bodyFull.shape, dtype=np.uint8)
            full[crop] = zoneMask
            segId = segmentationNode.GetSegmentation().AddEmptySegment(
                "", name, ZONE_COLORS[(sideKey, zone)]
            )
            slicer.util.updateSegmentBinaryLabelmapFromArray(
                full, segmentationNode, segId, workNode
            )

        dn = segmentationNode.GetDisplayNode()
        if dn is not None:
            dn.SetSegmentVisibility(bodyId, False)
        segmentationNode.CreateClosedSurfaceRepresentation()

        if not any(counts.values()):
            raise RuntimeError("All zones are empty - check the landmarks and the body mask.")
        progress("Zones done.")
        return counts

    # -- zone helpers ----------------------------------------------------

    @staticmethod
    def _halfSpaces(grid, body, box, planes):
        """Body voxels in box with (x - origin) . normal >= 0 for all planes."""
        (k0, k1), (j0, j1), (i0, i1) = box
        m = body[k0:k1, j0:j1, i0:i1].copy()
        for origin, normal in planes:
            if not m.any():
                break
            m &= grid.planeField(origin, normal, box) >= 0
        return m

    @staticmethod
    def _withinDistance(grid, m, box, A, B, radius):
        """Remove voxels farther than radius from segment A-B."""
        idx = np.argwhere(m)
        if len(idx) == 0:
            return m
        offset = np.array([b[0] for b in box])
        d = _distanceToSegment(grid.ras(idx + offset), A, B)
        far = idx[d > radius]
        m[far[:, 0], far[:, 1], far[:, 2]] = False
        return m

    @staticmethod
    def _assign(grid, label, cand, box, seedRas, code):
        (k0, k1), (j0, j1), (i0, i1) = box
        sub = label[k0:k1, j0:j1, i0:i1]
        cand &= sub == 0
        if not cand.any():
            return
        subGrid = _Grid(grid.M, cand.shape, offsetIJK=(i0, j0, k0))
        keep = _keepComponentNearest(cand, subGrid, seedRas)
        sub[keep] = code

    # ------------------------------------------------------------------
    # Step 4: quantification
    # ------------------------------------------------------------------

    @staticmethod
    def netDoseBq(syringeActivity, residualActivity, doseUnits, decayCorrect=False,
                  halfLifeMin=None, injectionToScanMin=0.0):
        dose = (syringeActivity - residualActivity) * DOSE_UNITS_TO_BQ[doseUnits]
        if dose <= 0:
            raise ValueError("Net injected activity must be > 0.")
        if decayCorrect:
            if not halfLifeMin or halfLifeMin <= 0:
                raise ValueError("A positive half-life is needed for decay correction.")
            dose *= math.exp(-math.log(2.0) * injectionToScanMin / halfLifeMin)
        return dose

    def quantify(self, segmentationNode, petVolumeNode, petUnits, syringeActivity,
                 residualActivity, doseUnits, bodyWeightG, decayCorrect=True,
                 halfLifeMin=109.77, injectionToScanMin=0.0, ischemicSide="L",
                 skinMarginMm=0.0, subjectID="", outputCsvPath=None,
                 appendToCsv=True, progressCallback=None):
        """Returns dict with 'stats' (list of row dicts), 'ratios', and table nodes."""
        import SimpleITK as sitk
        progress = self._progress(progressCallback)
        self._checkNoTransform(petVolumeNode, "PET")
        if bodyWeightG <= 0:
            raise ValueError("Body weight must be > 0 g.")
        doseBq = self.netDoseBq(syringeActivity, residualActivity, doseUnits,
                                decayCorrect, halfLifeMin, injectionToScanMin)
        toBqml = PET_UNITS_TO_BQML[petUnits]

        workNode = self.workingVolume(segmentationNode)
        if workNode is None:
            raise RuntimeError("Segmentation has no working geometry - run Steps 1 and 3 first.")
        workGrid = _Grid(_ijkToRas(workNode), slicer.util.arrayFromVolume(workNode).shape)
        voxelVolume = float(np.prod(workNode.GetSpacing()))
        petArr = slicer.util.arrayFromVolume(petVolumeNode).astype(np.float64)
        petGrid = _Grid(_ijkToRas(petVolumeNode), petArr.shape)

        stats = []
        for sideKey, sideName in SIDES:
            for zone in ZONE_NAMES:
                name = f"{sideName}_{zone}"
                segId = self._segmentIdByName(segmentationNode, name)
                if segId is None:
                    continue
                progress(f"Quantifying {name} …")
                mask = slicer.util.arrayFromSegmentBinaryLabelmap(
                    segmentationNode, segId, workNode).astype(bool)
                volume = float(mask.sum()) * voxelVolume
                if skinMarginMm > 0 and mask.any():
                    r = int(round(skinMarginMm / float(np.min(workNode.GetSpacing()))))
                    if r > 0:
                        img = sitk.BinaryErode(sitk.GetImageFromArray(mask.astype(np.uint8)), [r] * 3)
                        mask = sitk.GetArrayFromImage(img).astype(bool)
                values = petArr[self._maskOnGrid(mask, workGrid, petGrid)]

                row = {h: "" for h in STATS_HEADER}
                row.update({
                    "SubjectID": subjectID, "Segment": name, "Side": sideName,
                    "Zone": zone,
                    "Condition": "Ischemic" if sideKey == ischemicSide else "Control",
                    "Volume_mm3": volume, "PET_Voxels": int(values.size),
                    "PET_Units": petUnits,
                    "NetDose_MBq": doseBq / 1e6, "BodyWeight_g": bodyWeightG,
                })
                if values.size:
                    mean, sd, mx = values.mean(), values.std(), values.max()
                    pct = 100.0 * toBqml / doseBq            # %ID/g per image unit
                    suv = toBqml * bodyWeightG / doseBq      # SUV per image unit
                    row.update({
                        "Mean": mean, "StdDev": sd, "Median": float(np.median(values)),
                        "Min": values.min(), "Max": mx,
                        "%ID/g_Mean": mean * pct, "%ID/g_StdDev": sd * pct, "%ID/g_Max": mx * pct,
                        "SUV_Mean": mean * suv, "SUV_StdDev": sd * suv, "SUV_Max": mx * suv,
                    })
                else:
                    logging.warning("PADAnalysis2: %s does not overlap the PET volume", name)
                stats.append(row)

        if not stats:
            raise RuntimeError("No zone segments found - run Step 3 first.")

        # Ratios ischemic / control
        controlSide = "R" if ischemicSide == "L" else "L"
        bySeg = {r["Segment"]: r for r in stats}
        ratios = []
        for zone in ZONE_NAMES:
            isch = bySeg.get(f"{dict(SIDES)[ischemicSide]}_{zone}")
            ctrl = bySeg.get(f"{dict(SIDES)[controlSide]}_{zone}")
            if not isch or not ctrl:
                continue

            def ratio(key):
                try:
                    return float(isch[key]) / float(ctrl[key])
                except (ValueError, TypeError, ZeroDivisionError):
                    return float("nan")
            ratios.append({
                "SubjectID": subjectID, "Zone": zone,
                "IschemicSide": dict(SIDES)[ischemicSide],
                "ControlSide": dict(SIDES)[controlSide],
                "Ratio_Mean": ratio("Mean"), "Ratio_Max": ratio("Max"),
            })

        if outputCsvPath:
            self._writeCsv(outputCsvPath, STATS_HEADER, stats, appendToCsv)
            if ratios:
                base, _ = os.path.splitext(outputCsvPath)
                self._writeCsv(base + "_ratios.csv", RATIO_HEADER, ratios, appendToCsv)

        statsTable = self._makeTable("PAD2_Statistics", STATS_HEADER, stats)
        ratioTable = self._makeTable("PAD2_Ratios", RATIO_HEADER, ratios)
        progress("Quantification done.")
        return {"stats": stats, "ratios": ratios,
                "statsTable": statsTable, "ratioTable": ratioTable}

    @staticmethod
    def _maskOnGrid(mask, srcGrid, dstGrid):
        """Nearest-neighbour lookup of mask (srcGrid) at the voxel centres of
        dstGrid. Returns an index tuple into a dstGrid-shaped array."""
        idx = np.argwhere(mask)
        if len(idx) == 0:
            return (np.array([], int),) * 3
        lo, hi = idx.min(axis=0) - 0.5, idx.max(axis=0) + 0.5     # k, j, i
        corners = np.array([[k, j, i] for k in (lo[0], hi[0])
                            for j in (lo[1], hi[1]) for i in (lo[2], hi[2])])
        dstIjk = dstGrid.ijk(srcGrid.ras(corners))
        mn = np.floor(dstIjk.min(axis=0)).astype(int)
        mx = np.ceil(dstIjk.max(axis=0)).astype(int) + 1
        nk, nj, ni = dstGrid.shape
        i_r = np.arange(max(mn[0], 0), min(mx[0], ni))
        j_r = np.arange(max(mn[1], 0), min(mx[1], nj))
        k_r = np.arange(max(mn[2], 0), min(mx[2], nk))
        if not (len(i_r) and len(j_r) and len(k_r)):
            return (np.array([], int),) * 3
        kk, jj, ii = np.meshgrid(k_r, j_r, i_r, indexing="ij")
        kji = np.stack([kk.ravel(), jj.ravel(), ii.ravel()], axis=1)
        src = np.rint(srcGrid.ijk(dstGrid.ras(kji))).astype(int)[:, ::-1]   # k, j, i
        valid = np.all((src >= 0) & (src < np.array(srcGrid.shape)), axis=1)
        inside = np.zeros(len(kji), dtype=bool)
        s = src[valid]
        inside[valid] = mask[s[:, 0], s[:, 1], s[:, 2]]
        sel = kji[inside]
        return sel[:, 0], sel[:, 1], sel[:, 2]

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _writeCsv(path, header, rows, append):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        exists = os.path.isfile(path) and os.path.getsize(path) > 0
        mode = "a" if (append and exists) else "w"
        with open(path, mode, newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=header)
            if mode == "w":
                w.writeheader()
            w.writerows(rows)

    @staticmethod
    def _makeTable(name, header, rows):
        old = slicer.mrmlScene.GetFirstNodeByName(name)
        if old is not None:
            slicer.mrmlScene.RemoveNode(old)
        tableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", name)
        table = tableNode.GetTable()
        for h in header:
            numeric = rows and all(isinstance(r[h], (int, float, np.floating, np.integer))
                                   for r in rows)
            arr = vtk.vtkDoubleArray() if numeric else vtk.vtkStringArray()
            arr.SetName(h)
            for r in rows:
                v = r[h]
                arr.InsertNextValue(float(v) if numeric else
                                    (f"{v:.6g}" if isinstance(v, float) else str(v)))
            table.AddColumn(arr)
        tableNode.Modified()
        return tableNode

    @staticmethod
    def showTable(tableNode):
        lm = slicer.app.layoutManager()
        if lm is None or tableNode is None:
            return
        lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpTableView)
        sel = slicer.app.applicationLogic().GetSelectionNode()
        sel.SetActiveTableID(tableNode.GetID())
        slicer.app.applicationLogic().PropagateTableSelection()

    # ------------------------------------------------------------------
    # Misc helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _progress(callback):
        def progress(msg):
            logging.info("PADAnalysis2: %s", msg)
            if callback:
                callback(msg)
        return progress

    @staticmethod
    def _checkNoTransform(volumeNode, label):
        if volumeNode.GetParentTransformNode() is not None:
            raise RuntimeError(
                f"The {label} volume is under a transform. Harden the transform "
                "(Data module → right-click → Harden transform) first."
            )

    @staticmethod
    def workingVolume(segmentationNode):
        return segmentationNode.GetNodeReference(WORKING_VOLUME_REF)

    @staticmethod
    def _segmentIdByName(segmentationNode, name):
        seg = segmentationNode.GetSegmentation()
        for i in range(seg.GetNumberOfSegments()):
            sid = seg.GetNthSegmentID(i)
            if seg.GetSegment(sid).GetName() == name:
                return sid
        return None

    @staticmethod
    def _removeSegmentByName(segmentationNode, name):
        seg = segmentationNode.GetSegmentation()
        for i in range(seg.GetNumberOfSegments() - 1, -1, -1):
            sid = seg.GetNthSegmentID(i)
            if seg.GetSegment(sid).GetName() == name:
                seg.RemoveSegment(sid)


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

class PADAnalysis2Test(ScriptedLoadableModuleTest):
    """Self-test on a synthetic 'mouse' (torso, two legs, tail, bed)."""

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_FullPipeline()

    @staticmethod
    def _capsule(shape, A, B, radius):
        """Boolean (k, j, i) mask of a capsule between A and B (x, y, z = i, j, k)."""
        nk, nj, ni = shape
        kk, jj, ii = np.meshgrid(np.arange(nk), np.arange(nj), np.arange(ni), indexing="ij")
        P = np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1).astype(float)
        d = _distanceToSegment(P, np.asarray(A, float), np.asarray(B, float))
        return (d <= radius).reshape(shape)

    def test_FullPipeline(self):
        import tempfile
        self.delayDisplay("Building synthetic CT / PET …")
        # 1 mm grid, RAS == IJK. x = left-right, y = AP, z = head(+)/feet(-)
        shape = (110, 34, 70)                       # (z, y, x)
        ct = np.full(shape, -1000, dtype=np.int16)
        body = np.zeros(shape, dtype=bool)
        body[55:105, 6:28, 20:50] = True           # torso
        legs = {
            "L": dict(Hip=(42, 17, 58), Knee=(52, 17, 38), Ankle=(56, 17, 18), Toe=(57, 17, 6)),
            "R": dict(Hip=(28, 17, 58), Knee=(18, 17, 38), Ankle=(14, 17, 18), Toe=(13, 17, 6)),
        }
        for d in legs.values():
            body |= self._capsule(shape, d["Hip"], d["Knee"], 4)
            body |= self._capsule(shape, d["Knee"], d["Ankle"], 2.5)
            body |= self._capsule(shape, d["Ankle"], d["Toe"], 2)
        tail = self._capsule(shape, (35, 17, 56), (35, 17, 8), 1.5)
        body |= tail
        ct[body] = 40
        ct[0:3, :, :] = 200                        # "bed" plate, not touching the animal
        ctNode = slicer.util.addVolumeFromArray(ct, name="Test_CT")

        # PET on a coarser 2 mm grid: left half 2000 Bq/mL, right half 1000
        pet = np.zeros((55, 17, 35), dtype=np.float32)
        petBody = body[::2, ::2, ::2]
        pet[petBody] = 1000.0
        pet[:, :, 18:][petBody[:, :, 18:]] = 2000.0     # x >= 36 mm -> left leg
        petNode = slicer.util.addVolumeFromArray(pet, name="Test_PET")
        petNode.SetSpacing(2.0, 2.0, 2.0)

        logic = PADAnalysis2Logic()
        segNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentationNode")

        self.delayDisplay("Body mask …")
        n = logic.createBodyMask(ctNode, segNode, thresholdHU=-300,
                                 openingRadiusMm=0, workingSpacingMm=0)
        self.assertEqual(n, int(body.sum()))       # bed removed, animal intact

        self.delayDisplay("Landmarks + zones …")
        lmNodes = {}
        for s in ("L", "R"):
            node = logic.createLandmarkNode(s)
            self.assertFalse(logic.landmarksComplete(node))
            for i, name in enumerate(LANDMARKS):
                node.SetNthControlPointPosition(i, *legs[s][name])
            self.assertTrue(logic.landmarksComplete(node))
            lmNodes[s] = logic.landmarkPositions(node)
        counts = logic.generateZones(segNode, lmNodes)
        for name, c in counts.items():
            self.assertGreater(c, 0, f"{name} empty: {counts}")

        # The tail must not belong to any zone (check its distal part)
        lab = np.zeros(shape, dtype=bool)
        for name in counts:
            sid = logic._segmentIdByName(segNode, name)
            lab |= slicer.util.arrayFromSegmentBinaryLabelmap(segNode, sid, ctNode).astype(bool)
        tailOnly = tail.copy()
        tailOnly[50:, :, :] = False
        self.assertFalse((lab & tailOnly).any(), "tail voxels assigned to a zone")
        # Foot zones lie beyond the ankles
        footL = slicer.util.arrayFromSegmentBinaryLabelmap(
            segNode, logic._segmentIdByName(segNode, "Left_Foot"), ctNode)
        self.assertLess(np.argwhere(footL)[:, 0].max(), 20)

        self.delayDisplay("Quantification …")
        csvPath = os.path.join(tempfile.mkdtemp(), "pad2.csv")
        res = logic.quantify(segNode, petNode, petUnits="Bq/mL",
                             syringeActivity=10.5, residualActivity=0.5, doseUnits="MBq",
                             bodyWeightG=25.0, decayCorrect=False, ischemicSide="L",
                             subjectID="synthetic", outputCsvPath=csvPath)
        rows = {r["Segment"]: r for r in res["stats"]}
        self.assertAlmostEqual(rows["Right_Distal"]["Mean"], 1000.0, places=3)
        self.assertAlmostEqual(rows["Left_Distal"]["Mean"], 2000.0, places=3)
        # 2000 Bq/mL with 10 MBq net: %ID/g = 0.02, SUV = 2000*25/1e7 = 0.005
        self.assertAlmostEqual(rows["Left_Distal"]["%ID/g_Mean"], 0.02, places=6)
        self.assertAlmostEqual(rows["Left_Distal"]["SUV_Mean"], 0.005, places=6)
        ratios = {r["Zone"]: r for r in res["ratios"]}
        for zone in ("Distal", "Foot"):
            self.assertAlmostEqual(ratios[zone]["Ratio_Mean"], 2.0, places=3)
        self.assertTrue(os.path.isfile(csvPath))
        self.assertTrue(os.path.isfile(csvPath.replace(".csv", "_ratios.csv")))

        # Decay correction sanity: one half-life halves the dose
        d = logic.netDoseBq(10, 0, "MBq", True, 109.77, 109.77)
        self.assertAlmostEqual(d, 5e6, delta=1)
        self.delayDisplay("Test passed ✓")

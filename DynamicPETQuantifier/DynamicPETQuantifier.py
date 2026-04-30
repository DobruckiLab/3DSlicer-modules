"""
DynamicPETQuantifier - 3D Slicer Scripted Module
=================================================
Quantifies all frames of a dynamic PET sequence against a predefined
segmentation and exports results to a single CSV or Excel file.

Author: Generated for 3D Slicer (Python scripted module)
Slicer version: 5.x+
"""

import os
import csv
import logging
from datetime import datetime

import qt
import ctk
import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleWidget,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
)
from slicer.util import VTKObservationMixin


# ---------------------------------------------------------------------------
# Module descriptor
# ---------------------------------------------------------------------------

class DynamicPETQuantifier(ScriptedLoadableModule):
    """Module metadata shown in the Slicer module finder."""

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Dynamic PET Quantifier"
        self.parent.categories = ["Quantification"]
        self.parent.dependencies = ["Sequences", "SegmentStatistics"]
        self.parent.contributors = ["Wawrzyniec Dobrucki <br> Experimental Molecular Imaging Laboratory <br> University of Illinois at Urbana-Champaign"]
        self.parent.helpText = (
            "Quantifies every frame in a dynamic PET sequence using a "
            "pre-defined segmentation and exports results to CSV or Excel. "
            "<br><br>Select a Sequence Browser node that contains the PET "
            "frames, select the Segmentation node, choose which statistics "
            "to compute, and click <b>Run Quantification</b>."
            "<br><br><b>Proxy Volume:</b> The proxy volume is the scalar volume "
            "node that the Sequence Browser keeps updated as you scrub through "
            "frames. You can identify it in the Sequences module under the "
            "<i>Proxy node</i> column of the Synchronized nodes table."
            "<br><br><b>Excel export requires openpyxl.</b> If you see the "
            "warning <i>'openpyxl not available; writing CSV instead'</i>, "
            "install it via the Slicer Python Interactor "
            "(<b>View &rarr; Python Interactor</b>) using one of these commands:"
            "<br><br>&nbsp;&nbsp;<code>pip_install(\"openpyxl\")</code>"
            "<br><br>If that does not work, try:"
            "<br>&nbsp;&nbsp;<code>import subprocess, sys</code>"
            "<br>&nbsp;&nbsp;<code>subprocess.check_call([sys.executable, \"-m\", \"pip\", \"install\", \"openpyxl\"])</code>"
            "<br><br>No restart is required after installation."
        )
        self.parent.acknowledgementText = (
            ""
        )


# ---------------------------------------------------------------------------
# Widget (GUI)
# ---------------------------------------------------------------------------

class DynamicPETQuantifierWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """User interface for the module."""

    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None

    # ------------------------------------------------------------------
    # setup
    # ------------------------------------------------------------------

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = DynamicPETQuantifierLogic()

        # ── Input section ─────────────────────────────────────────────
        inputCollapsible = ctk.ctkCollapsibleButton()
        inputCollapsible.text = "Inputs"
        self.layout.addWidget(inputCollapsible)
        inputForm = qt.QFormLayout(inputCollapsible)

        # Sequence Browser selector
        self.sequenceBrowserSelector = slicer.qMRMLNodeComboBox()
        self.sequenceBrowserSelector.nodeTypes = ["vtkMRMLSequenceBrowserNode"]
        self.sequenceBrowserSelector.selectNodeUponCreation = True
        self.sequenceBrowserSelector.addEnabled = False
        self.sequenceBrowserSelector.removeEnabled = False
        self.sequenceBrowserSelector.noneEnabled = True
        self.sequenceBrowserSelector.showHidden = False
        self.sequenceBrowserSelector.setMRMLScene(slicer.mrmlScene)
        self.sequenceBrowserSelector.setToolTip(
            "Select the Sequence Browser that contains the dynamic PET frames."
        )
        inputForm.addRow("Sequence Browser:", self.sequenceBrowserSelector)

        # Proxy Volume selector (the PET scalar volume inside the sequence)
        self.proxyVolumeSelector = slicer.qMRMLNodeComboBox()
        self.proxyVolumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.proxyVolumeSelector.selectNodeUponCreation = False
        self.proxyVolumeSelector.addEnabled = False
        self.proxyVolumeSelector.removeEnabled = False
        self.proxyVolumeSelector.noneEnabled = True
        self.proxyVolumeSelector.showHidden = False
        self.proxyVolumeSelector.setMRMLScene(slicer.mrmlScene)
        self.proxyVolumeSelector.setToolTip(
            "Select the PET proxy volume node (scalar volume that is "
            "synchronized with the sequence browser)."
        )
        inputForm.addRow("PET Proxy Volume:", self.proxyVolumeSelector)

        # Segmentation selector
        self.segmentationSelector = slicer.qMRMLNodeComboBox()
        self.segmentationSelector.nodeTypes = ["vtkMRMLSegmentationNode"]
        self.segmentationSelector.selectNodeUponCreation = False
        self.segmentationSelector.addEnabled = False
        self.segmentationSelector.removeEnabled = False
        self.segmentationSelector.noneEnabled = True
        self.segmentationSelector.showHidden = False
        self.segmentationSelector.setMRMLScene(slicer.mrmlScene)
        self.segmentationSelector.setToolTip(
            "Select the segmentation that contains all organ/ROI segments."
        )
        inputForm.addRow("Segmentation:", self.segmentationSelector)

        # Time step / frame duration (optional, in seconds)
        self.frameDurationSpinBox = qt.QDoubleSpinBox()
        self.frameDurationSpinBox.setRange(0.0, 100000.0)
        self.frameDurationSpinBox.setValue(0.0)
        self.frameDurationSpinBox.setSuffix(" s")
        self.frameDurationSpinBox.setDecimals(2)
        self.frameDurationSpinBox.setToolTip(
            "Optional: uniform frame duration in seconds used to compute "
            "mid-frame times. Leave 0 to use frame index instead of time."
        )
        inputForm.addRow("Frame Duration (uniform):", self.frameDurationSpinBox)

        # ── Statistics options ─────────────────────────────────────────
        statsCollapsible = ctk.ctkCollapsibleButton()
        statsCollapsible.text = "Statistics to Compute"
        self.layout.addWidget(statsCollapsible)
        statsForm = qt.QFormLayout(statsCollapsible)

        self.cbCount   = self._makeCheckBox("Voxel Count",        True,  statsForm)
        self.cbVolume  = self._makeCheckBox("Volume (mm³)",        True,  statsForm)
        self.cbMin     = self._makeCheckBox("Minimum",             True,  statsForm)
        self.cbMax     = self._makeCheckBox("Maximum",             True,  statsForm)
        self.cbMean    = self._makeCheckBox("Mean",                True,  statsForm)
        self.cbStd     = self._makeCheckBox("Standard Deviation",  True,  statsForm)
        self.cbMedian  = self._makeCheckBox("Median",              False, statsForm)

        # ── Output section ─────────────────────────────────────────────
        outputCollapsible = ctk.ctkCollapsibleButton()
        outputCollapsible.text = "Output"
        self.layout.addWidget(outputCollapsible)
        outputForm = qt.QFormLayout(outputCollapsible)

        # Output format
        self.formatCombo = qt.QComboBox()
        self.formatCombo.addItems(["CSV (.csv)", "Excel (.xlsx)"])
        outputForm.addRow("Output Format:", self.formatCombo)

        # Output file path
        fileRow = qt.QHBoxLayout()
        self.outputPathEdit = qt.QLineEdit()
        self.outputPathEdit.setPlaceholderText("Choose output file path …")
        self.outputPathEdit.setToolTip("Full path for the output file.")
        browseBtn = qt.QPushButton("Browse…")
        browseBtn.clicked.connect(self.onBrowseOutput)
        fileRow.addWidget(self.outputPathEdit)
        fileRow.addWidget(browseBtn)
        outputForm.addRow("Output File:", fileRow)

        # Open after export
        self.openAfterExportCheckBox = qt.QCheckBox("Open file after export")
        self.openAfterExportCheckBox.checked = False
        outputForm.addRow("", self.openAfterExportCheckBox)

        # ── Run button ─────────────────────────────────────────────────
        self.runButton = qt.QPushButton("▶  Run Quantification")
        self.runButton.toolTip = "Quantify all frames and save results."
        self.runButton.enabled = False
        self.runButton.setStyleSheet(
            "QPushButton { background-color: #2ecc71; color: white; "
            "font-weight: bold; padding: 8px; border-radius: 4px; }"
            "QPushButton:disabled { background-color: #95a5a6; }"
            "QPushButton:hover { background-color: #27ae60; }"
        )
        self.layout.addWidget(self.runButton)

        # Progress bar
        self.progressBar = qt.QProgressBar()
        self.progressBar.setRange(0, 100)
        self.progressBar.setValue(0)
        self.progressBar.setVisible(False)
        self.layout.addWidget(self.progressBar)

        # Status label
        self.statusLabel = qt.QLabel("")
        self.statusLabel.setStyleSheet("color: #555; font-style: italic;")
        self.layout.addWidget(self.statusLabel)

        self.layout.addStretch()

        # Connections
        self.sequenceBrowserSelector.currentNodeChanged.connect(self.onInputChanged)
        self.proxyVolumeSelector.currentNodeChanged.connect(self.onInputChanged)
        self.segmentationSelector.currentNodeChanged.connect(self.onInputChanged)
        self.outputPathEdit.textChanged.connect(self.onInputChanged)
        self.runButton.clicked.connect(self.onRunClicked)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _makeCheckBox(self, label, default, layout):
        cb = qt.QCheckBox(label)
        cb.checked = default
        layout.addRow("", cb)
        return cb

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def onBrowseOutput(self):
        ext = "CSV (*.csv)" if self.formatCombo.currentIndex == 0 else "Excel (*.xlsx)"
        path = qt.QFileDialog.getSaveFileName(
            self.parent, "Save Output File", "", ext
        )
        if path:
            self.outputPathEdit.text = path

    def onInputChanged(self, *args):
        ready = (
            self.sequenceBrowserSelector.currentNode() is not None
            and self.proxyVolumeSelector.currentNode() is not None
            and self.segmentationSelector.currentNode() is not None
            and self.outputPathEdit.text.strip() != ""
        )
        self.runButton.enabled = ready

    def onRunClicked(self):
        browserNode    = self.sequenceBrowserSelector.currentNode()
        proxyVolume    = self.proxyVolumeSelector.currentNode()
        segmentationNode = self.segmentationSelector.currentNode()
        outputPath     = self.outputPathEdit.text.strip()
        frameDuration  = self.frameDurationSpinBox.value
        useExcel       = (self.formatCombo.currentIndex == 1)

        # Ensure correct extension
        if useExcel and not outputPath.lower().endswith(".xlsx"):
            outputPath += ".xlsx"
        elif not useExcel and not outputPath.lower().endswith(".csv"):
            outputPath += ".csv"
        self.outputPathEdit.text = outputPath

        statsFlags = {
            "voxel_count":  self.cbCount.checked,
            "volume_mm3":   self.cbVolume.checked,
            "minimum":      self.cbMin.checked,
            "maximum":      self.cbMax.checked,
            "mean":         self.cbMean.checked,
            "stdev":        self.cbStd.checked,
            "median":       self.cbMedian.checked,
        }

        self.runButton.enabled = False
        self.progressBar.setVisible(True)
        self.progressBar.setValue(0)
        self.statusLabel.text = "Initialising …"
        slicer.app.processEvents()

        try:
            self.logic.run(
                browserNode=browserNode,
                proxyVolume=proxyVolume,
                segmentationNode=segmentationNode,
                outputPath=outputPath,
                frameDurationSec=frameDuration,
                statsFlags=statsFlags,
                useExcel=useExcel,
                progressCallback=self.updateProgress,
            )
            self.statusLabel.text = f"✔ Done! Results saved to: {outputPath}"
            self.statusLabel.setStyleSheet("color: #27ae60; font-weight: bold;")
            if self.openAfterExportCheckBox.checked:
                qt.QDesktopServices.openUrl(qt.QUrl.fromLocalFile(outputPath))
        except Exception as e:
            slicer.util.errorDisplay(f"Quantification failed:\n{e}")
            self.statusLabel.text = f"✘ Error: {e}"
            self.statusLabel.setStyleSheet("color: #e74c3c; font-weight: bold;")
            logging.exception("DynamicPETQuantifier error")
        finally:
            self.progressBar.setVisible(False)
            self.onInputChanged()

    def updateProgress(self, value, message=""):
        self.progressBar.setValue(int(value))
        if message:
            self.statusLabel.text = message
        slicer.app.processEvents()

    def cleanup(self):
        pass


# ---------------------------------------------------------------------------
# Logic
# ---------------------------------------------------------------------------

class DynamicPETQuantifierLogic(ScriptedLoadableModuleLogic):
    """
    Core computation logic — no GUI dependencies.
    Can be used from Python scripts / Jupyter without the widget.
    """

    def run(
        self,
        browserNode,
        proxyVolume,
        segmentationNode,
        outputPath,
        frameDurationSec=0.0,
        statsFlags=None,
        useExcel=False,
        progressCallback=None,
    ):
        """
        Main entry point.

        Parameters
        ----------
        browserNode       : vtkMRMLSequenceBrowserNode
        proxyVolume       : vtkMRMLScalarVolumeNode  (proxy in the browser)
        segmentationNode  : vtkMRMLSegmentationNode
        outputPath        : str  — full path to output file
        frameDurationSec  : float — uniform frame duration in seconds
                            (0 = use frame index as "time")
        statsFlags        : dict[str, bool]  — which stats to include
        useExcel          : bool
        progressCallback  : callable(percent, message) | None
        """
        import SegmentStatistics
        import tempfile

        if statsFlags is None:
            statsFlags = {
                "voxel_count": True, "volume_mm3": True,
                "minimum": True, "maximum": True,
                "mean": True, "stdev": True, "median": False,
            }

        def progress(pct, msg=""):
            if progressCallback:
                progressCallback(pct, msg)

        # ── Identify sequence node for the proxy volume ────────────────
        sequenceNode = self._getSequenceNodeForProxy(browserNode, proxyVolume)
        if sequenceNode is None:
            raise RuntimeError(
                "Could not find a Sequence node linked to the selected proxy "
                "volume in the chosen Sequence Browser. "
                "Make sure the proxy volume belongs to the selected browser."
            )

        nFrames = sequenceNode.GetNumberOfDataNodes()
        if nFrames == 0:
            raise RuntimeError("The sequence contains 0 frames.")

        progress(2, f"Found {nFrames} frames. Collecting segment info ...")

        # ── Collect segment names in order ─────────────────────────────
        segmentation = segmentationNode.GetSegmentation()
        segmentIDs = [
            segmentation.GetNthSegmentID(i)
            for i in range(segmentation.GetNumberOfSegments())
        ]
        if not segmentIDs:
            raise RuntimeError("The segmentation contains no segments.")
        segmentNames = [segmentation.GetSegment(sid).GetName() for sid in segmentIDs]

        # ── Save current browser state ─────────────────────────────────
        savedFrame = browserNode.GetSelectedItemNumber()

        # ── Stat column labels we want to keep from the CSV ───────────
        # These match the column headers that SegmentStatistics writes.
        # We build a whitelist based on the user's statsFlags choices.
        STAT_CSV_SUBSTRINGS = {
            "voxel_count": ["voxel count", "voxelcount", "voxel_count", "number of voxels"],
            "volume_mm3":  ["volume mm3", "volume_mm3", "volume [mm3]", "volume[mm3]"],
            "minimum":     ["minimum", "min"],
            "maximum":     ["maximum", "max"],
            "mean":        ["mean"],
            "stdev":       ["standard deviation", "stdev", "std"],
            "median":      ["median"],
        }
        enabledStatKeys = [k for k, v in statsFlags.items() if v]

        # ── Process frame 0 to discover CSV column structure ──────────
        progress(3, "Probing frame 0 to detect column structure ...")
        browserNode.SetSelectedItemNumber(0)
        slicer.app.processEvents()

        frameCsvPath = self._computeFrameStats(
            segmentationNode, proxyVolume, tempfile.mktemp(suffix=".csv")
        )
        with open(frameCsvPath, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            probeHeaders = reader.fieldnames or []
            probeRows = list(reader)

        logging.info("DynamicPETQuantifier: CSV headers: %s", probeHeaders)

        # Map each enabled stat to the actual CSV column name(s)
        # One stat may appear as multiple columns (one per segment).
        # CSV rows are one-per-segment; columns are stat names.
        # Find which CSV column corresponds to each of our stat keys.
        resolvedCsvCols = {}  # our statKey -> actual CSV column header string
        for statKey in enabledStatKeys:
            substrings = STAT_CSV_SUBSTRINGS.get(statKey, [])
            for col in probeHeaders:
                colLower = col.strip().lower()
                if any(sub in colLower for sub in substrings):
                    resolvedCsvCols[statKey] = col
                    break
            if statKey not in resolvedCsvCols:
                logging.warning("DynamicPETQuantifier: no CSV column found for '%s'", statKey)

        # Find which CSV column identifies the segment name
        segNameCol = None
        for col in probeHeaders:
            if col.strip().lower() in ("segment", "name", "segment name"):
                segNameCol = col
                break
        if segNameCol is None and probeHeaders:
            segNameCol = probeHeaders[0]  # fall back to first column
        logging.info("DynamicPETQuantifier: resolved columns: %s", resolvedCsvCols)
        logging.info("DynamicPETQuantifier: segment name column: %s", segNameCol)

        # ── Build final output headers ─────────────────────────────────
        STAT_LABELS = {
            "voxel_count": "VoxelCount", "volume_mm3": "Volume_mm3",
            "minimum": "Min", "maximum": "Max",
            "mean": "Mean", "stdev": "StdDev", "median": "Median",
        }
        headers = ["Frame", "Time_s"]
        for segName in segmentNames:
            for statKey in enabledStatKeys:
                if statKey in resolvedCsvCols:
                    headers.append(f"{segName}_{STAT_LABELS[statKey]}")

        rows = []

        # ── Loop over all frames ───────────────────────────────────────
        for frameIdx in range(nFrames):
            pct = 5 + int(90 * frameIdx / nFrames)
            progress(pct, f"Processing frame {frameIdx + 1} / {nFrames} ...")

            browserNode.SetSelectedItemNumber(frameIdx)
            slicer.app.processEvents()

            timeSec = (
                (frameIdx + 0.5) * frameDurationSec
                if frameDurationSec > 0
                else float(frameIdx)
            )

            # Export stats for this frame to a temp CSV and read it back
            tmpCsv = tempfile.mktemp(suffix=".csv")
            frameCsvPath = self._computeFrameStats(
                segmentationNode, proxyVolume, tmpCsv
            )
            with open(frameCsvPath, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                frameData = {
                    r[segNameCol].strip(): r
                    for r in reader
                    if segNameCol and segNameCol in r
                }

            # Build row: one value per (segment × stat)
            row = [frameIdx, timeSec]
            for segName in segmentNames:
                segRow = frameData.get(segName, {})
                for statKey in enabledStatKeys:
                    if statKey in resolvedCsvCols:
                        rawVal = segRow.get(resolvedCsvCols[statKey], "")
                        # Convert to float where possible
                        try:
                            row.append(float(rawVal))
                        except (ValueError, TypeError):
                            row.append(rawVal)
            rows.append(row)

            # Clean up temp file
            try:
                os.remove(frameCsvPath)
            except OSError:
                pass

        # Restore original frame
        browserNode.SetSelectedItemNumber(savedFrame)

        progress(96, "Writing output file ...")
        self._writeOutput(headers, rows, outputPath, useExcel)
        progress(100, "Done.")

        logging.info(
            "DynamicPETQuantifier: wrote %d frames x %d columns to %s",
            nFrames, len(headers), outputPath,
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _computeFrameStats(self, segmentationNode, scalarVolumeNode, csvPath):
        """
        Run SegmentStatistics for the current frame and export to a CSV file.
        Returns the path to the written CSV.
        """
        import SegmentStatistics
        ssLogic = SegmentStatistics.SegmentStatisticsLogic()
        pNode = ssLogic.getParameterNode()

        # Set inputs
        pNode.SetParameter("Segmentation", segmentationNode.GetID())
        pNode.SetParameter("ScalarVolume", scalarVolumeNode.GetID())

        # Explicitly enable both plugins — they are disabled by default
        # when a fresh logic/parameter node is created
        pNode.SetParameter("LabelmapSegmentStatisticsPlugin.enabled", "True")
        pNode.SetParameter("ScalarVolumeSegmentStatisticsPlugin.enabled", "True")

        # Explicitly enable every individual measurement within each plugin
        for meas in ["voxel_count", "volume_mm3", "volume_cm3"]:
            pNode.SetParameter(
                f"LabelmapSegmentStatisticsPlugin.{meas}.enabled", "True"
            )
        for meas in ["voxel_count", "volume_mm3",
                     "min", "max", "mean", "median", "stdev"]:
            pNode.SetParameter(
                f"ScalarVolumeSegmentStatisticsPlugin.{meas}.enabled", "True"
            )

        ssLogic.computeStatistics()
        ssLogic.exportToCSVFile(csvPath)
        return csvPath

    def _getSequenceNodeForProxy(self, browserNode, proxyVolume):
        """Return the vtkMRMLSequenceNode whose proxy is proxyVolume."""
        proxyID = proxyVolume.GetID()
        for i in range(browserNode.GetNumberOfSynchronizedSequenceNodes()):
            seqNode = browserNode.GetSynchronizedSequenceNode(i)
            proxy   = browserNode.GetProxyNode(seqNode)
            if proxy and proxy.GetID() == proxyID:
                return seqNode
        # Also check the master sequence node
        masterSeq = browserNode.GetMasterSequenceNode()
        if masterSeq:
            proxy = browserNode.GetProxyNode(masterSeq)
            if proxy and proxy.GetID() == proxyID:
                return masterSeq
        return None

    def _writeOutput(self, headers, rows, outputPath, useExcel):
        """Write results to CSV or Excel."""
        os.makedirs(os.path.dirname(os.path.abspath(outputPath)), exist_ok=True)

        if useExcel:
            self._writeExcel(headers, rows, outputPath)
        else:
            self._writeCsv(headers, rows, outputPath)

    def _writeCsv(self, headers, rows, outputPath):
        with open(outputPath, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(rows)

    def _writeExcel(self, headers, rows, outputPath):
        """Write using openpyxl (ships with Slicer 5.x)."""
        try:
            import openpyxl
            from openpyxl.styles import Font, PatternFill, Alignment
        except ImportError:
            # Fall back to CSV if openpyxl is missing
            logging.warning(
                "openpyxl not available; writing CSV instead of Excel."
            )
            csvPath = outputPath.replace(".xlsx", ".csv")
            self._writeCsv(headers, rows, csvPath)
            return

        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "PET Quantification"

        # Header styling
        headerFont  = Font(bold=True, color="FFFFFF")
        headerFill  = PatternFill("solid", fgColor="2F5597")
        headerAlign = Alignment(horizontal="center", wrap_text=True)

        for col, h in enumerate(headers, start=1):
            cell = ws.cell(row=1, column=col, value=h)
            cell.font   = headerFont
            cell.fill   = headerFill
            cell.alignment = headerAlign

        # Data rows — alternate row shading
        fillEven = PatternFill("solid", fgColor="DCE6F1")
        for r, row in enumerate(rows, start=2):
            shade = fillEven if r % 2 == 0 else None
            for col, val in enumerate(row, start=1):
                cell = ws.cell(row=r, column=col, value=val)
                if shade:
                    cell.fill = shade

        # Auto-fit column widths (approximate)
        for col_cells in ws.columns:
            max_len = max(
                len(str(cell.value)) if cell.value is not None else 0
                for cell in col_cells
            )
            ws.column_dimensions[
                col_cells[0].column_letter
            ].width = min(max_len + 4, 40)

        # Freeze header row
        ws.freeze_panes = "A2"

        wb.save(outputPath)


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

class DynamicPETQuantifierTest(ScriptedLoadableModuleTest):
    """Basic self-test — runs without real PET data."""

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_LogicInstantiation()

    def test_LogicInstantiation(self):
        self.delayDisplay("Testing logic instantiation …")
        logic = DynamicPETQuantifierLogic()
        self.assertIsNotNone(logic)
        self.delayDisplay("OK — logic instantiated.")
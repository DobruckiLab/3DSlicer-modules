"""
WorkflowModule.py  –  3D Slicer Scripted Module
================================================
A visual workflow orchestrator that lets you define pipelines of processing
steps (Registration, Filtering, Segmentation, Quantification, etc.) with
pre-saved parameters, then run them on single volumes or entire folders.

Installation:
  1. Copy the SlicerWorkflow/ folder anywhere on disk.
  2. In Slicer: Edit → Application Settings → Modules → Additional module paths
     → add the SlicerWorkflow/ folder.
  3. Restart Slicer. The module appears under "Workflow" category.
"""

import os
import json
import logging
import qt
import ctk
import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleWidget,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
)

# ---------------------------------------------------------------------------
# Module descriptor
# ---------------------------------------------------------------------------

class WorkflowModule(ScriptedLoadableModule):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent.title       = "Workflow"
        self.parent.categories  = ["Utilities", "Workflow"]
        self.parent.dependencies = []
        self.parent.contributors = [
            "<b>Wawrzyniec Dobrucki</b><br>"
            "Experimental Molecular Imaging Laboratory (EMIL)<br>"
            "University of Illinois at Urbana-Champaign"
        ]
        self.parent.helpText    = (
            "<p>Build and run multi-step image-analysis pipelines with "
            "pre-configured parameters across single volumes or entire folders.</p>"
            "<p><b>Tips:</b></p>"
            "<ul>"
            "<li><b>Batch processing:</b> make sure the output directory is set before running, "
            "otherwise the Save Volume and statistics steps will silently skip.</li>"
            "<li><b>Step ordering:</b> segmentation steps (Threshold, Otsu) must come before "
            "Volume/Intensity Statistics, since those steps depend on a segmentation "
            "node being present in the pipeline state.</li>"
            "<li><b>Custom Python Script step:</b> a good way to add any quick one-off processing "
            "without needing a full new step entry.</li>"
            "<li><b>Workflow JSON files:</b> stored in <tt>~/.config/slicer.org/SlicerWorkflows/</tt> on Linux, "
            "so you can back them up or share them with colleagues directly.</li>"
            "</ul>"
        )
        self.parent.acknowledgementText = (
            "<p>Developed by:</p>"
            "<p><b>Wawrzyniec Dobrucki</b><br>"
            "Experimental Molecular Imaging Laboratory (EMIL)<br>"
            "University of Illinois at Urbana-Champaign</p>"
            "<p>This module was created to enable reproducible, parameter-controlled "
            "image analysis pipelines within 3D Slicer, supporting batch processing "
            "of preclinical and clinical imaging data.</p>"
        )


# ---------------------------------------------------------------------------
# Widget  (GUI)
# ---------------------------------------------------------------------------

class WorkflowModuleWidget(ScriptedLoadableModuleWidget):

    # ---- built-in step catalogue -------------------------------------------
    STEP_CATALOGUE = {
        "Registration": {
            "description": "Rigid / affine / deformable registration",
            "params": {
                "transformType": {
                    "type": "enum",
                    "values": ["Rigid", "Affine", "BSpline"],
                    "default": "Rigid",
                },
                "samplingPercentage": {
                    "type": "float",
                    "min": 0.01, "max": 1.0, "default": 0.25,
                },
                "numberOfIterations": {
                    "type": "int",
                    "min": 1, "max": 2000, "default": 200,
                },
                "fixedVolumeNode": {
                    "type": "node",
                    "nodeType": "vtkMRMLScalarVolumeNode",
                    "default": "",
                },
            },
        },
        "N4 Bias Correction": {
            "description": "MRI bias field correction (SimpleITK N4ITK)",
            "params": {
                "numberOfIterations": {
                    "type": "int", "min": 1, "max": 200, "default": 50,
                },
                "convergenceThreshold": {
                    "type": "float", "min": 1e-7, "max": 0.01, "default": 0.001,
                },
                "shrinkFactor": {
                    "type": "int", "min": 1, "max": 8, "default": 4,
                },
            },
        },
        "Gaussian Smoothing": {
            "description": "SimpleITK Gaussian blur / noise reduction",
            "params": {
                "sigma": {
                    "type": "float", "min": 0.1, "max": 10.0, "default": 1.0,
                },
            },
        },
        "Median Filtering": {
            "description": "Median filter (edge-preserving denoising)",
            "params": {
                "radius": {
                    "type": "int", "min": 1, "max": 10, "default": 1,
                },
            },
        },
        "Threshold Segmentation": {
            "description": "Create a segmentation via intensity thresholding",
            "params": {
                "lowerThreshold": {
                    "type": "float", "min": -9999, "max": 9999, "default": 100.0,
                },
                "upperThreshold": {
                    "type": "float", "min": -9999, "max": 9999, "default": 3000.0,
                },
                "segmentName": {
                    "type": "string", "default": "Segment_1",
                },
            },
        },
        "Otsu Segmentation": {
            "description": "Automatic Otsu threshold segmentation",
            "params": {
                "insideValue": {
                    "type": "int", "min": 0, "max": 1, "default": 1,
                },
                "outsideValue": {
                    "type": "int", "min": 0, "max": 1, "default": 0,
                },
            },
        },
        "Volume Statistics": {
            "description": "Compute volume / voxel statistics from a segmentation",
            "params": {
                "saveCSV": {
                    "type": "bool", "default": True,
                },
                "outputDir": {
                    "type": "path", "default": "",
                },
            },
        },
        "Intensity Statistics": {
            "description": "Mean / SD / min / max intensity inside a segment",
            "params": {
                "saveCSV": {
                    "type": "bool", "default": True,
                },
                "outputDir": {
                    "type": "path", "default": "",
                },
            },
        },
        "Save Volume": {
            "description": "Export the current volume to disk",
            "params": {
                "outputDir": {
                    "type": "path", "default": "",
                },
                "fileFormat": {
                    "type": "enum",
                    "values": ["nrrd", "nii.gz", "mha"],
                    "default": "nrrd",
                },
            },
        },
        "Custom Python Script": {
            "description": "Run an arbitrary Python snippet on the current volume node",
            "params": {
                "script": {
                    "type": "text",
                    "default": "# volumeNode and segmentationNode are available\nprint(volumeNode.GetName())",
                },
            },
        },
    }

    def setup(self):
        super().setup()
        self.logic = WorkflowModuleLogic()
        self._buildUI()
        self._connectSignals()
        self._refreshWorkflowList()

    # ------------------------------------------------------------------ UI --

    def _buildUI(self):
        """Build the complete module GUI."""
        # Use the layout Slicer's base class already created on self.parent.
        # Creating a second QVBoxLayout(self.parent) causes the Qt layout conflict error.
        layout = self.layout
        layout.setSpacing(6)
        layout.setContentsMargins(8, 8, 8, 8)

        # Wrap everything in a scroll area so nothing gets clipped
        scrollArea = qt.QScrollArea()
        scrollArea.setWidgetResizable(True)
        scrollWidget = qt.QWidget()
        scrollLayout = qt.QVBoxLayout(scrollWidget)
        scrollLayout.setSpacing(6)
        scrollLayout.setContentsMargins(4, 4, 4, 4)
        scrollArea.setWidget(scrollWidget)
        layout.addWidget(scrollArea)

        # All subsequent widgets are added to scrollLayout, not layout directly
        layout = scrollLayout

        # ── TOP BAR: workflow selector + management ──────────────────────
        topBar = qt.QHBoxLayout()
        self.workflowCombo = qt.QComboBox()
        self.workflowCombo.setSizePolicy(qt.QSizePolicy.Expanding, qt.QSizePolicy.Fixed)
        self.newWorkflowBtn  = self._iconBtn("New",    "New workflow",        "#27ae60")
        self.saveWorkflowBtn = self._iconBtn("Save",   "Save workflow",        "#2980b9")
        self.delWorkflowBtn  = self._iconBtn("Delete", "Delete workflow",      "#c0392b")
        self.importBtn       = self._iconBtn("Import", "Import workflow JSON", "#8e44ad")
        self.exportBtn       = self._iconBtn("Export", "Export workflow JSON", "#16a085")
        for w in [self.workflowCombo, self.newWorkflowBtn, self.saveWorkflowBtn,
                  self.delWorkflowBtn, self.importBtn, self.exportBtn]:
            topBar.addWidget(w)
        layout.addLayout(topBar)

        # ── PIPELINE EDITOR ──────────────────────────────────────────────
        pipelineBox = ctk.ctkCollapsibleButton()
        pipelineBox.text = "Pipeline Steps"
        pipelineBox.collapsed = False
        layout.addWidget(pipelineBox)
        pipelineLayout = qt.QVBoxLayout(pipelineBox)

        # Step list
        self.stepList = qt.QListWidget()
        self.stepList.setDragDropMode(qt.QAbstractItemView.InternalMove)
        self.stepList.setDefaultDropAction(qt.Qt.MoveAction)
        self.stepList.setMinimumHeight(140)
        self.stepList.setMaximumHeight(200)
        self.stepList.setToolTip("Drag rows to re-order. Click to edit parameters.")
        pipelineLayout.addWidget(self.stepList)

        # Step add / remove controls
        stepCtrl = qt.QHBoxLayout()
        self.stepTypeCombo = qt.QComboBox()
        for name in self.STEP_CATALOGUE:
            self.stepTypeCombo.addItem(name)
        self.addStepBtn = qt.QPushButton("Add Step")
        self.addStepBtn.setStyleSheet("background:#27ae60;color:white;font-weight:bold;padding:4px 12px;border-radius:4px;")
        self.removeStepBtn = qt.QPushButton("Remove Selected")
        self.removeStepBtn.setStyleSheet("background:#c0392b;color:white;padding:4px 12px;border-radius:4px;")
        self.moveUpBtn   = qt.QPushButton("▲")
        self.moveDownBtn = qt.QPushButton("▼")
        for w in [self.stepTypeCombo, self.addStepBtn, self.removeStepBtn,
                  self.moveUpBtn, self.moveDownBtn]:
            stepCtrl.addWidget(w)
        pipelineLayout.addLayout(stepCtrl)

        # Parameter panel (dynamic)
        self.paramGroupBox = ctk.ctkCollapsibleButton()
        self.paramGroupBox.text = "Step Parameters"
        self.paramGroupBox.collapsed = False
        pipelineLayout.addWidget(self.paramGroupBox)
        self.paramLayout = qt.QFormLayout(self.paramGroupBox)
        self._paramWidgets = {}   # label → widget

        # ── INPUT / BATCH ────────────────────────────────────────────────
        inputBox = ctk.ctkCollapsibleButton()
        inputBox.text = "Input"
        inputBox.collapsed = False
        layout.addWidget(inputBox)
        inputLayout = qt.QVBoxLayout(inputBox)

        modeRow = qt.QHBoxLayout()
        self.singleRadio = qt.QRadioButton("Single Volume (scene node)")
        self.batchRadio  = qt.QRadioButton("Batch Folder")
        self.singleRadio.setChecked(True)
        modeRow.addWidget(self.singleRadio)
        modeRow.addWidget(self.batchRadio)
        inputLayout.addLayout(modeRow)

        # single node selector
        self.singleNodeSelector = slicer.qMRMLNodeComboBox()
        self.singleNodeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.singleNodeSelector.selectNodeUponCreation = True
        self.singleNodeSelector.addEnabled = False
        self.singleNodeSelector.removeEnabled = False
        self.singleNodeSelector.noneEnabled  = True
        self.singleNodeSelector.setMRMLScene(slicer.mrmlScene)
        self.singleNodeSelector.setToolTip("Choose the volume to process")
        inputLayout.addWidget(self.singleNodeSelector)

        # batch folder
        self.batchWidget = qt.QWidget()
        batchLayout = qt.QHBoxLayout(self.batchWidget)
        batchLayout.setContentsMargins(0, 0, 0, 0)
        self.folderEdit = qt.QLineEdit()
        self.folderEdit.setPlaceholderText("Select input folder …")
        self.browseFolderBtn = qt.QPushButton("Browse…")
        batchLayout.addWidget(self.folderEdit)
        batchLayout.addWidget(self.browseFolderBtn)
        self.batchWidget.setVisible(False)
        inputLayout.addWidget(self.batchWidget)

        # file filter
        filterRow = qt.QHBoxLayout()
        filterRow.addWidget(qt.QLabel("File filter:"))
        self.fileFilterEdit = qt.QLineEdit("*.nrrd *.nii *.nii.gz *.mha *.mhd")
        filterRow.addWidget(self.fileFilterEdit)
        self.batchFilterWidget = qt.QWidget()
        self.batchFilterWidget.setLayout(filterRow)
        self.batchFilterWidget.setVisible(False)
        inputLayout.addWidget(self.batchFilterWidget)

        # ── OUTPUT ───────────────────────────────────────────────────────
        outputBox = ctk.ctkCollapsibleButton()
        outputBox.text = "Output"
        outputBox.collapsed = False
        layout.addWidget(outputBox)
        outputLayout = qt.QFormLayout(outputBox)

        self.outputDirEdit = qt.QLineEdit()
        self.outputDirEdit.setPlaceholderText("(optional) default output folder …")
        self.browseOutputBtn = qt.QPushButton("Browse…")
        outRow = qt.QHBoxLayout()
        outRow.addWidget(self.outputDirEdit)
        outRow.addWidget(self.browseOutputBtn)
        outputLayout.addRow("Output dir:", outRow)

        self.keepIntermediateCheck = qt.QCheckBox("Keep intermediate volumes in scene")
        self.keepIntermediateCheck.setChecked(False)
        outputLayout.addRow("", self.keepIntermediateCheck)

        # ── RUN ──────────────────────────────────────────────────────────
        self.runBtn = qt.QPushButton("▶  Run Workflow")
        self.runBtn.setStyleSheet(
            "background:qlineargradient(x1:0,y1:0,x2:1,y2:0,"
            "stop:0 #1abc9c,stop:1 #16a085);"
            "color:white;font-size:14px;font-weight:bold;"
            "padding:8px;border-radius:6px;"
        )
        layout.addWidget(self.runBtn)

        # progress bar + log
        self.progressBar = qt.QProgressBar()
        self.progressBar.setVisible(False)
        layout.addWidget(self.progressBar)

        self.logBox = qt.QPlainTextEdit()
        self.logBox.setReadOnly(True)
        self.logBox.setMaximumHeight(120)
        self.logBox.setStyleSheet("font-family:monospace;font-size:11px;")
        layout.addWidget(self.logBox)

        layout.addStretch(1)

    def _iconBtn(self, label, tip, color="#555"):
        btn = qt.QPushButton(label)
        btn.setToolTip(tip)
        btn.setStyleSheet(
            f"background:{color};color:white;font-size:11px;font-weight:bold;"
            "border-radius:4px;padding:3px 8px;"
        )
        return btn

    # --------------------------------------------------------- Signals ------

    def _connectSignals(self):
        self.newWorkflowBtn.clicked.connect(self._onNewWorkflow)
        self.saveWorkflowBtn.clicked.connect(self._onSaveWorkflow)
        self.delWorkflowBtn.clicked.connect(self._onDeleteWorkflow)
        self.importBtn.clicked.connect(self._onImport)
        self.exportBtn.clicked.connect(self._onExport)
        self.addStepBtn.clicked.connect(self._onAddStep)
        self.removeStepBtn.clicked.connect(self._onRemoveStep)
        self.moveUpBtn.clicked.connect(self._onMoveUp)
        self.moveDownBtn.clicked.connect(self._onMoveDown)
        self.stepList.currentRowChanged.connect(self._onStepSelected)
        self.singleRadio.toggled.connect(self._onModeToggled)
        self.batchRadio.toggled.connect(self._onModeToggled)
        self.browseFolderBtn.clicked.connect(self._onBrowseFolder)
        self.browseOutputBtn.clicked.connect(self._onBrowseOutput)
        self.runBtn.clicked.connect(self._onRun)
        self.workflowCombo.currentIndexChanged.connect(self._onWorkflowSelected)

    # --------------------------------------------------------- Step list ----

    def _onAddStep(self):
        stepName = self.stepTypeCombo.currentText
        idx = self.stepList.count
        label = f"{idx+1}. {stepName}"
        item = qt.QListWidgetItem(label)
        # store default params in item's UserRole
        defaults = {
            k: v["default"]
            for k, v in self.STEP_CATALOGUE[stepName]["params"].items()
        }
        item.setData(qt.Qt.UserRole, {"step": stepName, "params": defaults})
        self.stepList.addItem(item)
        self.stepList.setCurrentRow(self.stepList.count - 1)
        self._renumberSteps()

    def _onRemoveStep(self):
        row = self.stepList.currentRow
        if row >= 0:
            self.stepList.takeItem(row)
            self._renumberSteps()
            self._clearParamPanel()

    def _onMoveUp(self):
        row = self.stepList.currentRow
        if row > 0:
            item = self.stepList.takeItem(row)
            self.stepList.insertItem(row - 1, item)
            self.stepList.setCurrentRow(row - 1)
            self._renumberSteps()

    def _onMoveDown(self):
        row = self.stepList.currentRow
        if row < self.stepList.count - 1:
            item = self.stepList.takeItem(row)
            self.stepList.insertItem(row + 1, item)
            self.stepList.setCurrentRow(row + 1)
            self._renumberSteps()

    def _renumberSteps(self):
        for i in range(self.stepList.count):
            item = self.stepList.item(i)
            data = item.data(qt.Qt.UserRole)
            item.setText(f"{i+1}. {data['step']}")

    # --------------------------------------------------------- Param panel --

    def _onStepSelected(self, row):
        self._clearParamPanel()
        if row < 0:
            return
        item = self.stepList.item(row)
        if item is None:
            return
        data   = item.data(qt.Qt.UserRole)
        stepName = data["step"]
        params   = data["params"]
        catalogue = self.STEP_CATALOGUE[stepName]["params"]

        self.paramGroupBox.text = f"Parameters – {stepName}"

        for pname, pdef in catalogue.items():
            ptype = pdef["type"]
            val   = params.get(pname, pdef["default"])
            widget = None

            if ptype == "float":
                w = qt.QDoubleSpinBox()
                w.setMinimum(pdef.get("min", -1e9))
                w.setMaximum(pdef.get("max",  1e9))
                w.setSingleStep(0.01)
                w.setDecimals(4)
                w.setValue(float(val))
                widget = w
            elif ptype == "int":
                w = qt.QSpinBox()
                w.setMinimum(pdef.get("min", 0))
                w.setMaximum(pdef.get("max", 100000))
                w.setValue(int(val))
                widget = w
            elif ptype == "bool":
                w = qt.QCheckBox()
                w.setChecked(bool(val))
                widget = w
            elif ptype == "enum":
                w = qt.QComboBox()
                for v in pdef["values"]:
                    w.addItem(v)
                if val in pdef["values"]:
                    w.setCurrentIndex(pdef["values"].index(val))
                widget = w
            elif ptype == "string":
                w = qt.QLineEdit(str(val))
                widget = w
            elif ptype == "path":
                w = qt.QWidget()
                hl = qt.QHBoxLayout(w)
                hl.setContentsMargins(0, 0, 0, 0)
                le = qt.QLineEdit(str(val))
                le.setObjectName("pathEdit")
                btn = qt.QPushButton("…")
                btn.setFixedWidth(28)
                btn.clicked.connect(lambda _, le=le: le.setText(
                    qt.QFileDialog.getExistingDirectory(None, "Select folder")))
                hl.addWidget(le)
                hl.addWidget(btn)
                widget = w
            elif ptype == "node":
                w = slicer.qMRMLNodeComboBox()
                w.nodeTypes = [pdef.get("nodeType", "vtkMRMLScalarVolumeNode")]
                w.noneEnabled = True
                w.addEnabled = False
                w.removeEnabled = False
                w.setMRMLScene(slicer.mrmlScene)
                widget = w
            elif ptype == "text":
                w = qt.QPlainTextEdit()
                w.setPlainText(str(val))
                w.setMinimumHeight(80)
                widget = w

            if widget:
                self.paramLayout.addRow(pname + ":", widget)
                self._paramWidgets[pname] = (ptype, widget)
                # live save back to item on any change
                self._connectParamSave(ptype, widget, item, pname)

    def _connectParamSave(self, ptype, widget, item, pname):
        """Whenever a param widget changes, persist value back to list item."""
        def _save():
            data = item.data(qt.Qt.UserRole)
            if ptype == "float":
                data["params"][pname] = widget.value
            elif ptype == "int":
                data["params"][pname] = widget.value
            elif ptype == "bool":
                data["params"][pname] = widget.isChecked()
            elif ptype == "enum":
                data["params"][pname] = widget.currentText
            elif ptype == "string":
                data["params"][pname] = widget.text
            elif ptype == "path":
                le = widget.findChild(qt.QLineEdit, "pathEdit")
                if le:
                    data["params"][pname] = le.text
            elif ptype == "node":
                node = widget.currentNode()
                data["params"][pname] = node.GetID() if node else ""
            elif ptype == "text":
                data["params"][pname] = widget.toPlainText()
            item.setData(qt.Qt.UserRole, data)

        if ptype in ("float", "int"):
            widget.valueChanged.connect(lambda _: _save())
        elif ptype == "bool":
            widget.stateChanged.connect(lambda _: _save())
        elif ptype == "enum":
            widget.currentIndexChanged.connect(lambda _: _save())
        elif ptype in ("string", "path"):
            if ptype == "string":
                widget.textChanged.connect(lambda _: _save())
            else:
                le = widget.findChild(qt.QLineEdit, "pathEdit")
                if le:
                    le.textChanged.connect(lambda _: _save())
        elif ptype == "node":
            widget.currentNodeChanged.connect(lambda _: _save())
        elif ptype == "text":
            widget.textChanged.connect(_save)

    def _clearParamPanel(self):
        while self.paramLayout.rowCount() > 0:
            self.paramLayout.removeRow(0)
        self._paramWidgets.clear()

    # --------------------------------------------------------- Input mode ---

    def _onModeToggled(self):
        batch = self.batchRadio.isChecked()
        self.singleNodeSelector.setVisible(not batch)
        self.batchWidget.setVisible(batch)
        self.batchFilterWidget.setVisible(batch)

    def _onBrowseFolder(self):
        d = qt.QFileDialog.getExistingDirectory(None, "Select input folder")
        if d:
            self.folderEdit.setText(d)

    def _onBrowseOutput(self):
        d = qt.QFileDialog.getExistingDirectory(None, "Select output folder")
        if d:
            self.outputDirEdit.setText(d)

    # --------------------------------------------------------- Workflow IO --

    def _collectSteps(self):
        steps = []
        for i in range(self.stepList.count):
            item = self.stepList.item(i)
            steps.append(item.data(qt.Qt.UserRole))
        return steps

    def _loadSteps(self, steps):
        self.stepList.clear()
        self._clearParamPanel()
        for i, s in enumerate(steps):
            item = qt.QListWidgetItem(f"{i+1}. {s['step']}")
            item.setData(qt.Qt.UserRole, s)
            self.stepList.addItem(item)

    def _refreshWorkflowList(self):
        current = self.workflowCombo.currentText
        self.workflowCombo.blockSignals(True)
        self.workflowCombo.clear()
        for name in self.logic.listWorkflows():
            self.workflowCombo.addItem(name)
        idx = self.workflowCombo.findText(current)
        self.workflowCombo.setCurrentIndex(max(idx, 0))
        self.workflowCombo.blockSignals(False)
        # Manually trigger a load for whichever item ended up selected
        self._onWorkflowSelected(self.workflowCombo.currentIndex)

    def _onWorkflowSelected(self, idx):
        name = self.workflowCombo.currentText
        if not name:
            return
        wf = self.logic.loadWorkflow(name)
        if wf:
            self._loadSteps(wf.get("steps", []))

    def _askName(self, title):
        """Show a simple name-input dialog. Returns the entered string or None."""
        dialog = qt.QDialog()
        dialog.setWindowTitle(title)
        dialog.setMinimumWidth(300)
        vbox = qt.QVBoxLayout(dialog)
        vbox.addWidget(qt.QLabel("Workflow name:"))
        lineEdit = qt.QLineEdit()
        vbox.addWidget(lineEdit)
        # Plain QPushButtons — more reliable than QDialogButtonBox in Slicer
        btnRow = qt.QHBoxLayout()
        okBtn     = qt.QPushButton("OK")
        cancelBtn = qt.QPushButton("Cancel")
        btnRow.addWidget(okBtn)
        btnRow.addWidget(cancelBtn)
        vbox.addLayout(btnRow)
        okBtn.clicked.connect(lambda: dialog.accept())
        cancelBtn.clicked.connect(lambda: dialog.reject())
        # Also accept on Enter key in the line edit
        lineEdit.returnPressed.connect(dialog.accept)
        result = dialog.exec_()
        if result == 1:          # 1 == QDialog::Accepted in PythonQt
            return lineEdit.text.strip() or None
        return None

    def _onNewWorkflow(self):
        name = self._askName("New Workflow")
        if name:
            self.logic.saveWorkflow(name, {"steps": []})
            self._refreshWorkflowList()
            idx = self.workflowCombo.findText(name)
            if idx >= 0:
                self.workflowCombo.setCurrentIndex(idx)

    def _onSaveWorkflow(self):
        name = self.workflowCombo.currentText
        if not name:
            name = self._askName("Save Workflow")
            if not name:
                return
        self.logic.saveWorkflow(name, {"steps": self._collectSteps()})
        self._refreshWorkflowList()
        self._log(f"Workflow '{name}' saved.")

    def _onDeleteWorkflow(self):
        name = self.workflowCombo.currentText
        if not name:
            return
        if qt.QMessageBox.question(
            None, "Delete", f"Delete workflow '{name}'?",
            qt.QMessageBox.Yes | qt.QMessageBox.No
        ) == qt.QMessageBox.Yes:
            self.logic.deleteWorkflow(name)
            self._refreshWorkflowList()

    def _onImport(self):
        path, _ = qt.QFileDialog.getOpenFileName(
            None, "Import Workflow", "", "JSON (*.json)")
        if not path:
            return
        try:
            with open(path) as f:
                data = json.load(f)
            name = os.path.splitext(os.path.basename(path))[0]
            self.logic.saveWorkflow(name, data)
            self._refreshWorkflowList()
            self._log(f"Imported '{name}' from {path}")
        except Exception as e:
            self._log(f"Import error: {e}", error=True)

    def _onExport(self):
        name = self.workflowCombo.currentText
        if not name:
            self._log("No workflow selected.", error=True)
            return
        path, _ = qt.QFileDialog.getSaveFileName(
            None, "Export Workflow", f"{name}.json", "JSON (*.json)")
        if not path:
            return
        wf = {"steps": self._collectSteps()}
        with open(path, "w") as f:
            json.dump(wf, f, indent=2)
        self._log(f"Exported to {path}")

    # --------------------------------------------------------- Run ----------

    def _onRun(self):
        steps = self._collectSteps()
        if not steps:
            self._log("No steps defined.", error=True)
            return

        outputDir = self.outputDirEdit.text.strip() or None
        keepInter = self.keepIntermediateCheck.isChecked()

        if self.batchRadio.isChecked():
            folder = self.folderEdit.text.strip()
            if not folder or not os.path.isdir(folder):
                self._log("Invalid input folder.", error=True)
                return
            fileFilters = self.fileFilterEdit.text.split()
            files = []
            for root, _, fnames in os.walk(folder):
                for fn in fnames:
                    if any(fn.endswith(ext.lstrip("*")) for ext in fileFilters):
                        files.append(os.path.join(root, fn))
            if not files:
                self._log("No matching files found.", error=True)
                return
            self.progressBar.setVisible(True)
            self.progressBar.setMaximum(len(files))
            for i, fpath in enumerate(files):
                self._log(f"[{i+1}/{len(files)}] {os.path.basename(fpath)}")
                self.logic.runWorkflowOnFile(
                    fpath, steps, outputDir, keepInter, self._log)
                self.progressBar.setValue(i + 1)
                slicer.app.processEvents()
            self.progressBar.setVisible(False)
            self._log("Batch complete ✓")
        else:
            node = self.singleNodeSelector.currentNode()
            if node is None:
                self._log("No volume node selected.", error=True)
                return
            self.progressBar.setVisible(True)
            self.progressBar.setMaximum(len(steps))
            self.logic.runWorkflowOnNode(
                node, steps, outputDir, keepInter,
                self._log, self.progressBar)
            self.progressBar.setVisible(False)
            self._log("Workflow complete ✓")

    # --------------------------------------------------------- Log ----------

    def _log(self, msg, error=False):
        prefix = "❌ " if error else "✔ "
        self.logBox.appendPlainText(prefix + msg)
        logging.warning(msg) if error else logging.info(msg)
        slicer.app.processEvents()


# ---------------------------------------------------------------------------
# Logic
# ---------------------------------------------------------------------------

class WorkflowModuleLogic(ScriptedLoadableModuleLogic):

    @property
    def WORKFLOWS_DIR(self):
        # slicer.app.slicerHome/../<settings-dir>/SlicerWorkflows  is fragile.
        # Use the directory that contains the INI file, resolved properly.
        ini = slicer.app.userSettings().fileName()          # e.g. /home/user/.config/slicer.org/Slicer.ini
        base = os.path.dirname(ini)                         # /home/user/.config/slicer.org/
        d = os.path.join(base, "SlicerWorkflows")
        os.makedirs(d, exist_ok=True)
        return d

    def __init__(self):
        super().__init__()

    def _wfPath(self, name):
        return os.path.join(self.WORKFLOWS_DIR, name + ".json")

    def listWorkflows(self):
        return [
            os.path.splitext(f)[0]
            for f in sorted(os.listdir(self.WORKFLOWS_DIR))
            if f.endswith(".json")
        ]

    def saveWorkflow(self, name, data):
        with open(self._wfPath(name), "w") as f:
            json.dump(data, f, indent=2)

    def loadWorkflow(self, name):
        p = self._wfPath(name)
        if os.path.exists(p):
            with open(p) as f:
                return json.load(f)
        return None

    def deleteWorkflow(self, name):
        p = self._wfPath(name)
        if os.path.exists(p):
            os.remove(p)

    # ---------------------------------------------------------------- run --

    def runWorkflowOnFile(self, filePath, steps, outputDir, keepInter, logFn):
        """Load a file into scene, run pipeline, optionally save, then clean up."""
        loadedNodes = slicer.util.loadVolume(filePath, returnNode=True)
        if isinstance(loadedNodes, tuple):
            volumeNode = loadedNodes[1]
        else:
            volumeNode = loadedNodes
        if volumeNode is None:
            logFn(f"Could not load {filePath}", error=True)
            return
        self.runWorkflowOnNode(volumeNode, steps, outputDir, keepInter, logFn)
        if not keepInter:
            slicer.mrmlScene.RemoveNode(volumeNode)

    def runWorkflowOnNode(self, volumeNode, steps, outputDir,
                          keepInter, logFn, progressBar=None):
        """Execute each step in order on volumeNode."""
        currentVolume = volumeNode
        currentSeg    = None

        for i, stepData in enumerate(steps):
            stepName = stepData["step"]
            params   = stepData["params"]
            logFn(f"  Step {i+1}: {stepName}")
            try:
                result = self._executeStep(
                    stepName, params, currentVolume, currentSeg,
                    outputDir, logFn)
                # steps may return updated (volume, segmentation) tuple
                if isinstance(result, tuple):
                    if result[0] is not None:
                        currentVolume = result[0]
                    if result[1] is not None:
                        currentSeg = result[1]
                elif result is not None:
                    if hasattr(result, "IsA"):
                        if result.IsA("vtkMRMLScalarVolumeNode"):
                            currentVolume = result
                        elif result.IsA("vtkMRMLSegmentationNode"):
                            currentSeg = result
            except Exception as e:
                logFn(f"  ✗ {stepName} failed: {e}", error=True)
                logging.exception(e)
            if progressBar:
                progressBar.setValue(i + 1)
                slicer.app.processEvents()

    def _executeStep(self, stepName, params, volumeNode, segNode, outputDir, logFn):
        """Dispatch to the appropriate step implementation."""
        import SimpleITK as sitk
        import sitkUtils

        if stepName == "Gaussian Smoothing":
            sigma = float(params.get("sigma", 1.0))
            img   = sitkUtils.PullVolumeFromSlicer(volumeNode)
            img   = sitk.SmoothingRecursiveGaussian(img, sigma)
            outNode = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLScalarVolumeNode",
                volumeNode.GetName() + "_Gauss")
            sitkUtils.PushVolumeToSlicer(img, outNode)
            logFn(f"    σ={sigma}")
            return outNode, segNode

        elif stepName == "Median Filtering":
            radius = int(params.get("radius", 1))
            img    = sitkUtils.PullVolumeFromSlicer(volumeNode)
            img    = sitk.Median(img, [radius]*img.GetDimension())
            outNode = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLScalarVolumeNode",
                volumeNode.GetName() + "_Median")
            sitkUtils.PushVolumeToSlicer(img, outNode)
            return outNode, segNode

        elif stepName == "N4 Bias Correction":
            nIter  = int(params.get("numberOfIterations", 50))
            shrink = int(params.get("shrinkFactor", 4))
            thresh = float(params.get("convergenceThreshold", 0.001))
            img    = sitkUtils.PullVolumeFromSlicer(volumeNode)
            img    = sitk.Cast(img, sitk.sitkFloat32)
            shrunk = sitk.Shrink(img, [shrink]*img.GetDimension())
            mask   = sitk.OtsuThreshold(shrunk, 0, 1, 200)
            corrector = sitk.N4BiasFieldCorrectionImageFilter()
            corrector.SetMaximumNumberOfIterations([nIter]*4)
            corrector.SetConvergenceThreshold(thresh)
            corrector.Execute(shrunk, mask)
            log_bf = corrector.GetLogBiasFieldAsImage(img)
            corrected = img / sitk.Exp(log_bf)
            outNode = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLScalarVolumeNode",
                volumeNode.GetName() + "_N4")
            sitkUtils.PushVolumeToSlicer(corrected, outNode)
            return outNode, segNode

        elif stepName == "Registration":
            fixedID = params.get("fixedVolumeNode", "")
            fixedNode = (slicer.mrmlScene.GetNodeByID(fixedID)
                         if fixedID else None)
            if fixedNode is None:
                logFn("    No fixed volume – skipping registration", error=True)
                return volumeNode, segNode
            tType = params.get("transformType", "Rigid")
            sampling = float(params.get("samplingPercentage", 0.25))
            nIter    = int(params.get("numberOfIterations", 200))
            paramNode = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLBRAINSFitParametersNode")
            paramNode.SetFixedVolume(fixedNode.GetID())
            paramNode.SetMovingVolume(volumeNode.GetID())
            paramNode.SetTransformType(tType)
            paramNode.SetSamplingPercentage(sampling)
            paramNode.SetNumberOfIterations(nIter)
            outputTransformNode = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLLinearTransformNode",
                volumeNode.GetName() + "_reg_transform")
            paramNode.SetOutputTransform(outputTransformNode.GetID())
            slicer.cli.runSync(slicer.modules.brainsfit, paramNode)
            volumeNode.SetAndObserveTransformNodeID(
                outputTransformNode.GetID())
            slicer.vtkSlicerTransformLogic().hardenTransform(volumeNode)
            slicer.mrmlScene.RemoveNode(paramNode)
            return volumeNode, segNode

        elif stepName == "Threshold Segmentation":
            lo  = float(params.get("lowerThreshold", 100))
            hi  = float(params.get("upperThreshold", 3000))
            seg = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLSegmentationNode",
                volumeNode.GetName() + "_ThreshSeg")
            seg.CreateDefaultDisplayNodes()
            seg.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
            labelmap = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLLabelMapVolumeNode")
            import sitkUtils
            img = sitkUtils.PullVolumeFromSlicer(volumeNode)
            binaryImg = sitk.BinaryThreshold(
                img, lowerThreshold=lo, upperThreshold=hi,
                insideValue=1, outsideValue=0)
            sitkUtils.PushVolumeToSlicer(binaryImg, labelmap)
            slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(
                labelmap, seg)
            slicer.mrmlScene.RemoveNode(labelmap)
            return volumeNode, seg

        elif stepName == "Otsu Segmentation":
            inside  = int(params.get("insideValue", 1))
            outside = int(params.get("outsideValue", 0))
            img = sitkUtils.PullVolumeFromSlicer(volumeNode)
            binary = sitk.OtsuThreshold(img, outside, inside)
            seg = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLSegmentationNode",
                volumeNode.GetName() + "_OtsuSeg")
            seg.CreateDefaultDisplayNodes()
            seg.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
            labelmap = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLLabelMapVolumeNode")
            sitkUtils.PushVolumeToSlicer(binary, labelmap)
            slicer.modules.segmentations.logic().ImportLabelmapToSegmentationNode(
                labelmap, seg)
            slicer.mrmlScene.RemoveNode(labelmap)
            return volumeNode, seg

        elif stepName == "Volume Statistics":
            if segNode is None:
                logFn("    No segmentation node – skipping", error=True)
                return volumeNode, segNode
            import SegmentStatistics
            ssLogic = SegmentStatistics.SegmentStatisticsLogic()
            ssLogic.getParameterNode().SetParameter(
                "Segmentation", segNode.GetID())
            ssLogic.getParameterNode().SetParameter(
                "ScalarVolume", volumeNode.GetID())
            ssLogic.computeStatistics()
            if params.get("saveCSV", True):
                od = params.get("outputDir", "") or outputDir or ""
                if od:
                    csvPath = os.path.join(
                        od, volumeNode.GetName() + "_volumeStats.csv")
                    ssLogic.exportToCSVFile(csvPath)
                    logFn(f"    Saved {csvPath}")
            return volumeNode, segNode

        elif stepName == "Intensity Statistics":
            if segNode is None:
                logFn("    No segmentation – skipping", error=True)
                return volumeNode, segNode
            import SegmentStatistics
            ssLogic = SegmentStatistics.SegmentStatisticsLogic()
            ssLogic.getParameterNode().SetParameter(
                "Segmentation", segNode.GetID())
            ssLogic.getParameterNode().SetParameter(
                "ScalarVolume", volumeNode.GetID())
            ssLogic.computeStatistics()
            if params.get("saveCSV", True):
                od = params.get("outputDir", "") or outputDir or ""
                if od:
                    csvPath = os.path.join(
                        od, volumeNode.GetName() + "_intensityStats.csv")
                    ssLogic.exportToCSVFile(csvPath)
                    logFn(f"    Saved {csvPath}")
            return volumeNode, segNode

        elif stepName == "Save Volume":
            od = params.get("outputDir", "") or outputDir or ""
            fmt = params.get("fileFormat", "nrrd")
            if od:
                fname = volumeNode.GetName() + "." + fmt
                fpath = os.path.join(od, fname)
                slicer.util.saveNode(volumeNode, fpath)
                logFn(f"    Saved {fpath}")
            else:
                logFn("    No output dir – skipping Save Volume", error=True)
            return volumeNode, segNode

        elif stepName == "Custom Python Script":
            script = params.get("script", "")
            exec(script, {"volumeNode": volumeNode,
                          "segmentationNode": segNode,
                          "slicer": slicer,
                          "logFn": logFn})
            return volumeNode, segNode

        else:
            logFn(f"    Unknown step '{stepName}' – skipped", error=True)
            return volumeNode, segNode


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class WorkflowModuleTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear(0)

    def runTest(self):
        self.setUp()
        self.test_GaussianSmoothingStep()

    def test_GaussianSmoothingStep(self):
        self.delayDisplay("Creating test volume …")
        import numpy as np
        voxels = np.random.randint(0, 200, (64, 64, 64), dtype=np.int16)
        volumeNode = slicer.util.addVolumeFromArray(voxels)
        logic = WorkflowModuleLogic()
        steps = [{"step": "Gaussian Smoothing", "params": {"sigma": 1.5}}]
        logic.runWorkflowOnNode(volumeNode, steps, None, True, print)
        self.delayDisplay("Test passed ✓")
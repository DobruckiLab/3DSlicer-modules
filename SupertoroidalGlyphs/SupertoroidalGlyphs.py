import os
import csv
import logging

import numpy as np
import vtk
import qt
import ctk
import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin


PARAM_CHOICES = ["(none)", "color", "N1", "N2", "size"]
PARAM_TARGET_RANGE = {
    "color": None,
    "N1": (0.3, 4.0),
    "N2": (0.3, 4.0),
    "size": (0.4, 1.0),
}
PARAM_FALLBACK = {
    "N1": (0.3 + 4.0) / 2.0,
    "N2": (0.3 + 4.0) / 2.0,
    "size": (0.4 + 1.0) / 2.0,
}


class SupertoroidalGlyphs(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Supertoroidal Glyphs"
        self.parent.categories = ["Informatics"]
        self.parent.dependencies = []
        self.parent.contributors = ["Generated with Claude"]
        self.parent.helpText = """
Generate a supertoroidal-glyph field from co-registered imaging data and
display it superimposed on a CT volume. Select a CT plus up to four
co-registered volumes (e.g. PET), assign each to a glyph parameter (color, N1,
N2, or size), set a binning factor (block of NxNxN voxels averaged into one
glyph), and generate. Export the glyph table as CSV using the standalone
format: label, R_mm, A_mm, S_mm, param1_N1, param2_N2, param3_size, param4_color
"""
        self.parent.acknowledgementText = ""


class SupertoroidalGlyphsWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self.petSelectors = []
        self.paramCombos = []
        self.lastRows = None

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = SupertoroidalGlyphsLogic()

        volCollapsible = ctk.ctkCollapsibleButton()
        volCollapsible.text = "Volumes"
        self.layout.addWidget(volCollapsible)
        volForm = qt.QFormLayout(volCollapsible)

        self.ctSelector = self._makeVolumeSelector(
            "CT volume shown in 3D; glyphs are superimposed on it.")
        volForm.addRow("CT volume:", self.ctSelector)

        self.refSelector = self._makeVolumeSelector(
            "Volume whose voxel grid defines the glyph sampling grid.",
            noneEnabled=True)
        volForm.addRow("Reference grid volume:", self.refSelector)

        for idx in range(4):
            row = qt.QWidget()
            rowLayout = qt.QHBoxLayout(row)
            rowLayout.setContentsMargins(0, 0, 0, 0)
            sel = self._makeVolumeSelector(
                "Co-registered data volume %d (e.g. PET)." % (idx + 1),
                noneEnabled=True)
            combo = qt.QComboBox()
            combo.addItems(PARAM_CHOICES)
            combo.setToolTip("Which glyph parameter this volume drives.")
            combo.setCurrentIndex(min(idx + 1, len(PARAM_CHOICES) - 1))
            rowLayout.addWidget(sel, 2)
            rowLayout.addWidget(qt.QLabel(" encodes "))
            rowLayout.addWidget(combo, 1)
            self.petSelectors.append(sel)
            self.paramCombos.append(combo)
            volForm.addRow("Data volume %d:" % (idx + 1), row)

        binCollapsible = ctk.ctkCollapsibleButton()
        binCollapsible.text = "Binning"
        self.layout.addWidget(binCollapsible)
        binForm = qt.QFormLayout(binCollapsible)

        self.binSizeSlider = ctk.ctkSliderWidget()
        self.binSizeSlider.minimum = 1
        self.binSizeSlider.maximum = 64
        self.binSizeSlider.singleStep = 1
        self.binSizeSlider.decimals = 0
        self.binSizeSlider.value = 10
        self.binSizeSlider.suffix = " voxels"
        self.binSizeSlider.setToolTip(
            "Block edge length: this many voxels per axis are averaged into one "
            "glyph (e.g. 10 => 10x10x10 block).")
        binForm.addRow("Bin size:", self.binSizeSlider)

        self.aggCombo = qt.QComboBox()
        self.aggCombo.addItems(["mean (nonzero)", "mean (all)", "max"])
        binForm.addRow("Aggregation:", self.aggCombo)

        self.skipBackgroundCheck = qt.QCheckBox()
        self.skipBackgroundCheck.checked = True
        binForm.addRow("Skip background bins:", self.skipBackgroundCheck)

        self.bgThresholdSlider = ctk.ctkSliderWidget()
        self.bgThresholdSlider.minimum = 0
        self.bgThresholdSlider.maximum = 1000
        self.bgThresholdSlider.singleStep = 1
        self.bgThresholdSlider.value = 0
        self.bgThresholdSlider.setToolTip(
            "A reference voxel counts as foreground if its value exceeds this.")
        binForm.addRow("Background threshold:", self.bgThresholdSlider)

        optCollapsible = ctk.ctkCollapsibleButton()
        optCollapsible.text = "Display options"
        self.layout.addWidget(optCollapsible)
        optForm = qt.QFormLayout(optCollapsible)

        self.colorSelector = slicer.qMRMLNodeComboBox()
        self.colorSelector.nodeTypes = ["vtkMRMLColorTableNode", "vtkMRMLProceduralColorNode"]
        self.colorSelector.addEnabled = False
        self.colorSelector.removeEnabled = False
        self.colorSelector.noneEnabled = False
        self.colorSelector.showHidden = True
        self.colorSelector.setMRMLScene(slicer.mrmlScene)
        optForm.addRow("Color map:", self.colorSelector)
        self._selectColorNodeByName("Viridis")

        self.ringRadiusSlider = ctk.ctkSliderWidget()
        self.ringRadiusSlider.minimum = 0.5
        self.ringRadiusSlider.maximum = 30.0
        self.ringRadiusSlider.singleStep = 0.5
        self.ringRadiusSlider.value = 6.0
        self.ringRadiusSlider.suffix = " mm"
        optForm.addRow("Base ring radius:", self.ringRadiusSlider)

        self.crossFractionSlider = ctk.ctkSliderWidget()
        self.crossFractionSlider.minimum = 0.05
        self.crossFractionSlider.maximum = 1.0
        self.crossFractionSlider.singleStep = 0.05
        self.crossFractionSlider.value = 0.5
        optForm.addRow("Cross-section fraction:", self.crossFractionSlider)

        self.resolutionSlider = ctk.ctkSliderWidget()
        self.resolutionSlider.minimum = 8
        self.resolutionSlider.maximum = 80
        self.resolutionSlider.singleStep = 4
        self.resolutionSlider.decimals = 0
        self.resolutionSlider.value = 32
        self.resolutionSlider.setToolTip(
            "Surface tessellation per glyph. Lower if you have thousands of glyphs.")
        optForm.addRow("Resolution:", self.resolutionSlider)

        self.volumeRenderCheck = qt.QCheckBox()
        self.volumeRenderCheck.checked = True
        optForm.addRow("Show CT (volume render):", self.volumeRenderCheck)

        self.legendCheck = qt.QCheckBox()
        self.legendCheck.checked = True
        optForm.addRow("Show color legend:", self.legendCheck)

        self.applyButton = qt.QPushButton("Generate glyphs")
        self.applyButton.enabled = False
        self.layout.addWidget(self.applyButton)

        self.exportButton = qt.QPushButton("Export glyph table (CSV)...")
        self.exportButton.enabled = False
        self.layout.addWidget(self.exportButton)

        self.clearButton = qt.QPushButton("Clear glyphs")
        self.layout.addWidget(self.clearButton)

        self.statusLabel = qt.QLabel("")
        self.statusLabel.wordWrap = True
        self.layout.addWidget(self.statusLabel)
        self.layout.addStretch(1)

        self.applyButton.connect("clicked(bool)", self.onApply)
        self.exportButton.connect("clicked(bool)", self.onExport)
        self.clearButton.connect("clicked(bool)", self.onClear)
        self.ctSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateApplyState)
        for sel in self.petSelectors:
            sel.connect("currentNodeChanged(vtkMRMLNode*)", self.updateApplyState)
        self.updateApplyState()

    def _makeVolumeSelector(self, tooltip, noneEnabled=False):
        sel = slicer.qMRMLNodeComboBox()
        sel.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        sel.addEnabled = False
        sel.removeEnabled = False
        sel.noneEnabled = noneEnabled
        sel.showHidden = False
        sel.setMRMLScene(slicer.mrmlScene)
        sel.setToolTip(tooltip)
        return sel

    def _selectColorNodeByName(self, name):
        node = slicer.util.getFirstNodeByName(name, "vtkMRMLColorTableNode")
        if node is not None:
            self.colorSelector.setCurrentNode(node)

    def _assignments(self):
        out = []
        used = set()
        for sel, combo in zip(self.petSelectors, self.paramCombos):
            vol = sel.currentNode()
            param = combo.currentText
            if vol is None or param == "(none)":
                continue
            if param in used:
                raise ValueError(
                    "Parameter '%s' is assigned to more than one volume." % param)
            used.add(param)
            out.append((vol, param))
        return out

    def updateApplyState(self, *args):
        hasCT = self.ctSelector.currentNode() is not None
        hasData = any(sel.currentNode() is not None for sel in self.petSelectors)
        self.applyButton.enabled = hasCT and hasData

    def onApply(self):
        try:
            self.statusLabel.text = "Generating (sampling volumes and building glyphs)..."
            slicer.app.processEvents()
            assignments = self._assignments()
            aggMap = {"mean (nonzero)": "mean_nonzero",
                      "mean (all)": "mean_all", "max": "max"}
            refNode = self.refSelector.currentNode()
            if refNode is None and assignments:
                refNode = assignments[0][0]
            rows, n = self.logic.run(
                ctNode=self.ctSelector.currentNode(),
                referenceNode=refNode,
                assignments=assignments,
                binSize=int(self.binSizeSlider.value),
                aggregation=aggMap[self.aggCombo.currentText],
                skipBackground=bool(self.skipBackgroundCheck.checked),
                backgroundThreshold=float(self.bgThresholdSlider.value),
                colorNode=self.colorSelector.currentNode(),
                baseRingRadius=float(self.ringRadiusSlider.value),
                crossFraction=float(self.crossFractionSlider.value),
                resolution=int(self.resolutionSlider.value),
                showVolumeRendering=bool(self.volumeRenderCheck.checked),
                showLegend=bool(self.legendCheck.checked),
            )
            self.lastRows = rows
            self.exportButton.enabled = bool(rows)
            self.statusLabel.text = (
                "Done. Created %d glyphs from %d assigned volume(s)."
                % (n, len(assignments)))
        except Exception as e:
            slicer.util.errorDisplay("Failed: %s" % str(e))
            self.statusLabel.text = "Error: %s" % str(e)
            logging.exception("SupertoroidalGlyphs generation failed")

    def onExport(self):
        if not self.lastRows:
            slicer.util.errorDisplay("No glyph table to export. Generate glyphs first.")
            return
        path = qt.QFileDialog.getSaveFileName(
            None, "Save glyph table", "glyphs.csv", "CSV files (*.csv)")
        if not path:
            return
        try:
            self.logic.writeTable(self.lastRows, path)
            self.statusLabel.text = "Saved %d rows to %s" % (len(self.lastRows), path)
        except Exception as e:
            slicer.util.errorDisplay("Failed to save CSV: %s" % str(e))

    def onClear(self):
        removed = self.logic.clear()
        self.lastRows = None
        self.exportButton.enabled = False
        self.statusLabel.text = "Removed %d glyph nodes." % removed


class SupertoroidalGlyphsLogic(ScriptedLoadableModuleLogic):

    FOLDER_NAME = "Supertoroidal Glyphs"
    COLOR_ARRAY_NAME = "glyph_color"

    @staticmethod
    def _ijkToRas(volumeNode):
        m = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(m)
        return np.array([[m.GetElement(r, c) for c in range(4)] for r in range(4)])

    @staticmethod
    def _aggregate(block, mode):
        if block.size == 0:
            return None
        if mode == "max":
            return float(block.max())
        if mode == "mean_all":
            return float(block.mean())
        nz = block[block != 0]
        return float(nz.mean()) if nz.size else None

    def _sampleAtRas(self, arr, invAffine, rasCenter, half, mode):
        ijk = invAffine @ np.array([rasCenter[0], rasCenter[1], rasCenter[2], 1.0])
        ci, cj, ck = ijk[0], ijk[1], ijk[2]
        K, J, I = arr.shape
        i0, i1 = int(round(ci - half)), int(round(ci + half))
        j0, j1 = int(round(cj - half)), int(round(cj + half))
        k0, k1 = int(round(ck - half)), int(round(ck + half))
        i0, j0, k0 = max(i0, 0), max(j0, 0), max(k0, 0)
        i1, j1, k1 = min(i1, I - 1), min(j1, J - 1), min(k1, K - 1)
        if i0 > i1 or j0 > j1 or k0 > k1:
            return None
        block = arr[k0:k1 + 1, j0:j1 + 1, i0:i1 + 1]
        return self._aggregate(block, mode)

    def sampleGlyphs(self, referenceNode, assignments, binSize, aggregation,
                     skipBackground, backgroundThreshold):
        refArr = slicer.util.arrayFromVolume(referenceNode)
        refAffine = self._ijkToRas(referenceNode)
        Kr, Jr, Ir = refArr.shape
        half = binSize / 2.0
        volInfo = []
        for vol, param in assignments:
            arr = slicer.util.arrayFromVolume(vol)
            inv = np.linalg.inv(self._ijkToRas(vol))
            volInfo.append((param, arr, inv))
        rows = []
        gid = 0
        for k0 in range(0, Kr, binSize):
            for j0 in range(0, Jr, binSize):
                for i0 in range(0, Ir, binSize):
                    block = refArr[k0:min(k0 + binSize, Kr),
                                   j0:min(j0 + binSize, Jr),
                                   i0:min(i0 + binSize, Ir)]
                    if skipBackground and not np.any(block > backgroundThreshold):
                        continue
                    ci = i0 + block.shape[2] / 2.0
                    cj = j0 + block.shape[1] / 2.0
                    ck = k0 + block.shape[0] / 2.0
                    ras = refAffine @ np.array([ci, cj, ck, 1.0])
                    R, A, S = float(ras[0]), float(ras[1]), float(ras[2])
                    values = {}
                    for param, arr, inv in volInfo:
                        values[param] = self._sampleAtRas(arr, inv, (R, A, S), half, aggregation)
                    gid += 1
                    rows.append({"label": "glyph_%04d" % gid,
                                 "R": R, "A": A, "S": S, "values": values})
        return rows

    @staticmethod
    def _mapParam(rawList, param):
        target = PARAM_TARGET_RANGE[param]
        present = [v for v in rawList if v is not None]
        if not present:
            return [None] * len(rawList)
        vmin, vmax = min(present), max(present)
        span = (vmax - vmin) or 1.0
        if target is None:
            return list(rawList)
        lo, hi = target
        return [None if v is None else lo + (v - vmin) / span * (hi - lo) for v in rawList]

    def buildTable(self, sampledRows, assignments):
        assignedParams = [p for _, p in assignments]
        mapped = {}
        for param in assignedParams:
            raw = [r["values"].get(param) for r in sampledRows]
            mapped[param] = self._mapParam(raw, param)
        rows = []
        for idx, r in enumerate(sampledRows):
            out = {"label": r["label"], "R": r["R"], "A": r["A"], "S": r["S"]}
            for gp in ("N1", "N2", "size"):
                if gp in mapped:
                    v = mapped[gp][idx]
                    out[gp] = PARAM_FALLBACK[gp] if v is None else v
                else:
                    out[gp] = PARAM_FALLBACK[gp]
            if "color" in mapped:
                out["color"] = mapped["color"][idx]
            else:
                out["color"] = 0.0
            rows.append(out)
        return rows

    def writeTable(self, rows, path):
        header = ["label", "R_mm", "A_mm", "S_mm",
                  "param1_N1", "param2_N2", "param3_size", "param4_color"]
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            for r in rows:
                color = r["color"] if r["color"] is not None else 0.0
                w.writerow([r["label"], round(r["R"], 4), round(r["A"], 4),
                            round(r["S"], 4), round(r["N1"], 4), round(r["N2"], 4),
                            round(r["size"], 4), round(color, 4)])

    def readTable(self, path):
        cols = {"label": "label", "R": "R_mm", "A": "A_mm", "S": "S_mm",
                "N1": "param1_N1", "N2": "param2_N2",
                "size": "param3_size", "color": "param4_color"}
        rows = []
        with open(path, "r") as f:
            reader = csv.DictReader(f)
            missing = [v for v in cols.values() if v not in reader.fieldnames]
            if missing:
                raise ValueError("CSV missing columns: %s" % missing)
            for r in reader:
                rows.append({"label": r[cols["label"]].strip(),
                             "R": float(r[cols["R"]]), "A": float(r[cols["A"]]),
                             "S": float(r[cols["S"]]), "N1": float(r[cols["N1"]]),
                             "N2": float(r[cols["N2"]]), "size": float(r[cols["size"]]),
                             "color": float(r[cols["color"]])})
        return rows

    def _supertorusPolyData(self, n1, n2, ringRadius, crossFraction, resolution, colorValue):
        sup = vtk.vtkParametricSuperToroid()
        sup.SetN1(max(n1, 1e-3))
        sup.SetN2(max(n2, 1e-3))
        sup.SetRingRadius(ringRadius)
        sup.SetCrossSectionRadius(crossFraction * ringRadius)
        src = vtk.vtkParametricFunctionSource()
        src.SetParametricFunction(sup)
        src.SetUResolution(resolution)
        src.SetVResolution(resolution)
        src.SetScalarModeToNone()
        src.Update()
        poly = src.GetOutput()
        npts = poly.GetNumberOfPoints()
        arr = vtk.vtkFloatArray()
        arr.SetName(self.COLOR_ARRAY_NAME)
        arr.SetNumberOfTuples(npts)
        for i in range(npts):
            arr.SetValue(i, colorValue)
        poly.GetPointData().SetScalars(arr)
        return poly

    @staticmethod
    def _placePolyData(poly, R, A, S):
        tf = vtk.vtkTransform()
        tf.Translate(R, A, S)
        tpd = vtk.vtkTransformPolyDataFilter()
        tpd.SetInputData(poly)
        tpd.SetTransform(tf)
        tpd.Update()
        return tpd.GetOutput()

    def _folderItem(self):
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        sceneItem = shNode.GetSceneItemID()
        item = shNode.GetItemChildWithName(sceneItem, self.FOLDER_NAME)
        if not item:
            item = shNode.CreateFolderItem(sceneItem, self.FOLDER_NAME)
        return shNode, item

    def clear(self):
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        sceneItem = shNode.GetSceneItemID()
        folderItem = shNode.GetItemChildWithName(sceneItem, self.FOLDER_NAME)
        if not folderItem:
            return 0
        childIds = vtk.vtkIdList()
        shNode.GetItemChildren(folderItem, childIds, True)
        nodes = []
        for i in range(childIds.GetNumberOfIds()):
            node = shNode.GetItemDataNode(childIds.GetId(i))
            if node is not None:
                nodes.append(node)
        for node in nodes:
            slicer.mrmlScene.RemoveNode(node)
        shNode.RemoveItem(folderItem)
        return len(nodes)

    def renderGlyphs(self, rows, colorNode, baseRingRadius, crossFraction,
                     resolution, showLegend):
        colorVals = [r["color"] for r in rows if r["color"] is not None]
        if colorVals:
            cmin, cmax = min(colorVals), max(colorVals)
        else:
            cmin, cmax = 0.0, 1.0
        if cmax <= cmin:
            cmax = cmin + 1.0
        if colorNode is None:
            colorNode = slicer.util.getFirstNodeByName("Viridis", "vtkMRMLColorTableNode")
        if colorNode is None:
            colorNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLColorTableNode")
        self.clear()
        shNode, folderItem = self._folderItem()
        modelsLogic = slicer.modules.models.logic()
        created = []
        for r in rows:
            ringRadius = baseRingRadius * max(r["size"], 1e-2)
            colorValue = r["color"] if r["color"] is not None else cmin
            poly = self._supertorusPolyData(r["N1"], r["N2"], ringRadius,
                                            crossFraction, resolution, colorValue)
            poly = self._placePolyData(poly, r["R"], r["A"], r["S"])
            modelNode = modelsLogic.AddModel(poly)
            modelNode.SetName(r["label"])
            disp = modelNode.GetDisplayNode()
            disp.SetScalarVisibility(True)
            disp.SetActiveScalarName(self.COLOR_ARRAY_NAME)
            disp.SetAndObserveColorNodeID(colorNode.GetID())
            disp.SetScalarRangeFlag(0)
            disp.SetScalarRange(cmin, cmax)
            disp.SetOpacity(1.0)
            disp.SetBackfaceCulling(0)
            disp.SetVisibility2D(True)
            modelItem = shNode.GetItemByDataNode(modelNode)
            shNode.SetItemParent(modelItem, folderItem)
            created.append(modelNode)
        if showLegend and created:
            try:
                legend = slicer.modules.colors.logic().AddDefaultColorLegendDisplayNode(created[0])
                legend.SetTitleText("color")
                legend.SetLabelFormat("%.1f")
            except Exception as e:
                logging.warning("Color legend not added: %s", e)
        return created

    @staticmethod
    def _showVolumeRendering(ctNode):
        try:
            vrLogic = slicer.modules.volumerendering.logic()
            displayNode = vrLogic.GetFirstVolumeRenderingDisplayNode(ctNode)
            if displayNode is None:
                displayNode = vrLogic.CreateDefaultVolumeRenderingNodes(ctNode)
            displayNode.SetVisibility(True)
            preset = vrLogic.GetPresetByName("CT-Chest-Contrast-Enhanced")
            if preset is not None:
                displayNode.GetVolumePropertyNode().Copy(preset)
        except Exception as e:
            logging.warning("Volume rendering setup skipped: %s", e)

    @staticmethod
    def _resetThreeDView():
        try:
            lm = slicer.app.layoutManager()
            for i in range(lm.threeDViewCount):
                v = lm.threeDWidget(i).threeDView()
                v.resetFocalPoint()
                v.resetCamera()
        except Exception as e:
            logging.warning("Could not reset 3D view: %s", e)

    def run(self, ctNode, referenceNode, assignments, binSize=10,
            aggregation="mean_nonzero", skipBackground=True,
            backgroundThreshold=0.0, colorNode=None, baseRingRadius=6.0,
            crossFraction=0.5, resolution=32, showVolumeRendering=True,
            showLegend=True):
        if ctNode is None:
            raise ValueError("No CT volume selected.")
        if not assignments:
            raise ValueError("Assign at least one data volume to a glyph parameter.")
        if referenceNode is None:
            referenceNode = assignments[0][0]
        if binSize < 1:
            raise ValueError("Bin size must be >= 1.")
        sampled = self.sampleGlyphs(referenceNode, assignments, binSize,
                                    aggregation, skipBackground, backgroundThreshold)
        if not sampled:
            raise ValueError("No glyphs produced. Every bin was background; "
                             "lower the threshold or turn off 'skip background'.")
        rows = self.buildTable(sampled, assignments)
        self.renderGlyphs(rows, colorNode, baseRingRadius, crossFraction,
                          resolution, showLegend)
        if showVolumeRendering:
            self._showVolumeRendering(ctNode)
        self._resetThreeDView()
        return rows, len(rows)


class SupertoroidalGlyphsTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_SampleAndBuild()

    def test_SampleAndBuild(self):
        self.delayDisplay("Starting the test")
        ctArr = np.zeros((30, 40, 40), dtype=np.float32)
        ctArr[5:25, 8:32, 8:32] = 200
        ct = slicer.util.addVolumeFromArray(ctArr, np.diag([2., 2., 2., 1.]))
        ct.SetName("CT")
        petArr = np.zeros((30, 40, 40), dtype=np.float32)
        petArr[5:25, 8:32, 8:32] = np.linspace(1, 9, 20 * 24 * 24).reshape(20, 24, 24)
        pet1 = slicer.util.addVolumeFromArray(petArr, np.diag([2., 2., 2., 1.]))
        pet1.SetName("PET1")
        pet2 = slicer.util.addVolumeFromArray(petArr[::-1].copy(), np.diag([2., 2., 2., 1.]))
        pet2.SetName("PET2")
        logic = SupertoroidalGlyphsLogic()
        rows, n = logic.run(ctNode=ct, referenceNode=pet1,
                            assignments=[(pet1, "color"), (pet2, "N1")],
                            binSize=8, resolution=12)
        self.assertGreater(n, 0)
        for r in rows:
            for key in ("R", "A", "S", "N1", "N2", "size", "color"):
                self.assertIn(key, r)
        import tempfile
        tmp = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
        tmp.close()
        logic.writeTable(rows, tmp.name)
        back = logic.readTable(tmp.name)
        self.assertEqual(len(back), len(rows))
        self.delayDisplay("Test passed: %d glyphs, CSV round-trip OK" % n)

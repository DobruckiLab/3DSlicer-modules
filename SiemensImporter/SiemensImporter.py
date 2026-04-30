"""
SiemensImporter — 3D Slicer scripted module

Reads Siemens Inveon PET/CT datasets (.img + .img.hdr) and converts
every frame to a NIfTI file (.nii.gz), optionally loading the result
into the Slicer scene.

Installation
------------
In Slicer: Edit → Application Settings → Modules → Additional module paths
Add the directory that contains this file, then restart Slicer.
The module appears under Utilities → Siemens Inveon Importer.
"""

import logging
import os

import numpy as np
import qt
import ctk
import vtk
import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin


# ── Module metadata ───────────────────────────────────────────────────────────

class SiemensImporter(ScriptedLoadableModule):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent.title = "Siemens Inveon Importer"
        self.parent.categories = ["Utilities"]
        self.parent.dependencies = []
        self.parent.contributors = ["Wawrzyniec Dobrucki <br> Experimental Molecular Imaging Laboratory <br> University of Illinois at Urbana-Champaign"]
        self.parent.helpText = (
            "Import Siemens Inveon PET/CT datasets (.img + .img.hdr) "
            "and export each frame as a NIfTI file (.nii.gz).<br><br>"
            "Select the .img.hdr file, choose an output directory, "
            "configure the scaling options, then click Import."
        )
        self.parent.acknowledgementText = ""


# ── Widget ────────────────────────────────────────────────────────────────────

class SiemensImporterWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None

    def setup(self):
        super().setup()
        self.logic = SiemensImporterLogic()

        # ── Input ─────────────────────────────────────────────────────────
        inputBox = ctk.ctkCollapsibleButton()
        inputBox.text = "Input"
        self.layout.addWidget(inputBox)
        inputForm = qt.QFormLayout(inputBox)

        self.hdrSelector = ctk.ctkPathLineEdit()
        self.hdrSelector.filters = ctk.ctkPathLineEdit.Files
        self.hdrSelector.nameFilters = ["Siemens Inveon headers (*.hdr)"]
        self.hdrSelector.setToolTip("Select a Siemens Inveon .img.hdr file")
        inputForm.addRow("Header file (.img.hdr):", self.hdrSelector)

        self.infoLabel = qt.QLabel("—")
        self.infoLabel.setWordWrap(True)
        inputForm.addRow("Detected:", self.infoLabel)

        # ── Output ────────────────────────────────────────────────────────
        outputBox = ctk.ctkCollapsibleButton()
        outputBox.text = "Output"
        self.layout.addWidget(outputBox)
        outputForm = qt.QFormLayout(outputBox)

        self.outputDirButton = ctk.ctkDirectoryButton()
        outputForm.addRow("Output directory:", self.outputDirButton)

        self.applyScaleCheck = qt.QCheckBox(
            "Apply per-frame scale factor (recommended for PET)"
        )
        self.applyScaleCheck.setChecked(True)
        outputForm.addRow("", self.applyScaleCheck)

        self.applyCalibCheck = qt.QCheckBox(
            "Apply calibration factor — PET only (converts to activity units)"
        )
        self.applyCalibCheck.setChecked(True)
        outputForm.addRow("", self.applyCalibCheck)

        self.loadIntoSlicerCheck = qt.QCheckBox("Load result into Slicer scene")
        self.loadIntoSlicerCheck.setChecked(True)
        outputForm.addRow("", self.loadIntoSlicerCheck)

        # ── Run ───────────────────────────────────────────────────────────
        self.importButton = qt.QPushButton("Import")
        self.importButton.toolTip = "Parse header and convert frames to NIfTI"
        self.importButton.enabled = False
        self.layout.addWidget(self.importButton)

        self.statusLabel = qt.QPlainTextEdit()
        self.statusLabel.setReadOnly(True)
        self.statusLabel.setMaximumHeight(160)
        self.layout.addWidget(self.statusLabel)

        self.layout.addStretch(1)

        self.hdrSelector.currentPathChanged.connect(self._onHdrChanged)
        self.importButton.clicked.connect(self._onImport)

    def _onHdrChanged(self, path):
        if not path or not os.path.isfile(path):
            self.infoLabel.setText("—")
            self.importButton.enabled = False
            return
        try:
            hdr = self.logic.parseHeader(path)
            modality_name = {0: "PET", 1: "CT", 2: "SPECT"}.get(
                hdr.get("modality", -1), "Unknown"
            )
            unit_name = {1: "nCi/cc", 2: "Bq/cc", 3: "HU"}.get(
                hdr.get("calibration_units", 0), ""
            )
            self.infoLabel.setText(
                f"{modality_name}  |  "
                f"{hdr['x_dimension']} × {hdr['y_dimension']} × {hdr['z_dimension']} vox  |  "
                f"{hdr.get('pixel_size_x', 0):.3f} × "
                f"{hdr.get('pixel_size_y', 0):.3f} × "
                f"{hdr.get('pixel_size_z', 0):.3f} mm  |  "
                f"{hdr.get('total_frames', 1)} frame(s)"
                + (f"  |  {unit_name}" if unit_name else "")
            )
            self.outputDirButton.directory = os.path.dirname(path)
            self.importButton.enabled = True
        except Exception as exc:
            self.infoLabel.setText(f"Parse error: {exc}")
            self.importButton.enabled = False

    def _onImport(self):
        self.importButton.enabled = False
        self.statusLabel.setPlainText("Working…")
        slicer.app.processEvents()
        try:
            paths = self.logic.importDataset(
                hdrPath=self.hdrSelector.currentPath,
                outputDir=self.outputDirButton.directory,
                applyScale=self.applyScaleCheck.isChecked(),
                applyCalib=self.applyCalibCheck.isChecked(),
                loadIntoSlicer=self.loadIntoSlicerCheck.isChecked(),
            )
            self.statusLabel.setPlainText("Done.\n" + "\n".join(paths))
        except Exception as exc:
            import traceback
            self.statusLabel.setPlainText(f"Error: {exc}\n{traceback.format_exc()}")
        finally:
            self.importButton.enabled = True


# ── Logic ─────────────────────────────────────────────────────────────────────

class SiemensImporterLogic(ScriptedLoadableModuleLogic):
    """Stateless helper — all public methods are safe to call from scripts."""

    # data_type code → numpy dtype (see Siemens Inveon header spec)
    _DTYPE_MAP = {
        1: np.dtype("uint8"),
        2: np.dtype("<i2"),   # 16-bit int, little-endian
        3: np.dtype("<i4"),   # 32-bit int, little-endian
        4: np.dtype("<f4"),   # 32-bit float, little-endian (most PET images)
        5: np.dtype(">f4"),   # 32-bit float, big-endian
        6: np.dtype(">i2"),   # 16-bit int, big-endian
        7: np.dtype(">i4"),   # 32-bit int, big-endian
    }

    def parseHeader(self, hdrPath):
        """Parse a Siemens Inveon .img.hdr file.

        The header is plain-text key/value pairs separated by the token
        ``end_of_header``.  The first section holds global parameters; each
        subsequent section is a per-frame subheader.

        Returns a dict of global parameters with an additional ``'frames'``
        key whose value is a list of per-frame dicts.
        """
        sections, current = [], {}
        with open(hdrPath, "r", errors="replace") as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                if line == "end_of_header":
                    sections.append(current)
                    current = {}
                    continue
                key, _, rest = line.partition(" ")
                current[key] = rest.strip()
        if current:
            sections.append(current)

        header = dict(sections[0]) if sections else {}

        def _cast_int(d, *keys):
            for k in keys:
                if k in d:
                    try:
                        d[k] = int(d[k])
                    except (ValueError, TypeError):
                        pass

        def _cast_float(d, *keys):
            for k in keys:
                if k in d:
                    try:
                        d[k] = float(d[k])
                    except (ValueError, TypeError):
                        pass

        _cast_int(header,
                  "modality", "data_type", "total_frames", "time_frames",
                  "x_dimension", "y_dimension", "z_dimension", "w_dimension",
                  "calibration_units", "acquisition_mode", "model",
                  "number_of_bed_positions")
        _cast_float(header,
                    "pixel_size_x", "pixel_size_y", "pixel_size_z", "pixel_size",
                    "calibration_factor", "isotope_half_life",
                    "isotope_branching_fraction", "dose")

        frame_subheaders = []
        for sec in sections[1:]:
            fh_dict = dict(sec)
            _cast_int(fh_dict, "frame", "gate", "bed",
                      "event_type", "energy_window", "detector_panel")
            _cast_float(fh_dict,
                        "scale_factor", "frame_start", "frame_duration",
                        "deadtime_correction", "decay_correction",
                        "bed_offset", "vertical_bed_offset")
            # data_file_pointer is stored as two 32-bit ints: "high_part low_part"
            if "data_file_pointer" in fh_dict:
                parts = fh_dict["data_file_pointer"].split()
                if len(parts) >= 2:
                    fh_dict["data_file_offset"] = (int(parts[0]) << 32) | int(parts[1])
                elif parts:
                    fh_dict["data_file_offset"] = int(parts[0])
                else:
                    fh_dict["data_file_offset"] = 0
            frame_subheaders.append(fh_dict)

        header["frames"] = frame_subheaders
        return header

    def readFrame(self, imgPath, header, frameIdx=0):
        """Read one frame from a Siemens Inveon binary .img file.

        Voxels in the file are stored X-fastest (row-major in Fortran sense),
        so the returned array shape is ``(nz, ny, nx)`` — the [K, J, I] order
        that Slicer's ``updateVolumeFromArray`` expects.
        """
        dtype = self._DTYPE_MAP.get(header.get("data_type"))
        if dtype is None:
            raise ValueError(f"Unsupported data_type code: {header.get('data_type')}")

        nx = header["x_dimension"]
        ny = header["y_dimension"]
        nz = header["z_dimension"]
        n_voxels = nx * ny * nz

        frames = header.get("frames", [])
        if frameIdx < len(frames):
            byte_offset = frames[frameIdx].get("data_file_offset", 0)
        else:
            # Fallback: frames are packed sequentially
            byte_offset = frameIdx * n_voxels * dtype.itemsize

        with open(imgPath, "rb") as f:
            f.seek(byte_offset)
            raw = f.read(n_voxels * dtype.itemsize)

        return np.frombuffer(raw, dtype=dtype).reshape((nz, ny, nx))

    def importDataset(self, hdrPath, outputDir, applyScale=True,
                      applyCalib=True, loadIntoSlicer=True):
        """Convert all frames of a Siemens Inveon dataset to NIfTI.

        Parameters
        ----------
        hdrPath : str
            Path to the .img.hdr file.
        outputDir : str
            Directory where .nii.gz files will be written.
        applyScale : bool
            Multiply voxel values by the per-frame ``scale_factor``.
            Should be True for PET; scale_factor is usually 1.0 for CT.
        applyCalib : bool
            Multiply by ``calibration_factor`` (PET only).
            Together with applyScale this converts raw values to nCi/cc
            (or Bq/cc, depending on calibration_units).
        loadIntoSlicer : bool
            Keep the volume nodes in the scene after export.

        Returns
        -------
        list of str
            Paths of the written NIfTI files (one per frame).
        """
        header = self.parseHeader(hdrPath)

        # Binary image lives alongside the header: strip .hdr → .img
        imgPath = os.path.splitext(hdrPath)[0]
        if not os.path.isfile(imgPath):
            raise FileNotFoundError(f"Binary image file not found: {imgPath}")

        modality = header.get("modality", -1)
        total_frames = header.get("total_frames", 1)
        px = header.get("pixel_size_x", 1.0)
        py = header.get("pixel_size_y", 1.0)
        pz = header.get("pixel_size_z", 1.0)
        nx, ny, nz = header["x_dimension"], header["y_dimension"], header["z_dimension"]
        calib_factor = header.get("calibration_factor", 1.0)

        # Validate that the binary file is large enough for the declared volume.
        dtype = self._DTYPE_MAP.get(header.get("data_type"))
        if dtype is None:
            raise ValueError(f"Unsupported data_type code: {header.get('data_type')}")
        expected_bytes = nx * ny * nz * total_frames * dtype.itemsize
        actual_bytes = os.path.getsize(imgPath)
        if actual_bytes < expected_bytes:
            acq_mode = header.get("acquisition_mode", -1)
            if acq_mode == 9:
                raise ValueError(
                    f"This file contains raw CT projection data (acquisition_mode=9), "
                    f"not a reconstructed CT image.\n\n"
                    f"The header describes a {nx}×{ny}×{nz} reconstruction volume "
                    f"({expected_bytes:,} bytes) but the file is only {actual_bytes:,} bytes "
                    f"({actual_bytes / 1e6:.1f} MB).\n\n"
                    f"The Feldkamp-reconstructed CT volume was produced on the Siemens workstation "
                    f"but was not exported with this dataset. To obtain the reconstructed image:\n"
                    f"  1. Re-run Feldkamp reconstruction in Siemens IRW on the original .cat file, or\n"
                    f"  2. Use a compatible tool (e.g. RTK, TIGRE) with the projection data."
                )
            raise ValueError(
                f"File size mismatch: the header declares a {nx}×{ny}×{nz} volume "
                f"({expected_bytes:,} bytes) but {os.path.basename(imgPath)} is only "
                f"{actual_bytes:,} bytes. The file may be truncated."
            )

        # Base name for output files: strip .hdr (e.g. "name.pet.img" → "name.pet.img")
        base_name = os.path.splitext(os.path.basename(hdrPath))[0]
        os.makedirs(outputDir, exist_ok=True)
        output_paths = []

        for frame_idx in range(total_frames):
            data = self.readFrame(imgPath, header, frame_idx)
            frame_info = (header["frames"][frame_idx]
                          if frame_idx < len(header["frames"]) else {})

            # Apply scaling only where needed (avoids unnecessary float32 copy)
            needs_scale = applyScale and "scale_factor" in frame_info
            needs_calib = applyCalib and modality == 0  # calibration only for PET
            if needs_scale or needs_calib:
                data = data.astype(np.float32)
                if needs_scale:
                    data *= float(frame_info["scale_factor"])
                if needs_calib:
                    data *= calib_factor

            suffix = f"_frame{frame_idx:03d}" if total_frames > 1 else ""
            out_path = os.path.join(outputDir, f"{base_name}{suffix}.nii.gz")
            node_name = f"{base_name}{suffix}"

            volumeNode = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLScalarVolumeNode", node_name
            )

            # data shape is (nz, ny, nx) = [K, J, I] as required by Slicer
            slicer.util.updateVolumeFromArray(volumeNode, data)
            volumeNode.SetSpacing(px, py, pz)

            # Direction matrix: LPS → RAS requires negating X and Y axes.
            # The Siemens Inveon scanner uses an LPS-like coordinate system
            # (no reliable patient orientation in preclinical data).
            dirMatrix = vtk.vtkMatrix4x4()
            dirMatrix.Identity()
            dirMatrix.SetElement(0, 0, -1.0)  # I → -R (LPS X → RAS)
            dirMatrix.SetElement(1, 1, -1.0)  # J → -A (LPS Y → RAS)
            volumeNode.SetIJKToRASDirectionMatrix(dirMatrix)

            # Origin: place volume center at RAS (0, 0, 0).
            # With the direction matrix above, voxel (0,0,0) must sit at:
            #   ox = +px*(nx-1)/2,  oy = +py*(ny-1)/2,  oz = -pz*(nz-1)/2
            volumeNode.SetOrigin(
                 px * (nx - 1) / 2.0,
                 py * (ny - 1) / 2.0,
                -pz * (nz - 1) / 2.0,
            )

            # Apply a PET-appropriate colormap
            if modality == 0:
                displayNode = slicer.mrmlScene.AddNewNodeByClass(
                    "vtkMRMLScalarVolumeDisplayNode"
                )
                displayNode.SetAndObserveColorNodeID("vtkMRMLColorTableNodePET-Heat")
                displayNode.AutoWindowLevelOn()
                volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())

            slicer.util.saveNode(volumeNode, out_path)
            output_paths.append(out_path)
            logging.info("SiemensImporter: saved %s", out_path)

            if not loadIntoSlicer:
                slicer.mrmlScene.RemoveNode(volumeNode)

        return output_paths


# ── Tests ─────────────────────────────────────────────────────────────────────

class SiemensImporterTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_parseHeader()
        self.test_readFrame()

    def test_parseHeader(self):
        """Verify header parsing with a minimal synthetic header."""
        import tempfile
        content = (
            "version 002.000\n"
            "modality 0\n"
            "data_type 4\n"
            "x_dimension 64\n"
            "y_dimension 64\n"
            "z_dimension 32\n"
            "pixel_size_x 0.776\n"
            "pixel_size_y 0.776\n"
            "pixel_size_z 0.796\n"
            "total_frames 2\n"
            "calibration_factor 4150000.0\n"
            "calibration_units 1\n"
            "end_of_header\n"
            "frame 0\n"
            "scale_factor 0.001\n"
            "frame_start 0\n"
            "frame_duration 60\n"
            "data_file_pointer 0 0\n"
            "end_of_header\n"
            "frame 1\n"
            "scale_factor 0.002\n"
            "frame_start 60\n"
            "frame_duration 60\n"
            "data_file_pointer 0 524288\n"  # 64*64*32*4 bytes
            "end_of_header\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".hdr", delete=False) as f:
            f.write(content)
            tmp = f.name
        try:
            logic = SiemensImporterLogic()
            hdr = logic.parseHeader(tmp)
            assert hdr["modality"] == 0, "modality"
            assert hdr["x_dimension"] == 64, "x_dimension"
            assert hdr["total_frames"] == 2, "total_frames"
            assert len(hdr["frames"]) == 2, "frame count"
            assert hdr["frames"][0]["data_file_offset"] == 0, "frame 0 offset"
            assert hdr["frames"][1]["data_file_offset"] == 524288, "frame 1 offset"
            self.delayDisplay("test_parseHeader passed")
        finally:
            os.unlink(tmp)

    def test_readFrame(self):
        """Write a tiny synthetic .img file and verify readFrame round-trips."""
        import tempfile
        nx, ny, nz = 4, 4, 4
        rng = np.random.default_rng(0)
        expected = rng.standard_normal((nz, ny, nx)).astype(np.float32)

        with tempfile.NamedTemporaryFile(suffix=".img", delete=False) as f:
            f.write(expected.tobytes())
            img_path = f.name
        try:
            logic = SiemensImporterLogic()
            header = {
                "data_type": 4,
                "x_dimension": nx,
                "y_dimension": ny,
                "z_dimension": nz,
                "frames": [{"data_file_offset": 0}],
            }
            result = logic.readFrame(img_path, header, frameIdx=0)
            assert result.shape == (nz, ny, nx), f"shape mismatch: {result.shape}"
            assert np.allclose(result, expected), "data mismatch"
            self.delayDisplay("test_readFrame passed")
        finally:
            os.unlink(img_path)

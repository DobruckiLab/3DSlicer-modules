# Slicer Scripted Module Development Guide

Practical pitfalls and solutions learned from building the Radiomics Heterogeneity module.

---

## 1. Module Loading & Structure

### Problem: "Failed to instantiate scripted pythonqt class"
**Cause:** Slicer discovers every `*Module.py` file and tries to instantiate it as a loadable module.

**Solution:** Rename sibling source files with a non-`Module` suffix:
- `MyModuleLogic.py` → `MyModuleLogic_`
- `MyModuleGUI.py` → `MyModuleGUI_`

### Problem: "Warning, there is no UI for the module"
**Cause:** The widget class doesn't follow Slicer's naming convention, or the module doesn't use `ScriptedLoadableModuleWidget`.

**Solution:** Slicer discovers the widget class by the naming convention `<ModuleName>Widget` that inherits from `ScriptedLoadableModuleWidget`. The module class only handles metadata (title, categories, help text).

```python
# Module class — metadata only
class MyModule(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "..."
        parent.categories = ["Quantification"]
        # ... no setup() needed here

# Widget class — all UI logic
class MyModuleWidget(ScriptedLoadableModuleWidget):
    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)  # Always call first!
        self.layout = self.parent.layout()  # Get the widget's layout
        # ... add your UI elements to self.layout
```

### Problem: `ScriptedLoadableModuleWidget` not found
**Cause:** Not imported.

**Solution:**
```python
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleWidget,
)
```

---

## 2. Qt API — Common Pitfalls

### `ctk.ctkCollapsibleButton` doesn't accept label in constructor
**Wrong:** `ctk.ctkCollapsibleButton("Parameters")`
**Right:**
```python
btn = ctk.ctkCollapsibleButton()
btn.text = "Parameters"
```

### Qt properties are camelCase, not snake_case
**Wrong:** `widget.current_index`, `widget.read_only = True`, `widget.visible = False`
**Right:** `widget.currentIndex`, `widget.setFlags(flags & ~Qt.ItemIsEditable)`, `widget.setVisible(False)` / `widget.hide()`

### Qt method calls use camelCase
**Wrong:** `table.insert_row()`, `table.resize_columns_to_contents()`, `table.columnCount = 2`
**Right:** `table.insertRow()`, `table.resizeColumnsToContents()`, `table.setColumnCount(2)`

### `QProgressBar.visible` doesn't work — use `show()`/`hide()`
**Wrong:** `progressBar.visible = True`
**Right:** `progressBar.show()` / `progressBar.hide()`

### `QComboBox.currentText` and `currentIndex` are properties (no parentheses)
**Wrong:** `comboBox.currentIndex()`, `comboBox.currentText()`
**Right:** `comboBox.currentIndex`, `comboBox.currentText`

---

## 3. Slicer API — What Works and What Doesn't

### `slicer.util.array()` vs `slicer.util.arrayFromVolume()`
- `slicer.util.array(name)` expects a **string** (node name pattern), not a VTK node object
- `slicer.util.arrayFromVolume(node)` accepts a **VTK node object**

### `slicer.util.saveFileName()` doesn't exist
**Wrong:** `slicer.util.saveFileName(filter="CSV (*.csv)")`
**Right:** Use `qt.QFileDialog.getSaveFileName()`:
```python
path = qt.QFileDialog.getSaveFileName(
    None, "Export", "/path/to/default/file.csv", "CSV Files (*.csv)"
)
if path and path[0]:
    filepath = path[0]  # Qt returns a tuple (path, filter)
```

### `slicer.util.getNodesByClass()` returns a list, not a dict
**Right:**
```python
nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
# nodes is already a list — use it directly
```

### Scene node iteration
Use `GetFirstNodeByClass()` / `GetNextNodeByClass()` or `slicer.util.getNodesByClass()`.
The old `scene.GetNodeCollection().GetUniqueIDList()` pattern is **deprecated** and doesn't exist in newer Slicer versions.

### `slicer.vtkMRMLScene.NodeAddedEvent` is correct
This works in current Slicer:
```python
scene.AddObserver(slicer.vtkMRMLScene.NodeAddedEvent, self._onSceneChanged)
```

---

## 4. Segmentation API

### Getting the first segment
```python
segmentation = segment_node.GetSegmentation()
segment_id = segmentation.GetNthSegmentID(0)  # Use GetNthSegmentID, NOT GetSegmentIdBySegmentNumber
```

### Extracting binary labelmap from a segment
```python
# Requires reference volume node for proper geometry
mask_array = slicer.util.arrayFromSegmentBinaryLabelmap(
    segment_node, segment_id, referenceVolumeNode=volume_node
)
```

### `vtkSegmentationConverter.ExtractBinaryLabelmap` doesn't exist as a static method
Don't try to use it. Use `slicer.util.arrayFromSegmentBinaryLabelmap()` instead.

### `vtkSegmentationConverter` is in `vtkSegmentationCore`, not `vtk`
If you ever need it directly:
```python
import vtkSegmentationCore
converter = vtkSegmentationCore.vtkSegmentationConverter
```

---

## 5. NumPy / Array Indexing

### `np.where()` on 3D arrays returns `(x_coords, y_coords, z_coords)` for `[x, y, z]`-indexed arrays
**Wrong:**
```python
nz, ny, nx = np.where(mask)  # Unpacks as (x, y, z) but variables say (z, y, x)
z0, y0, x0 = nz[idx], ny[idx], nx[idx]  # Wrong!
```
**Right:**
```python
x_coords, y_coords, z_coords = np.where(mask)
x0, y0, z0 = x_coords[idx], y_coords[idx], z_coords[idx]
```

---

## 6. Feature Computation Structure

For a typical Slicer scripted module, structure as:

```
Module/
  MyModule.py          # Single file: module class + widget class + logic
  README.md            # Documentation
  Resources/
    Icons/
      MyModule.png     # Module icon (optional)
```

- Single file per module (no `__init__.py`, no sibling `*Module.py` files)
- All logic inlined or in sibling files with non-`Module` suffix
- Widget class named `<ModuleName>Widget` inheriting from `ScriptedLoadableModuleWidget`

---

## 7. Checklist Before Testing

- [ ] No `*Module.py` sibling files (rename with `_` suffix)
- [ ] No `__init__.py` in module directory
- [ ] No `VTK_OBSERVE_MODIFICATIONS` import (not available in all Slicer builds)
- [ ] Widget class named `<ModuleName>Widget` inheriting `ScriptedLoadableModuleWidget`
- [ ] `ScriptedLoadableModuleWidget.setup(self)` called first in `setup()`
- [ ] `self.layout = self.parent.layout()` to get the widget layout
- [ ] `ctk.ctkCollapsibleButton` — set `.text` property, not constructor argument
- [ ] Qt properties: `currentIndex`, `currentText`, `setColumnCount`, `insertRow`, etc.
- [ ] Use `slicer.util.arrayFromVolume()` not `slicer.util.array()` for VTK nodes
- [ ] Use `slicer.util.arrayFromSegmentBinaryLabelmap()` with `referenceVolumeNode`
- [ ] `np.where()` unpacks as `(x, y, z)` for `[x, y, z]`-indexed arrays
- [ ] Progress bar: use `show()`/`hide()`, not `.visible`

# 🧫 CellPax

**Standardised 96-Well Plate Fluorescence Quantifier for ImageJ/Fiji** 🌿

[![Version](https://img.shields.io/badge/Version-2.0.0-brightgreen.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![ImageJ](https://img.shields.io/badge/ImageJ-Fiji%201.54+-blue.svg)](https://fiji.sc/)
[![Language](https://img.shields.io/badge/Language-Jython-yellow.svg)](https://www.jython.org/)

---

CellPax quantifies fluorescence across an entire 96-well plate in one pass. It lays a fixed, standardised grid of 96 circular ROIs over the wells using SBS/ANSI plate geometry &mdash; **no thresholding, no spot detection** &mdash; and measures the mean intensity inside every well, negative controls included. Built for tobacco BY-2 cell packs in agroinfiltration screens, but works on any 96-well fluorescence plate. 🎯

## ✨ Features

- 📐 **Geometry-anchored grid** &mdash; 96 ROIs placed from the SBS/ANSI plate standard, not from signal, so empty/negative wells are measured too
- 🤖 **Auto plate detection** &mdash; finds the plate outline automatically, with a manual fallback
- 🖐️ **Manual well selector** &mdash; pick visible wells on an 8&times;12 layout and anchor the grid from two reference points when corners are cropped
- 🔁 **Batch mode** &mdash; point at a folder and process every image hands-free
- 👁️ **Visual QC overlay** &mdash; the 96-ROI grid stays on the image and is saved as a TIFF so you can verify alignment
- 📊 **Clean CSV output** &mdash; one row per well: image, well, mean intensity
- 🎚️ **Reference-box mode** &mdash; measure only a sub-region of the plate when you didn't use the full grid

## 📋 Requirements

- [Fiji](https://fiji.sc/) (ImageJ 1.54 or later)
- Jython interpreter (included with Fiji by default)

## 🚀 Installation

1. Download [`CellPax.py`](CellPax.py) from this repository
2. Place it anywhere on your computer (e.g. your Fiji `plugins/` folder, or any convenient location)
3. That's it &mdash; no additional dependencies required! 🎉

## 🔧 Usage

1. **Open** your 96-well plate fluorescence image in Fiji
2. **Run** CellPax via `Plugins > Macros > Run...` and select `CellPax.py`
3. **Fill in the parameter dialog** and click **OK**
4. **Confirm the grid** &mdash; the 96 ROIs are drawn on the image; check the circles line up with the wells before measuring
5. **Collect results** &mdash; a per-well CSV and an overlay TIFF are written next to your image

### 🗺️ Workflow

```
📂 Open Image
     │
     ▼
📐 Plate Detection   (auto outline  ──►  or manual well selector)
     │
     ▼
🎯 Grid Placement    (96 ROIs from SBS/ANSI geometry)
     │
     ▼
👁️ Confirm Grid      (optional: verify alignment, uncheck for batch)
     │
     ▼
📊 Measure Wells     (mean intensity per ROI)
     │
     ▼
💾 CSV + Overlay TIFF
```

## ⚙️ Parameters

| Parameter | Default | Description |
|---|---|---|
| **Input Image or Folder** | &mdash; | A single image, or a folder of images for batch processing |
| **Output Folder** | *(blank)* | Where results are saved &mdash; blank = next to the input image |
| **Plate Detection Mode** | `auto` | `auto` finds the plate outline; `manual` opens the well selector |
| **ROI Circle Diameter** | `0.70` | Circle diameter as a fraction of the 9&nbsp;mm well spacing (0.70 ≈ 6.3&nbsp;mm, sits inside each well) |
| **Measurement Order** | `row` | `row`: A1…A12, B1… &nbsp;·&nbsp; `column`: A1, B1…H1, A2… |
| **Confirm well grid** | `true` | Pause to verify grid alignment before measuring (uncheck for hands-free batch) |
| **Save overlay image** | `true` | Save a TIFF with the ROI grid burned in as a QC record |
| **Measurement Region** | `all` | `all`: every well · `reference_box`: only wells inside the two reference wells (manual mode) |

## 📐 Plate Detection

**Auto mode** *(recommended when the full plate is visible)* &mdash; the script blurs and auto-thresholds a duplicate to find the plate outline, draws a green box, and asks you to **Accept** or switch to manual. If detection looks wrong it falls back to manual automatically. Works best with strong contrast against a dark, non-reflective background.

**Manual mode** *(when corner wells are cropped or empty)* &mdash; pick the visible/usable wells on an 8&times;12 layout (with Select-All / Rows / Columns helpers), then click the centres of the two reference wells when prompted. The full 96-well grid is extrapolated from those two points using the known 9&nbsp;mm spacing. Selections don't need to be rectangular &mdash; the reference corners come from the bounding box of whatever you pick.

## 📤 Output

Saved next to the input image (or in the Output Folder if set):

| File | Contents |
|---|---|
| `<name>_well_measurements.csv` | One row per well: **Image**, **Well** (A1–H12), **Mean** intensity |
| `<name>_overlay.tif` | Original image with all 96 ROIs drawn (cyan = measured, grey = skipped) |

```
Image,Well,Mean
plate_001,A1,142.831200
plate_001,A2,12.004100
...
```

In `reference_box` mode the CSV is named with the reference labels, e.g. `plate_001_well_measurements_C1-E12.csv`.

## 🔬 Measurement Notes

Mean intensity is measured on a 32-bit floating-point duplicate of the input. For **RGB images**, the conversion uses luminance weighting (0.299R + 0.587G + 0.114B) &mdash; if your signal is in one channel (e.g. GFP), split channels first via `Image > Color > Split Channels` before running.

## 🧩 Adapting to Other Plate Formats

All plate geometry lives in the `CONFIG` block at the top of `CellPax.py`. To switch to a **384-well plate**, change only these values:

```python
WELL_A1_X_MM    = 12.13
WELL_A1_Y_MM    =  8.99
WELL_SPACING_MM =  4.50
NUM_ROWS        = 16
NUM_COLS        = 24
ROW_LABELS      = list("ABCDEFGHIJKLMNOP")
```

No other code changes are needed.

## 👤 Author

**Nick Eilmann**
- 📧 Email: nme122@ic.ac.uk
- 🐙 GitHub: [@ne1899](https://github.com/ne1899)

## 📝 Citation

If you use CellPax in your research, please cite:

> Eilmann N. (2026). CellPax: Standardised 96-Well Plate Fluorescence Quantifier for ImageJ/Fiji. https://github.com/ne1899/CellPax

## 📄 License

This project is licensed under the MIT License &mdash; see the [LICENSE](LICENSE) file for details.

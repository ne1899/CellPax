"""
WellPlateQuant.py
=================
Fiji/ImageJ Jython script for high-throughput quantification of fluorescence
in 96-well plates containing tobacco BY-2 cell packs used for agroinfiltration.

Approach
--------
A fixed, standardised grid of 96 circular ROIs is placed over the well
positions using SBS/ANSI plate geometry.  No spot detection, no thresholding
for ROI placement  -  just geometry anchored to the plate outline.

Negative-control wells (zero signal) are measured automatically alongside
positive wells because the grid is drawn regardless of signal.

Usage
-----
In Fiji:  Plugins > Macros > Run... -> select WellPlateQuant.py
A parameter dialog appears; fill in fields, click OK.

Single-image mode : point Input at an image file.
Batch mode        : point Input at a folder; all images in it are processed.

Output (saved next to each input image unless Output Folder is set)
------
  <name>_well_measurements.csv   per-well mean intensity table
  <name>_overlay.tif             original image with 96 ROIs drawn + labelled

Plate geometry reference : ANSI/SLAS 1-2004 through 4-2004 (SBS standard).

Author  : Nick Eilmann <nme122@ic.ac.uk>
Version : 2.0.0
Copyright (c) 2026 Nick Eilmann
License : MIT
"""

# -- Fiji Script Parameters ----------------------------------------------------
# Each #@ line becomes one widget in the auto-generated dialog.

#@ File    (label="Input Image or Folder",                                                           style="open")                       input_path
#@ File    (label="Output Folder  (leave blank = same folder as input)",                             style="directory", required=false)  output_dir
#@ String  (label="Plate Detection Mode",  choices={"auto","manual"},   value="auto")                detection_mode
#@ Double  (label="ROI Circle Diameter  (fraction of 9 mm well spacing, default 0.70)",              value=0.70, min=0.10, max=1.00, stepSize=0.05) roi_diameter_fraction
#@ String  (label="Measurement Order",     choices={"row","column"},    value="row")                 measurement_order
#@ Boolean (label="Confirm well grid for each image before measuring",  value=True)                  confirm_grid
#@ Boolean (label="Save overlay image (.tif)",                          value=True)                  save_overlay
#@ String  (label="Measurement Region",            choices={"all","reference_box"}, value="all")     measure_region

# -- Imports -------------------------------------------------------------------
from ij import IJ, ImagePlus, WindowManager
from ij.plugin import Duplicator
from ij.plugin.filter import ParticleAnalyzer
from ij.measure import ResultsTable, Measurements
from ij.plugin.frame import RoiManager
from ij.gui import (OvalRoi, Roi, TextRoi, Overlay,
                    NonBlockingGenericDialog, GenericDialog, WaitForUserDialog)
from ij.io import FileSaver
from java.awt import (Color, Font, GridLayout, BorderLayout, FlowLayout,
                      Dimension, Insets)
from java.awt.event import ActionListener
from java.lang import Double as JDouble
from javax.swing import (JDialog, JPanel, JToggleButton, JButton, JLabel,
                         JCheckBox, JOptionPane, BorderFactory, SwingConstants)
import os, math, sys

# -- CONFIG --------------------------------------------------------------------
# All plate geometry follows the SBS/ANSI 96-well standard (mm).
# To adapt to 384-well plates change only the values in this block:
#   WELL_A1_X_MM = 12.13, WELL_A1_Y_MM = 8.99, WELL_SPACING_MM = 4.50
#   NUM_ROWS = 16, NUM_COLS = 24, ROW_LABELS = list("ABCDEFGHIJKLMNOP")

PLATE_WIDTH_MM   = 127.76   # outer plate length (long axis), mm
PLATE_HEIGHT_MM  =  85.48   # outer plate width  (short axis), mm

WELL_A1_X_MM     =  14.38   # A1 centre distance from plate LEFT edge, mm
WELL_A1_Y_MM     =  11.24   # A1 centre distance from plate TOP  edge, mm
WELL_SPACING_MM  =   9.00   # centre-to-centre pitch (X and Y identical), mm

NUM_ROWS   = 8
NUM_COLS   = 12
ROW_LABELS = list("ABCDEFGH")

# -- Auto-detection parameters -------------------------------------------------
DETECT_BLUR_SIGMA       = 5.0        # Gaussian blur before thresholding
DETECT_THRESHOLD_METHOD = "Triangle" # try "Otsu" if this fails
DETECT_MIN_AREA_FRAC    = 0.25       # plate must cover >= 25% of image area
DETECT_ASPECT_RATIO     = PLATE_WIDTH_MM / PLATE_HEIGHT_MM  # ~1.495
DETECT_ASPECT_TOLERANCE = 0.45       # accept +/-45% deviation

# -- Measurement region --------------------------------------------------------
# Controls which wells are measured and written to the CSV.
#
#   "all"           : all 96 wells are measured regardless of detection mode.
#                     CSV is named  <name>_well_measurements.csv
#
#   "reference_box" : only wells that fall inside the inclusive rectangle
#                     defined by the two reference wells (top-left to
#                     bottom-right) are measured.  In auto-detection mode this
#                     has no effect and silently falls back to "all".
#                     CSV is named  <name>_well_measurements_<TL>-<BR>.csv
#                     e.g. plate_001_well_measurements_C1-E12.csv
#
# The overlay always shows all 96 ROIs; skipped wells are drawn in grey so
# alignment can be verified across the whole plate.
MEASURE_REGION = "all"   # overridden by the script parameter above

# -- Overlay visual style ------------------------------------------------------
OVERLAY_ROI_COLOR     = Color.CYAN           # selected wells
OVERLAY_SKIPPED_COLOR = Color.RED            # skipped wells
OVERLAY_LABEL_COLOR   = Color.YELLOW
OVERLAY_STROKE_WIDTH  = 1.5
OVERLAY_FONT_SIZE     = 11   # decrease for denser plates (e.g. 384-well)

# -- Batch processing ----------------------------------------------------------
IMAGE_EXTENSIONS = (".tif", ".tiff", ".png", ".jpg", ".jpeg", ".bmp")


# ==============================================================================
# Helper functions
# ==============================================================================

def well_label(row_idx, col_idx):
    """Return standard label e.g. well_label(0,0)='A1', well_label(7,11)='H12'."""
    return ROW_LABELS[row_idx] + str(col_idx + 1)


def all_wells():
    """Return set of all (row, col) index tuples."""
    s = set()
    for r in range(NUM_ROWS):
        for c in range(NUM_COLS):
            s.add((r, c))
    return s


def well_order(order):
    """
    Return list of (row, col) tuples in the requested order.
    order='row'    -> A1, A2 ... A12, B1 ... H12
    order='column' -> A1, B1 ... H1,  A2 ... H12
    """
    wells = []
    if order == "row":
        for r in range(NUM_ROWS):
            for c in range(NUM_COLS):
                wells.append((r, c))
    else:
        for c in range(NUM_COLS):
            for r in range(NUM_ROWS):
                wells.append((r, c))
    return wells


def compute_centers_from_box(px, py, pw, ph):
    """
    Compute 96 well-centre pixel coords from a detected plate bounding box.

    Converts mm distances (from SBS spec) to pixels using the ratio of
    bounding-box size to known plate outer dimensions.

    Parameters
    ----------
    px, py : top-left corner of the plate bounding box (pixels)
    pw, ph : width and height of the bounding box (pixels)

    Returns
    -------
    dict { (row_idx, col_idx) : (cx_px, cy_px) }
    """
    sx = float(pw) / PLATE_WIDTH_MM    # px per mm, X axis
    sy = float(ph) / PLATE_HEIGHT_MM   # px per mm, Y axis

    centers = {}
    for r in range(NUM_ROWS):
        for c in range(NUM_COLS):
            x_mm = WELL_A1_X_MM + c * WELL_SPACING_MM
            y_mm = WELL_A1_Y_MM + r * WELL_SPACING_MM
            centers[(r, c)] = (int(round(px + x_mm * sx)),
                               int(round(py + y_mm * sy)))
    return centers


def compute_centers_from_two_refs(r1, c1, px1, py1, r2, c2, px2, py2):
    """
    Compute all 96 well-centre pixel coords from two reference wells at known
    pixel positions.

    The two reference wells can be ANY pair on the plate -- they do NOT have to
    be A1 and H12.  The only requirement is that r2 > r1 and c2 > c1 so that
    the pixel-per-well step sizes in X and Y are both well-defined.

    Math
    ----
    dx_per_col = (px2 - px1) / (c2 - c1)   [pixels per one column step]
    dy_per_row = (py2 - py1) / (r2 - r1)   [pixels per one row step]

    For any well (r, c):
        x = px1 + (c - c1) * dx_per_col
        y = py1 + (r - r1) * dy_per_row

    Note: compute_centers_from_box is equivalent to calling this with
    r1=0, c1=0 at the A1 position and r2=7, c2=11 at the H12 position.

    Parameters
    ----------
    r1, c1 : row/col index of reference well 1 (top-left of the selection)
    px1,py1: pixel centre of reference well 1
    r2, c2 : row/col index of reference well 2 (bottom-right of the selection)
    px2,py2: pixel centre of reference well 2

    Returns
    -------
    dict { (row_idx, col_idx) : (cx_px, cy_px) }
    """
    dx_per_col = float(px2 - px1) / (c2 - c1)
    dy_per_row = float(py2 - py1) / (r2 - r1)

    centers = {}
    for r in range(NUM_ROWS):
        for c in range(NUM_COLS):
            centers[(r, c)] = (int(round(px1 + (c - c1) * dx_per_col)),
                               int(round(py1 + (r - r1) * dy_per_row)))
    return centers


# ==============================================================================
# Well selector dialog (Swing)
# ==============================================================================

def show_well_selector():
    """
    Show a modal Swing dialog with an 8x12 grid of toggle buttons.

    The user clicks wells to mark which ones are visible in their image.
    Preset buttons allow fast selection of all, none, specific rows, or
    specific columns.

    Validation (enforced before OK closes the dialog):
    - At least 2 wells must be selected.
    - The selection must span at least 2 rows (so row spacing is derivable).
    - The selection must span at least 2 columns (so column spacing is derivable).

    Returns
    -------
    set of (row_idx, col_idx) tuples, or None if the user cancelled.
    """
    SEL_COLOR   = Color.CYAN
    UNSEL_COLOR = Color(220, 220, 220)
    BTN_W       = 46
    BTN_H       = 28
    HDR_FONT    = Font("SansSerif", Font.BOLD,  10)
    BTN_FONT    = Font("SansSerif", Font.PLAIN,  9)

    # Mutable state shared across all inner listeners.
    # Using lists so inner classes can modify via index assignment (Jython 2.7).
    selected   = [[False] * NUM_COLS for _ in range(NUM_ROWS)]
    buttons    = [[None]  * NUM_COLS for _ in range(NUM_ROWS)]
    result     = [None]       # set to selected-set on OK, or stays None on Cancel
    dialog_ref = [None]       # filled once the JDialog is created below

    # --- apply state to a collection of (r,c) pairs -------------------------
    def _set_wells(pairs, state):
        for (r, c) in pairs:
            selected[r][c] = state
            if buttons[r][c] is not None:
                buttons[r][c].setSelected(state)
                buttons[r][c].setBackground(SEL_COLOR if state else UNSEL_COLOR)

    # --- toggle listener: one instance per button ---------------------------
    class ToggleListener(ActionListener):
        def __init__(self, r, c):
            self._r = r
            self._c = c
        def actionPerformed(self, event):
            btn = event.getSource()
            selected[self._r][self._c] = btn.isSelected()
            btn.setBackground(SEL_COLOR if btn.isSelected() else UNSEL_COLOR)

    # --- preset listeners ---------------------------------------------------
    class SelectAllListener(ActionListener):
        def actionPerformed(self, event):
            _set_wells(all_wells(), True)

    class ClearAllListener(ActionListener):
        def actionPerformed(self, event):
            _set_wells(all_wells(), False)

    class SelectRowsListener(ActionListener):
        def actionPerformed(self, event):
            panel = JPanel(GridLayout(NUM_ROWS, 1, 2, 2))
            checks = []
            for r in range(NUM_ROWS):
                cb = JCheckBox("Row " + ROW_LABELS[r])
                checks.append(cb)
                panel.add(cb)
            ret = JOptionPane.showConfirmDialog(
                dialog_ref[0], panel, "Select rows",
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE
            )
            if ret == JOptionPane.OK_OPTION:
                for r in range(NUM_ROWS):
                    if checks[r].isSelected():
                        _set_wells([(r, c) for c in range(NUM_COLS)], True)

    class SelectColsListener(ActionListener):
        def actionPerformed(self, event):
            panel = JPanel(GridLayout(2, 6, 2, 2))
            checks = []
            for c in range(NUM_COLS):
                cb = JCheckBox("Col " + str(c + 1))
                checks.append(cb)
                panel.add(cb)
            ret = JOptionPane.showConfirmDialog(
                dialog_ref[0], panel, "Select columns",
                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE
            )
            if ret == JOptionPane.OK_OPTION:
                for c in range(NUM_COLS):
                    if checks[c].isSelected():
                        _set_wells([(r, c) for r in range(NUM_ROWS)], True)

    # --- OK: validate then close -------------------------------------------
    class OKListener(ActionListener):
        def actionPerformed(self, event):
            sel = set()
            for r in range(NUM_ROWS):
                for c in range(NUM_COLS):
                    if selected[r][c]:
                        sel.add((r, c))

            if len(sel) < 2:
                JOptionPane.showMessageDialog(
                    dialog_ref[0],
                    "Please select at least 2 wells.",
                    "WellPlateQuant", JOptionPane.WARNING_MESSAGE
                )
                return

            rows_in_sel = [r for (r, c) in sel]
            cols_in_sel = [c for (r, c) in sel]

            if max(rows_in_sel) == min(rows_in_sel):
                JOptionPane.showMessageDialog(
                    dialog_ref[0],
                    "Selection must span at least 2 rows.\n"
                    "Row spacing cannot be derived from a single row.",
                    "WellPlateQuant", JOptionPane.WARNING_MESSAGE
                )
                return

            if max(cols_in_sel) == min(cols_in_sel):
                JOptionPane.showMessageDialog(
                    dialog_ref[0],
                    "Selection must span at least 2 columns.\n"
                    "Column spacing cannot be derived from a single column.",
                    "WellPlateQuant", JOptionPane.WARNING_MESSAGE
                )
                return

            result[0] = sel
            dialog_ref[0].dispose()

    class CancelListener(ActionListener):
        def actionPerformed(self, event):
            result[0] = None
            dialog_ref[0].dispose()

    # --- build the JDialog --------------------------------------------------
    dlg = JDialog()
    dlg.setTitle("WellPlateQuant  -  Select Visible Wells")
    dlg.setModal(True)
    dlg.setResizable(True)
    dlg.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE)
    dialog_ref[0] = dlg

    content = dlg.getContentPane()
    content.setLayout(BorderLayout(4, 4))

    instr = JLabel(
        "<html>Click wells that are <b>visible and usable</b> in your image"
        "  (cyan = selected).  Use presets below for common patterns."
        "  Click OK when done.</html>"
    )
    instr.setBorder(BorderFactory.createEmptyBorder(8, 10, 4, 10))
    content.add(instr, BorderLayout.NORTH)

    # 9 rows (header + 8 data) x 13 cols (header + 12 data)
    grid_panel = JPanel(GridLayout(NUM_ROWS + 1, NUM_COLS + 1, 2, 2))
    grid_panel.setBorder(BorderFactory.createEmptyBorder(2, 10, 2, 10))

    # top-left corner cell
    grid_panel.add(JLabel(""))

    # column headers
    for c in range(NUM_COLS):
        lbl = JLabel(str(c + 1), SwingConstants.CENTER)
        lbl.setFont(HDR_FONT)
        grid_panel.add(lbl)

    # row headers + toggle buttons
    for r in range(NUM_ROWS):
        row_lbl = JLabel(ROW_LABELS[r], SwingConstants.CENTER)
        row_lbl.setFont(HDR_FONT)
        grid_panel.add(row_lbl)

        for c in range(NUM_COLS):
            btn = JToggleButton(well_label(r, c))
            btn.setBackground(UNSEL_COLOR)
            btn.setPreferredSize(Dimension(BTN_W, BTN_H))
            btn.setFont(BTN_FONT)
            btn.setMargin(Insets(0, 0, 0, 0))
            btn.setOpaque(True)   # required on macOS for background to show
            btn.addActionListener(ToggleListener(r, c))
            buttons[r][c] = btn
            grid_panel.add(btn)

    content.add(grid_panel, BorderLayout.CENTER)

    # south panel: presets on left, OK/Cancel on right
    south = JPanel(BorderLayout(4, 4))
    south.setBorder(BorderFactory.createEmptyBorder(2, 10, 10, 10))

    preset_panel = JPanel(FlowLayout(FlowLayout.LEFT, 4, 2))
    btn_all  = JButton("Select All")
    btn_none = JButton("Clear All")
    btn_rows = JButton("Select Rows...")
    btn_cols = JButton("Select Columns...")
    btn_all.addActionListener(SelectAllListener())
    btn_none.addActionListener(ClearAllListener())
    btn_rows.addActionListener(SelectRowsListener())
    btn_cols.addActionListener(SelectColsListener())
    for b in (btn_all, btn_none, btn_rows, btn_cols):
        preset_panel.add(b)
    south.add(preset_panel, BorderLayout.CENTER)

    ok_panel = JPanel(FlowLayout(FlowLayout.RIGHT, 6, 2))
    btn_ok     = JButton("OK")
    btn_cancel = JButton("Cancel")
    btn_ok.addActionListener(OKListener())
    btn_cancel.addActionListener(CancelListener())
    ok_panel.add(btn_ok)
    ok_panel.add(btn_cancel)
    south.add(ok_panel, BorderLayout.SOUTH)

    content.add(south, BorderLayout.SOUTH)

    dlg.pack()
    dlg.setLocationRelativeTo(None)
    dlg.setVisible(True)   # blocks the script thread until dialog is disposed

    return result[0]


# ==============================================================================
# Plate detection
# ==============================================================================

def auto_detect_plate(imp):
    """
    Attempt to locate the plate rectangle by thresholding the image.

    Algorithm
    ---------
    1. Duplicate -> 8-bit grayscale.
    2. Gaussian blur (DETECT_BLUR_SIGMA px) so the plate body is a uniform blob.
    3. Auto-threshold (DETECT_THRESHOLD_METHOD) -> binary mask.
    4. Fill holes + two dilate passes to close mask gaps.
    5. Analyze Particles to find the plate body (min area = DETECT_MIN_AREA_FRAC).
    6. Accept the largest region if its aspect ratio is within DETECT_ASPECT_TOLERANCE.
    7. Retry on the inverted image (dark plate on bright background).

    Returns
    -------
    (x, y, w, h) bounding box in pixels, or None if detection fails.
    """
    W = imp.getWidth()
    H = imp.getHeight()
    min_area_px = int(W * H * DETECT_MIN_AREA_FRAC)

    def _try_detect(source_imp, inverted):
        temp = Duplicator().run(source_imp)
        if temp.getType() != ImagePlus.GRAY8:
            IJ.run(temp, "8-bit", "")
        if inverted:
            IJ.run(temp, "Invert", "")
        IJ.run(temp, "Gaussian Blur...", "sigma=" + str(DETECT_BLUR_SIGMA))
        IJ.run(temp, "Auto Threshold", "method=" + DETECT_THRESHOLD_METHOD + " white")
        IJ.run(temp, "Fill Holes", "")
        IJ.run(temp, "Dilate", "")
        IJ.run(temp, "Dilate", "")

        rm = RoiManager(False)
        pa = ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.SHOW_NONE,
            Measurements.AREA,
            ResultsTable(),
            min_area_px,
            JDouble.POSITIVE_INFINITY,
            0.0, 1.0
        )
        pa.setRoiManager(rm)
        pa.analyze(temp)
        temp.close()

        if rm.getCount() == 0:
            return None

        best_b = None
        best_a = 0
        for i in range(rm.getCount()):
            b = rm.getRoi(i).getBounds()
            a = b.width * b.height
            if a > best_a:
                best_a = a
                best_b = b
        return best_b

    for inverted in (False, True):
        b = _try_detect(imp, inverted)
        if b is None:
            continue
        aspect = float(b.width) / max(b.height, 1)
        if abs(aspect - DETECT_ASPECT_RATIO) <= DETECT_ASPECT_TOLERANCE:
            tag = " (inverted)" if inverted else ""
            IJ.log("  Auto-detect%s: bbox x=%d y=%d w=%d h=%d  aspect=%.3f"
                   % (tag, b.x, b.y, b.width, b.height, aspect))
            return (b.x, b.y, b.width, b.height)
        else:
            IJ.log("  Auto-detect%s: aspect %.3f outside tolerance -- skipping."
                   % (" (inverted)" if inverted else "", aspect))

    return None


def manual_detect_plate(imp):
    """
    Full manual registration workflow:
      1. Show well-selector dialog so the user marks which wells are visible.
      2. Derive the top-left and bottom-right reference wells from the
         bounding box of the selection.
      3. Prompt the user to click those two well centres on the image.
      4. Return everything needed to compute the full grid.

    Reference-corner derivation
    ---------------------------
    top-left  reference = (min_row, min_col) of the selected set
    bottom-right reference = (max_row, max_col) of the selected set

    This works for any shape of selection -- rectangular, scattered, partial
    rows/columns -- as long as the bounding box spans >=1 row and >=1 column
    (validated in show_well_selector before it returns).

    Returns
    -------
    (r1, c1, px1, py1, r2, c2, px2, py2, selected_wells) or None if cancelled.
    """
    IJ.showStatus("WellPlateQuant: select visible wells...")
    selected_wells = show_well_selector()
    if selected_wells is None:
        return None

    # Derive the two reference wells from the selection bounding box.
    rows_in_sel = sorted(set(r for (r, c) in selected_wells))
    cols_in_sel = sorted(set(c for (r, c) in selected_wells))
    r1, c1 = rows_in_sel[0],  cols_in_sel[0]    # top-left
    r2, c2 = rows_in_sel[-1], cols_in_sel[-1]   # bottom-right

    ref1_lbl = well_label(r1, c1)
    ref2_lbl = well_label(r2, c2)
    IJ.log("  Reference wells derived from selection: %s (top-left), %s (bottom-right)"
           % (ref1_lbl, ref2_lbl))

    IJ.setTool("point")
    imp.show()
    WindowManager.setCurrentWindow(imp.getWindow())

    WaitForUserDialog(
        "WellPlateQuant  -  Step 1 / 2",
        "Point tool is active.  Click the CENTRE of well " + ref1_lbl + ".\n\n"
        "Click OK when done."
    ).show()

    roi1 = imp.getRoi()
    if roi1 is None:
        IJ.error("WellPlateQuant", "No point selected for well " + ref1_lbl + ".  Aborting.")
        return None
    b1  = roi1.getBounds()
    px1 = b1.x
    py1 = b1.y
    imp.setRoi(None)

    WaitForUserDialog(
        "WellPlateQuant  -  Step 2 / 2",
        "Now click the CENTRE of well " + ref2_lbl + ".\n\n"
        "Click OK when done."
    ).show()

    roi2 = imp.getRoi()
    if roi2 is None:
        IJ.error("WellPlateQuant", "No point selected for well " + ref2_lbl + ".  Aborting.")
        return None
    b2  = roi2.getBounds()
    px2 = b2.x
    py2 = b2.y
    imp.setRoi(None)

    IJ.log("  %s=(%d,%d)  %s=(%d,%d)" % (ref1_lbl, px1, py1, ref2_lbl, px2, py2))
    return (r1, c1, px1, py1, r2, c2, px2, py2, selected_wells)


# ==============================================================================
# ROI / overlay
# ==============================================================================

def build_well_overlay(well_centers, radius_px, selected_wells=None):
    """
    Build a Fiji Overlay with 96 circular ROIs.

    Appearance
    ----------
    selected_wells=None (auto mode, all measured):
        All 96 wells drawn in cyan with yellow labels.
    selected_wells=set (manual mode):
        Wells IN the set   -> cyan circle + yellow label (measured).
        Wells NOT in the set -> grey circle, thinner stroke, no label (skipped).

    Parameters
    ----------
    well_centers    : dict { (row, col) : (cx, cy) }
    radius_px       : float  circle radius in pixels
    selected_wells  : set of (row, col) or None

    Returns
    -------
    Overlay
    """
    overlay  = Overlay()
    font     = Font("SansSerif", Font.BOLD, OVERLAY_FONT_SIZE)
    r_int    = int(radius_px)
    show_all = (selected_wells is None)

    for (row, col), (cx, cy) in sorted(well_centers.items()):
        x0   = cx - r_int
        y0   = cy - r_int
        diam = 2 * r_int
        oval = OvalRoi(x0, y0, diam, diam)
        lbl  = well_label(row, col)
        oval.setName(lbl)

        is_sel = show_all or ((row, col) in selected_wells)
        if is_sel:
            oval.setStrokeColor(OVERLAY_ROI_COLOR)
            oval.setStrokeWidth(OVERLAY_STROKE_WIDTH)
            overlay.add(oval)
            text = TextRoi(float(x0 + 2), float(y0 + 2), lbl)
            text.setFont(font)
            text.setStrokeColor(OVERLAY_LABEL_COLOR)
            overlay.add(text)
        else:
            oval.setStrokeColor(OVERLAY_SKIPPED_COLOR)
            oval.setStrokeWidth(1.0)
            overlay.add(oval)
            # no label on skipped wells to keep the image readable

    return overlay


# ==============================================================================
# Measurement
# ==============================================================================

def measure_well(ip, cx, cy, radius_px):
    """
    Return the mean pixel intensity inside a circular ROI.

    Parameters
    ----------
    ip        : ImageProcessor (32-bit)
    cx, cy    : int  centre coordinates in pixels
    radius_px : float  circle radius in pixels

    Returns
    -------
    float  mean pixel intensity
    """
    r    = int(radius_px)
    oval = OvalRoi(cx - r, cy - r, 2 * r, 2 * r)
    ip.setRoi(oval)
    s = ip.getStatistics()
    return float(s.mean)


# ==============================================================================
# Main per-image pipeline
# ==============================================================================

def process_image(imp, out_dir, det_mode, roi_frac, meas_order,
                  do_confirm, do_overlay, meas_region):
    """
    Full analysis pipeline for a single ImagePlus.

    Steps
    -----
    1.  Detect (or register) the plate to get well centres.
    2.  Compute ROI radius from the measured pixel pitch.
    3.  Populate the ROI Manager with 96 labelled OvalRois.
    4.  Optionally show a grid-preview dialog for the user to confirm.
    5.  Measure wells (all 96, or only those inside the reference box).
    6.  Show the Fiji ResultsTable.
    7.  Save CSV.
    8.  Optionally save a flattened overlay TIFF.
    9.  Leave the overlay visible on the image for instant inspection.

    Parameters
    ----------
    imp         : ImagePlus (must already be shown)
    out_dir     : str   where output files are written
    det_mode    : str   "auto" | "manual"
    roi_frac    : float circle diameter as fraction of well spacing
    meas_order  : str   "row" | "column"
    do_confirm  : bool  show grid-preview dialog before measuring
    do_overlay  : bool  save overlay TIFF
    meas_region : str   "all" = measure all 96 wells
                        "reference_box" = measure only wells inside the
                        rectangle spanned by the two reference wells
    """
    # -- Output base name --------------------------------------------------
    fi = imp.getOriginalFileInfo()
    if fi is not None and fi.fileName:
        base = os.path.splitext(fi.fileName)[0]
    else:
        title = imp.getTitle()
        base  = title.rsplit(".", 1)[0] if "." in title else title

    IJ.log("WellPlateQuant -- Processing: " + base)

    # -- Step 1: Plate detection / registration ----------------------------
    well_centers = None
    ref_box      = None   # set of (r,c) inside reference rectangle, or None
    ref_r1 = ref_c1 = ref_r2 = ref_c2 = None

    if det_mode == "auto":
        IJ.showStatus("WellPlateQuant: detecting plate (auto)...")
        box = auto_detect_plate(imp)

        if box is not None:
            px, py, pw, ph = box
            aspect = float(pw) / max(ph, 1)

            prev = Overlay()
            r_rect = Roi(px, py, pw, ph)
            r_rect.setStrokeColor(Color.GREEN)
            r_rect.setStrokeWidth(2)
            prev.add(r_rect)
            imp.setOverlay(prev)
            imp.updateAndDraw()

            dlg = NonBlockingGenericDialog("WellPlateQuant  -  Plate Detection")
            dlg.addMessage(
                "Green rectangle = auto-detected plate region.\n"
                "  x=%-5d  y=%-5d  w=%-5d  h=%d\n"
                "  Aspect ratio: %.3f  (expected %.3f)\n\n"
                "Accept to place well grid, or switch to Manual mode\n"
                "if the rectangle does not match the plate." %
                (px, py, pw, ph, aspect, DETECT_ASPECT_RATIO)
            )
            dlg.setOKLabel("Accept")
            dlg.setCancelLabel("Use Manual Mode")
            dlg.showDialog()
            imp.setOverlay(None)
            imp.updateAndDraw()

            if dlg.wasOKed():
                well_centers = compute_centers_from_box(px, py, pw, ph)
                if meas_region == "reference_box":
                    IJ.log("  Note: 'reference_box' has no effect in auto mode -- measuring all 96.")
            else:
                det_mode = "manual"
        else:
            IJ.log("  Auto-detection found no plate -- switching to manual mode.")
            IJ.showMessage("WellPlateQuant",
                           "Automatic plate detection failed for:\n" + base +
                           "\n\nYou will be asked to select wells and click reference points.")
            det_mode = "manual"

    if det_mode == "manual" or well_centers is None:
        IJ.showStatus("WellPlateQuant: manual plate registration...")
        pts = manual_detect_plate(imp)
        if pts is None:
            IJ.log("  Manual registration cancelled -- skipping: " + base)
            return
        ref_r1, ref_c1, px1, py1, ref_r2, ref_c2, px2, py2, _sel = pts
        well_centers = compute_centers_from_two_refs(ref_r1, ref_c1, px1, py1,
                                                     ref_r2, ref_c2, px2, py2)
        # Build the inclusive rectangle of wells between the two reference wells
        ref_box = set()
        for r in range(ref_r1, ref_r2 + 1):
            for c in range(ref_c1, ref_c2 + 1):
                ref_box.add((r, c))

    # -- Step 2: ROI radius from pixel pitch -------------------------------
    # Average the A1->A2 (column) and A1->B1 (row) step distances.
    c00 = well_centers[(0, 0)]
    c01 = well_centers[(0, 1)]
    c10 = well_centers[(1, 0)]
    dx_col  = math.sqrt((c01[0]-c00[0])**2 + (c01[1]-c00[1])**2)
    dy_row  = math.sqrt((c10[0]-c00[0])**2 + (c10[1]-c00[1])**2)
    spacing = (dx_col + dy_row) / 2.0
    radius  = roi_frac * spacing / 2.0
    IJ.log("  Pixel pitch: col=%.1f px  row=%.1f px | ROI radius=%.1f px  (frac=%.2f)"
           % (dx_col, dy_row, radius, roi_frac))

    # -- Step 3: Populate ROI Manager --------------------------------------
    rm = RoiManager.getInstance()
    if rm is None:
        rm = RoiManager()
    rm.reset()

    r_int = int(radius)
    for (row, col), (cx, cy) in sorted(well_centers.items()):
        oval = OvalRoi(cx - r_int, cy - r_int, 2 * r_int, 2 * r_int)
        oval.setName(well_label(row, col))
        rm.addRoi(oval)

    # Determine the active measurement set for overlay + filtering.
    # ref_box is only populated in manual mode; auto mode always measures all 96.
    wells_to_measure = ref_box if (meas_region == "reference_box" and ref_box is not None) else None
    grid_overlay     = build_well_overlay(well_centers, radius, wells_to_measure)

    # -- Step 4: Optional grid preview -------------------------------------
    if do_confirm:
        imp.setOverlay(grid_overlay)
        imp.updateAndDraw()

        confirm = NonBlockingGenericDialog("WellPlateQuant  -  Grid Preview: " + base)
        confirm.addMessage(
            "Cyan circles = wells to be measured.\n"
            "Grey circles  = wells that will be skipped (reference_box mode).\n"
            "Yellow labels = well IDs.\n\n"
            "Zoom in to verify circles sit inside each well.\n\n"
            "Click 'Measure' to proceed, or 'Cancel' to skip this image."
        )
        confirm.setOKLabel("Measure")
        confirm.setCancelLabel("Cancel")
        confirm.showDialog()

        imp.setOverlay(None)
        imp.updateAndDraw()

        if confirm.wasCanceled():
            IJ.log("  Measurement cancelled by user for: " + base)
            return

    # -- Step 5: Measure ---------------------------------------------------
    # 32-bit duplicate so measurements are bit-depth-agnostic.
    # For RGB images the 32-bit conversion uses luminance weighting.
    imp_work = Duplicator().run(imp)
    IJ.run(imp_work, "32-bit", "")
    ip = imp_work.getProcessor()

    ordered  = well_order(meas_order)
    csv_rows = []

    for (row, col) in ordered:
        if wells_to_measure is not None and (row, col) not in wells_to_measure:
            continue
        cx, cy = well_centers[(row, col)]
        mean   = measure_well(ip, cx, cy, radius)
        csv_rows.append({
            "well": well_label(row, col),
            "mean": mean,
        })

    imp_work.close()

    n_measured = len(csv_rows)
    IJ.log("  Measured %d wells." % n_measured)

    # -- Step 6: Fiji ResultsTable -----------------------------------------
    rt = ResultsTable()
    for i, rd in enumerate(csv_rows):
        rt.incrementCounter()
        rt.setValue("Well", i, rd["well"])
        rt.setValue("Mean", i, rd["mean"])
    rt.show("WellPlateQuant  -  " + base)

    # -- Step 7: Save CSV --------------------------------------------------
    if meas_region == "reference_box" and ref_r1 is not None:
        region_tag = "_" + well_label(ref_r1, ref_c1) + "-" + well_label(ref_r2, ref_c2)
    else:
        region_tag = ""
    csv_path = os.path.join(out_dir, base + "_well_measurements" + region_tag + ".csv")
    with open(csv_path, "w") as fh:
        fh.write("Image,Well,Mean\n")
        for rd in csv_rows:
            fh.write("%s,%s,%.6f\n" % (base, rd["well"], rd["mean"]))
    IJ.log("  CSV  -> " + csv_path)

    # -- Step 8: Save overlay TIFF -----------------------------------------
    if do_overlay:
        imp.setOverlay(grid_overlay)
        flat = imp.flatten()
        if flat is not None:
            overlay_path = os.path.join(out_dir, base + "_overlay.tif")
            ok = FileSaver(flat).saveAsTiff(overlay_path)
            flat.close()
            if ok:
                IJ.log("  TIF  -> " + overlay_path)
            else:
                IJ.log("  WARNING: overlay save failed -> " + overlay_path)
        else:
            IJ.log("  WARNING: flatten() returned null, overlay not saved.")

    # Leave overlay visible on the image for instant inspection
    imp.setOverlay(grid_overlay)
    imp.updateAndDraw()

    IJ.showStatus("WellPlateQuant: done  -  " + base)
    IJ.log("WellPlateQuant -- Done: " + base)


# ==============================================================================
# Entry point  -  single image or batch
# ==============================================================================

_input        = str(input_path)
_outdir       = str(output_dir) if (output_dir is not None and
                                    str(output_dir) not in ("null", "")) else None
_det          = str(detection_mode)
_frac         = float(roi_diameter_fraction)
_order        = str(measurement_order)
_confirm      = bool(confirm_grid)
_overlay      = bool(save_overlay)
_meas_region  = str(measure_region)

if os.path.isdir(_input):
    # -- Batch mode ----------------------------------------------------------
    IJ.log("WellPlateQuant: Batch mode  -  " + _input)
    files = sorted([f for f in os.listdir(_input)
                    if f.lower().endswith(IMAGE_EXTENSIONS)])
    if not files:
        IJ.error("WellPlateQuant", "No supported image files found in:\n" + _input)
        sys.exit()

    IJ.log("  Found %d image(s)." % len(files))
    for fname in files:
        fpath = os.path.join(_input, fname)
        imp   = IJ.openImage(fpath)
        if imp is None:
            IJ.log("  Cannot open: " + fname + "  -  skipped.")
            continue
        imp.show()
        eff_out = _outdir if _outdir else _input
        process_image(imp, eff_out, _det, _frac, _order,
                      _confirm, _overlay, _meas_region)
        imp.close()

    IJ.log("WellPlateQuant: Batch complete.")

else:
    # -- Single-image mode ---------------------------------------------------
    imp = IJ.openImage(_input)
    if imp is None:
        IJ.error("WellPlateQuant", "Cannot open image:\n" + _input)
        sys.exit()
    imp.show()

    img_dir = os.path.dirname(_input)
    eff_out = _outdir if _outdir else (img_dir if img_dir else os.getcwd())
    process_image(imp, eff_out, _det, _frac, _order,
                  _confirm, _overlay, _meas_region)

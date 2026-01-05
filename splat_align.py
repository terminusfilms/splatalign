#!/usr/bin/env python3
"""
SplatAlign - ICP Alignment for 3D Gaussian Splatting PLY Files

Aligns two 3DGS captures of the same location taken at different times.
Outputs a 4x4 transformation matrix AND optionally exports an aligned PLY
with the transform baked in (ready to use alongside the primary in any viewer).

Based on proven alignment pipeline achieving ~5cm median error.

Usage:
    python splat_align.py                    # Launch GUI
    python splat_align.py --cli primary.ply secondary.ply  # Command line
    python splat_align.py --cli primary.ply secondary.ply --bake  # With baked PLY
"""

import numpy as np
from plyfile import PlyData, PlyElement
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
import json
import sys
import os
from pathlib import Path
from datetime import datetime


def load_ply_fast(filepath, sample_size=300000, z_percentile_max=40, log_func=print):
    """
    Fast PLY loading with ground filtering.
    
    Ground filtering (bottom 40% of Z) excludes seasonal vegetation
    for more stable alignment between captures.
    """
    log_func(f"Loading {Path(filepath).name}...")
    plydata = PlyData.read(filepath)
    v = plydata['vertex'].data

    total = len(v)
    log_func(f"  Total points: {total:,}")

    # Sample indices first
    np.random.seed(42)
    idx = np.random.choice(total, min(sample_size * 2, total), replace=False)

    # Extract XYZ for sampled points
    x = v['x'][idx].astype(np.float64)
    y = v['y'][idx].astype(np.float64)
    z = v['z'][idx].astype(np.float64)

    # Filter to ground/low structure (bottom percentile of Z)
    z_thresh = np.percentile(z, z_percentile_max)
    log_func(f"  Z threshold (p{z_percentile_max}): {z_thresh:.2f}")

    mask = z <= z_thresh
    x, y, z = x[mask], y[mask], z[mask]

    # Limit to sample_size
    if len(x) > sample_size:
        keep = np.random.choice(len(x), sample_size, replace=False)
        x, y, z = x[keep], y[keep], z[keep]

    points = np.column_stack([x, y, z])
    log_func(f"  Ground points: {len(points):,}")
    return points


def icp_multiscale(source, target, log_func=print):
    """
    Multi-scale ICP with outlier rejection.
    
    Runs ICP at progressively finer distance thresholds:
    5m → 1m → 0.3m → 0.1m
    """

    def run_icp(src, tgt, max_iter, tol, outlier_thresh):
        tree = cKDTree(tgt)
        R_total = np.eye(3)
        t_total = np.zeros(3)
        pts = src.copy()
        prev_err = float('inf')

        for i in range(max_iter):
            dist, idx = tree.query(pts, k=1)

            # Outlier rejection
            inlier = dist < outlier_thresh
            if np.sum(inlier) < 100:
                break

            err = np.mean(dist[inlier])
            if abs(prev_err - err) < tol:
                break
            prev_err = err

            # ICP step on inliers
            p_src = pts[inlier]
            p_tgt = tgt[idx[inlier]]

            c_src = p_src.mean(axis=0)
            c_tgt = p_tgt.mean(axis=0)

            H = (p_src - c_src).T @ (p_tgt - c_tgt)
            U, S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:
                Vt[-1] *= -1
                R = Vt.T @ U.T

            t = c_tgt - R @ c_src
            pts = (R @ pts.T).T + t
            R_total = R @ R_total
            t_total = R @ t_total + t

        T = np.eye(4)
        T[:3, :3] = R_total
        T[:3, 3] = t_total
        return T, err

    log_func("\n=== Stage 1: Coarse (5m threshold) ===")
    T1, e1 = run_icp(source, target, 100, 1e-5, 5.0)
    src1 = (T1[:3,:3] @ source.T).T + T1[:3,3]
    log_func(f"  Error: {e1:.4f}m")

    log_func("\n=== Stage 2: Medium (1m threshold) ===")
    T2, e2 = run_icp(src1, target, 150, 1e-6, 1.0)
    src2 = (T2[:3,:3] @ src1.T).T + T2[:3,3]
    log_func(f"  Error: {e2:.4f}m")

    log_func("\n=== Stage 3: Fine (0.3m threshold) ===")
    T3, e3 = run_icp(src2, target, 200, 1e-7, 0.3)
    src3 = (T3[:3,:3] @ src2.T).T + T3[:3,3]
    log_func(f"  Error: {e3:.4f}m")

    log_func("\n=== Stage 4: Ultra-fine (0.1m threshold) ===")
    T4, e4 = run_icp(src3, target, 200, 1e-8, 0.1)
    log_func(f"  Error: {e4:.4f}m")

    # Combine transforms
    T_final = T4 @ T3 @ T2 @ T1

    # Final stats
    src_final = (T_final[:3,:3] @ source.T).T + T_final[:3,3]
    tree = cKDTree(target)
    d, _ = tree.query(src_final, k=1)
    inlier = d < 0.5

    stats = {
        'mean_m': float(np.mean(d[inlier])),
        'median_m': float(np.median(d[inlier])),
        'p90_m': float(np.percentile(d[inlier], 90)),
        'p95_m': float(np.percentile(d[inlier], 95)),
        'inlier_pct': float(100 * np.mean(inlier))
    }

    return T_final, stats


def apply_transform_to_ply(input_path, output_path, T, log_func=print):
    """
    Apply 4x4 transform to PLY file, baking alignment into vertex data.
    
    Transforms:
    - Positions (x, y, z) — vectorized, fast
    - Rotations (rot_0, rot_1, rot_2, rot_3 quaternion) — vectorized
    
    Preserves all other 3DGS properties (colors, scales, spherical harmonics).
    """
    log_func(f"\nBaking transform into PLY...")
    log_func(f"  Reading {Path(input_path).name}...")
    
    plydata = PlyData.read(input_path)
    vertex = plydata['vertex']
    n_points = len(vertex.data)
    
    # Get property names
    prop_names = [p.name for p in vertex.properties]
    log_func(f"  Properties: {len(prop_names)} ({', '.join(prop_names[:5])}...)")
    log_func(f"  Transforming {n_points:,} gaussians...")
    
    # Extract rotation matrix and translation
    R = T[:3, :3]
    t = T[:3, 3]
    
    # --- Transform positions (vectorized) ---
    log_func(f"  [1/3] Transforming positions...")
    x = np.array(vertex['x'], dtype=np.float64)
    y = np.array(vertex['y'], dtype=np.float64)
    z = np.array(vertex['z'], dtype=np.float64)
    positions = np.column_stack([x, y, z])
    new_positions = (R @ positions.T).T + t
    
    # --- Transform rotations (vectorized) ---
    has_rotations = all(f'rot_{i}' in prop_names for i in range(4))
    new_quats = None
    
    if has_rotations:
        log_func(f"  [2/3] Transforming rotations...")
        # 3DGS quaternion order: w, x, y, z (rot_0=w, rot_1=x, rot_2=y, rot_3=z)
        qw = np.array(vertex['rot_0'], dtype=np.float64)
        qx = np.array(vertex['rot_1'], dtype=np.float64)
        qy = np.array(vertex['rot_2'], dtype=np.float64)
        qz = np.array(vertex['rot_3'], dtype=np.float64)
        
        # scipy uses x, y, z, w order
        quats_xyzw = np.column_stack([qx, qy, qz, qw])
        
        # Batch rotation transform
        rot_orig = Rotation.from_quat(quats_xyzw)
        rot_R = Rotation.from_matrix(R)
        rot_new = rot_R * rot_orig
        new_quats_xyzw = rot_new.as_quat()  # Returns Nx4 [x, y, z, w]
    else:
        log_func(f"  [2/3] No rotations to transform")
    
    # --- Build new vertex array ---
    log_func(f"  [3/3] Writing {n_points:,} gaussians...")
    
    # Copy original data
    new_data = np.copy(vertex.data)
    
    # Update positions
    new_data['x'] = new_positions[:, 0].astype(np.float32)
    new_data['y'] = new_positions[:, 1].astype(np.float32)
    new_data['z'] = new_positions[:, 2].astype(np.float32)
    
    # Update rotations if present
    if has_rotations and new_quats_xyzw is not None:
        new_data['rot_0'] = new_quats_xyzw[:, 3].astype(np.float32)  # w
        new_data['rot_1'] = new_quats_xyzw[:, 0].astype(np.float32)  # x
        new_data['rot_2'] = new_quats_xyzw[:, 1].astype(np.float32)  # y
        new_data['rot_3'] = new_quats_xyzw[:, 2].astype(np.float32)  # z
    
    # Create new PLY element
    new_element = PlyElement.describe(new_data, 'vertex')
    
    # Preserve other elements (if any)
    elements = [new_element]
    for elem in plydata.elements:
        if elem.name != 'vertex':
            elements.append(elem)
    
    # Write as binary for efficiency
    log_func(f"  Saving...")
    PlyData(elements, text=False).write(output_path)
    
    file_size = os.path.getsize(output_path) / (1024 * 1024)
    log_func(f"  Saved: {output_path}")
    log_func(f"  Size: {file_size:.1f} MB")
    
    return output_path


def run_alignment(primary_path, secondary_path, output_dir=None, export_aligned_ply=False, log_func=print):
    """
    Main alignment function.
    
    Args:
        primary_path: Reference PLY (stays unchanged)
        secondary_path: PLY to align
        output_dir: Where to save outputs
        export_aligned_ply: If True, also export secondary with transform baked in
        log_func: Logging function
    
    Returns:
        Tuple of (output_data dict, json_path, aligned_ply_path or None)
    """
    if output_dir is None:
        output_dir = Path(secondary_path).parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    log_func("=" * 60)
    log_func("SplatAlign - ICP Alignment for 3DGS")
    log_func("=" * 60)
    
    # Load points
    log_func("\nLOADING FILES")
    log_func("-" * 40)
    primary_pts = load_ply_fast(primary_path, log_func=log_func)
    secondary_pts = load_ply_fast(secondary_path, log_func=log_func)
    
    log_func(f"\nPrimary Z range: [{primary_pts[:,2].min():.1f}, {primary_pts[:,2].max():.1f}]")
    log_func(f"Secondary Z range: [{secondary_pts[:,2].min():.1f}, {secondary_pts[:,2].max():.1f}]")
    
    # Run ICP
    log_func("\nRUNNING ICP ALIGNMENT")
    log_func("-" * 40)
    T, stats = icp_multiscale(secondary_pts, primary_pts, log_func=log_func)
    
    # Results
    log_func("\n" + "=" * 60)
    log_func("RESULTS")
    log_func("=" * 60)
    log_func(f"\nAlignment Quality:")
    log_func(f"  Mean error:   {stats['mean_m']*100:.2f} cm")
    log_func(f"  Median error: {stats['median_m']*100:.2f} cm")
    log_func(f"  90th %ile:    {stats['p90_m']*100:.2f} cm")
    log_func(f"  95th %ile:    {stats['p95_m']*100:.2f} cm")
    log_func(f"  Inliers:      {stats['inlier_pct']:.1f}%")
    
    # Column-major for PlayCanvas/WebGL
    col_major = T.T.flatten().tolist()
    
    log_func("\n4x4 Matrix (row-major):")
    for row in T:
        log_func(f"  [{row[0]:12.8f}, {row[1]:12.8f}, {row[2]:12.8f}, {row[3]:12.8f}]")
    
    log_func("\nColumn-major (PlayCanvas/WebGL):")
    # Format nicely for copy-paste
    log_func(f"  {col_major}")
    
    # Save transform JSON
    primary_name = Path(primary_path).stem
    secondary_name = Path(secondary_path).stem
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    json_filename = f"alignment_{secondary_name}_to_{primary_name}_{timestamp}.json"
    json_path = output_dir / json_filename
    
    output_data = {
        "description": f"Alignment: {secondary_name} -> {primary_name}",
        "created": datetime.now().isoformat(),
        "primary_file": str(primary_path),
        "secondary_file": str(secondary_path),
        "quality": {
            "mean_cm": round(stats['mean_m'] * 100, 2),
            "median_cm": round(stats['median_m'] * 100, 2),
            "p90_cm": round(stats['p90_m'] * 100, 2),
            "p95_cm": round(stats['p95_m'] * 100, 2),
            "inlier_pct": round(stats['inlier_pct'], 1)
        },
        "matrix_row_major": T.tolist(),
        "matrix_column_major_flat": col_major,
    }
    
    with open(json_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    log_func(f"\nSaved transform: {json_path}")
    
    # Optionally export aligned PLY
    aligned_ply_path = None
    if export_aligned_ply:
        aligned_filename = f"{secondary_name}_aligned.ply"
        aligned_ply_path = output_dir / aligned_filename
        apply_transform_to_ply(secondary_path, aligned_ply_path, T, log_func=log_func)
        output_data["aligned_ply"] = str(aligned_ply_path)
    
    return output_data, json_path, aligned_ply_path


# ============================================================
# PyQt6 GUI
# ============================================================

def launch_gui():
    """Launch the PyQt6 GUI."""
    try:
        from PyQt6.QtWidgets import (
            QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
            QLabel, QLineEdit, QPushButton, QTextEdit, QFileDialog,
            QGroupBox, QProgressBar, QMessageBox, QCheckBox
        )
        from PyQt6.QtCore import Qt, QThread, pyqtSignal
        from PyQt6.QtGui import QFont
    except ImportError:
        print("Error: PyQt6 not installed. Install with: pip install PyQt6")
        print("Or use CLI mode: python splat_align.py --cli primary.ply secondary.ply")
        sys.exit(1)
    
    class AlignmentWorker(QThread):
        """Background thread for alignment."""
        log_signal = pyqtSignal(str)
        finished_signal = pyqtSignal(bool, str, object, str)  # Added aligned_ply path
        
        def __init__(self, primary, secondary, output_dir, export_ply):
            super().__init__()
            self.primary = primary
            self.secondary = secondary
            self.output_dir = output_dir
            self.export_ply = export_ply
        
        def run(self):
            try:
                result, json_path, aligned_ply = run_alignment(
                    self.primary,
                    self.secondary,
                    self.output_dir or None,
                    export_aligned_ply=self.export_ply,
                    log_func=lambda msg: self.log_signal.emit(msg)
                )
                self.finished_signal.emit(True, str(json_path), result, str(aligned_ply) if aligned_ply else "")
            except Exception as e:
                import traceback
                self.finished_signal.emit(False, str(e) + "\n" + traceback.format_exc(), None, "")
    
    class SplatAlignWindow(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("SplatAlign - 3DGS Temporal Alignment")
            self.setMinimumSize(800, 700)
            self.last_result = None
            self.worker = None
            self.init_ui()
        
        def init_ui(self):
            central = QWidget()
            self.setCentralWidget(central)
            layout = QVBoxLayout(central)
            layout.setContentsMargins(30, 30, 30, 30)
            layout.setSpacing(15)
            
            # Title
            title = QLabel("SplatAlign")
            title.setFont(QFont("Helvetica", 32, QFont.Weight.Bold))
            title.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(title)
            
            subtitle = QLabel("ICP Alignment for 3D Gaussian Splatting")
            subtitle.setFont(QFont("Helvetica", 14))
            subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
            subtitle.setStyleSheet("color: #666; margin-bottom: 10px;")
            layout.addWidget(subtitle)
            
            # Input files group
            input_group = QGroupBox("Input Files")
            input_layout = QVBoxLayout(input_group)
            
            # Primary file row
            primary_row = QHBoxLayout()
            primary_row.addWidget(QLabel("Primary PLY (reference):"))
            self.primary_edit = QLineEdit()
            self.primary_edit.setPlaceholderText("Select the reference scan...")
            primary_row.addWidget(self.primary_edit, 1)
            primary_btn = QPushButton("Browse...")
            primary_btn.clicked.connect(lambda: self.browse_file(self.primary_edit))
            primary_row.addWidget(primary_btn)
            input_layout.addLayout(primary_row)
            
            # Secondary file row
            secondary_row = QHBoxLayout()
            secondary_row.addWidget(QLabel("Secondary PLY (to align):"))
            self.secondary_edit = QLineEdit()
            self.secondary_edit.setPlaceholderText("Select the scan to align...")
            secondary_row.addWidget(self.secondary_edit, 1)
            secondary_btn = QPushButton("Browse...")
            secondary_btn.clicked.connect(lambda: self.browse_file(self.secondary_edit))
            secondary_row.addWidget(secondary_btn)
            input_layout.addLayout(secondary_row)
            
            layout.addWidget(input_group)
            
            # Output group
            output_group = QGroupBox("Output Options")
            output_layout = QVBoxLayout(output_group)
            
            # Output directory row
            dir_row = QHBoxLayout()
            dir_row.addWidget(QLabel("Output directory:"))
            self.output_edit = QLineEdit()
            self.output_edit.setPlaceholderText("Same as secondary file...")
            dir_row.addWidget(self.output_edit, 1)
            output_btn = QPushButton("Browse...")
            output_btn.clicked.connect(self.browse_output)
            dir_row.addWidget(output_btn)
            output_layout.addLayout(dir_row)
            
            # Export options
            self.export_ply_checkbox = QCheckBox("Export aligned PLY (bake transform into file)")
            self.export_ply_checkbox.setChecked(True)
            self.export_ply_checkbox.setToolTip(
                "Creates a new PLY file with the alignment transform baked in.\n"
                "This file can be loaded directly alongside the primary in any viewer."
            )
            output_layout.addWidget(self.export_ply_checkbox)
            
            # Help text
            help_label = QLabel(
                "💡 The aligned PLY can be loaded directly alongside your primary scan in any 3DGS viewer or editor."
            )
            help_label.setStyleSheet("color: #666; font-size: 11px; margin-top: 5px;")
            help_label.setWordWrap(True)
            output_layout.addWidget(help_label)
            
            layout.addWidget(output_group)
            
            # Buttons
            btn_layout = QHBoxLayout()
            btn_layout.addStretch()
            
            self.run_btn = QPushButton("▶  Run Alignment")
            self.run_btn.setMinimumWidth(150)
            self.run_btn.setMinimumHeight(40)
            self.run_btn.setFont(QFont("Helvetica", 12, QFont.Weight.Bold))
            self.run_btn.clicked.connect(self.run_alignment)
            btn_layout.addWidget(self.run_btn)
            
            self.copy_btn = QPushButton("📋  Copy Transform")
            self.copy_btn.setMinimumWidth(150)
            self.copy_btn.setMinimumHeight(40)
            self.copy_btn.setEnabled(False)
            self.copy_btn.clicked.connect(self.copy_transform)
            btn_layout.addWidget(self.copy_btn)
            
            btn_layout.addStretch()
            layout.addLayout(btn_layout)
            
            # Progress bar
            self.progress = QProgressBar()
            self.progress.setRange(0, 0)  # Indeterminate
            self.progress.setVisible(False)
            layout.addWidget(self.progress)
            
            # Log output
            log_group = QGroupBox("Log")
            log_layout = QVBoxLayout(log_group)
            self.log_text = QTextEdit()
            self.log_text.setReadOnly(True)
            self.log_text.setFont(QFont("Courier", 11))
            self.log_text.setMinimumHeight(180)
            log_layout.addWidget(self.log_text)
            layout.addWidget(log_group, 1)
        
        def browse_file(self, line_edit):
            path, _ = QFileDialog.getOpenFileName(
                self, "Select PLY File", "",
                "PLY Files (*.ply);;All Files (*)"
            )
            if path:
                line_edit.setText(path)
                if not self.output_edit.text():
                    self.output_edit.setText(str(Path(path).parent))
        
        def browse_output(self):
            path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
            if path:
                self.output_edit.setText(path)
        
        def log(self, message):
            self.log_text.append(message)
            self.log_text.verticalScrollBar().setValue(
                self.log_text.verticalScrollBar().maximum()
            )
        
        def run_alignment(self):
            if not self.primary_edit.text():
                QMessageBox.warning(self, "Error", "Please select a primary PLY file.")
                return
            if not self.secondary_edit.text():
                QMessageBox.warning(self, "Error", "Please select a secondary PLY file.")
                return
            
            self.log_text.clear()
            self.run_btn.setEnabled(False)
            self.copy_btn.setEnabled(False)
            self.progress.setVisible(True)
            
            self.worker = AlignmentWorker(
                self.primary_edit.text(),
                self.secondary_edit.text(),
                self.output_edit.text(),
                self.export_ply_checkbox.isChecked()
            )
            self.worker.log_signal.connect(self.log)
            self.worker.finished_signal.connect(self.on_complete)
            self.worker.start()
        
        def on_complete(self, success, message, result, aligned_ply):
            self.progress.setVisible(False)
            self.run_btn.setEnabled(True)
            
            if success:
                self.last_result = result
                self.copy_btn.setEnabled(True)
                self.log(f"\n✅ Alignment complete!")
                
                # Build success message
                msg = f"Alignment complete!\n\nMedian error: {result['quality']['median_cm']:.2f} cm\n\n"
                msg += f"Transform saved to:\n{message}\n"
                if aligned_ply:
                    msg += f"\nAligned PLY saved to:\n{aligned_ply}"
                
                QMessageBox.information(self, "Success", msg)
            else:
                QMessageBox.critical(self, "Error", f"Alignment failed:\n{message}")
                self.log(f"\n❌ Error: {message}")
        
        def copy_transform(self):
            if self.last_result:
                transform = self.last_result.get('matrix_column_major_flat', [])
                clipboard = QApplication.clipboard()
                clipboard.setText(json.dumps(transform))
                self.log("\n📋 Transform copied to clipboard!")
                QMessageBox.information(
                    self, "Copied",
                    "Column-major transform copied to clipboard.\n\nPaste into your config.js transforms object."
                )
    
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    window = SplatAlignWindow()
    window.show()
    sys.exit(app.exec())


# ============================================================
# CLI
# ============================================================

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--cli":
        # Parse CLI args
        if len(sys.argv) >= 4:
            primary = sys.argv[2]
            secondary = sys.argv[3]
            output_dir = None
            export_ply = False
            
            # Check for optional args
            for arg in sys.argv[4:]:
                if arg == "--bake":
                    export_ply = True
                elif not arg.startswith("-"):
                    output_dir = arg
            
            run_alignment(primary, secondary, output_dir, export_aligned_ply=export_ply)
        else:
            print("Usage: splat_align.py --cli PRIMARY.ply SECONDARY.ply [OUTPUT_DIR] [--bake]")
            print("")
            print("Options:")
            print("  --bake    Export aligned PLY with transform baked in")
            sys.exit(1)
    else:
        launch_gui()


if __name__ == "__main__":
    main()

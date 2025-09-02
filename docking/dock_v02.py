#dock_v02.py
import os, stat
import threading
from threading import Semaphore
import queue
import glob
import subprocess
import time
from pathlib import Path
import psutil
import uuid
import gc
import xml.etree.ElementTree as ET
import numpy as np
import shutil  
from threading import Barrier
from config import LIGANDS_DIR, DOCKING_DIR, ANALYSIS_DIR, VINA_DIR, AUTODOCK_GPU_DIR, MACRO_MOL_DIR, DB_PATH, GPU_TYPE, RESULTS_DIR, AUTODOCK_GPU_DIR, NUMWI
from db_manager import DockingDatabaseManager

def _ensure_vina_exec():
    for exe in ("vina", "vina_split"):
        path = os.path.join(VINA_DIR, "bin", exe)
        try:
            st = os.stat(path)
            # if owner-x isn’t already set, add it (resulting in at least 0o755)
            if not (st.st_mode & stat.S_IXUSR):
                os.chmod(path, st.st_mode | stat.S_IXUSR)
        except FileNotFoundError:
            print(f"Warning: {exe} not found at {path}")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
BINARY_PATH = os.path.join(AUTODOCK_GPU_DIR, "bin", "autodock_gpu_" + NUMWI + "wi")
COMPILER_SCRIPT = os.path.join(os.path.dirname(__file__), "autodock-gpu-compiler.sh")

BACKEND = (GPU_TYPE or "CPU").upper()
IS_CUDA = BACKEND == "CUDA"
IS_OPENCL = BACKEND == "OPENCL"
IS_GPU = IS_CUDA or IS_OPENCL

def list_nvidia_gpus():
    """Return a list of GPU indices [0,1,...]. Falls back to [0] if unknown."""
    if shutil.which("nvidia-smi") is None:
        return [0]
    try:
        out = subprocess.check_output(["nvidia-smi", "-L"], text=True)
        idxs = []
        for line in out.strip().splitlines():
            # lines look like: "GPU 0: NVIDIA GeForce RTX 5070 (UUID: ...)"
            if line.startswith("GPU "):
                num = line.split()[1].rstrip(":")
                idxs.append(int(num))
        return idxs or [0]
    except Exception:
        return [0]

GPU_IDS = list_nvidia_gpus()
# BoundedSemaphore protects against accidental double-release
GPU_SLOTS_PER_DEV = int(os.environ.get("GPU_SLOTS_PER_DEV", "1"))
GPU_SEMAPHORES = {i: threading.BoundedSemaphore(GPU_SLOTS_PER_DEV) for i in GPU_IDS}
_GPU_RR_LOCK = threading.Lock()
_GPU_RR = 0

def next_gpu_id():
    """Round-robin GPU selection; safe even if GPU_IDS == [0] fallback."""
    global _GPU_RR
    with _GPU_RR_LOCK:
        if not GPU_IDS:
            # extremely defensive; shouldn't happen because list_nvidia_gpus returns [0] on failure
            return 0
        gid = GPU_IDS[_GPU_RR % len(GPU_IDS)]
        _GPU_RR = (_GPU_RR + 1) % (10_000_000)  # avoid unbounded growth
        return gid

class DockingProcessor:

    def __init__(self):
        self.FILES = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")
        print(f"[INIT] ligands discovered: {len(self.FILES)} in {LIGANDS_DIR}")

        # NOTE: barrier/event were not used; keeping them None avoids accidental waits later
        self.barrier = None
        self.event = None

        # Initialize database manager and queue for async writes
        cpu_count = os.cpu_count() or 1
        cpus_per_vina = 2
        max_parallel = max(1, cpu_count // cpus_per_vina)
        print(f"[INIT] allowing up to {max_parallel} concurrent Vina runs")
        self.vina_sem = Semaphore(max_parallel)

        self.db_manager = DockingDatabaseManager(DB_PATH)
        self.db_queue = queue.Queue()
        self.db_thread = threading.Thread(target=self._db_worker, daemon=True)
        self.db_thread.start()


    def _db_worker(self):
        print("DB worker started")
        buffer = []
        seen_none = False

        while True:
            try:
                print("DB worker waiting for item...")
                item = self.db_queue.get(timeout=10)
                print(f"DB worker received: {item}")

                if item == "INIT":
                    print("Received warm-up INIT. Touching DB...")
                    self.db_manager.insert_bulk([])  # Safe no-op
                    self.db_queue.task_done()
                    continue

                elif item == "SHUTDOWN":
                    print("Shutdown signal received. Flushing and exiting DB thread.")
                    if buffer:
                        self.db_manager.insert_bulk(buffer)
                    self.db_queue.task_done()
                    break
                # Item is valid
                added = 0
                if isinstance(item, list):
                    for rec in item:
                        if isinstance(rec, tuple) and len(rec) == 5:
                            buffer.append(rec); added += 1
                        else:
                            print(f"[DB] skipping bad sub-record: {rec}")
                elif isinstance(item, tuple) and len(item) == 5:
                    buffer.append(item); added = 1
                else:
                    print(f"[DB] skipping malformed record: {item}")

                self.db_queue.task_done()  # Very important!

                if len(buffer) >= 500:
                    try:
                        print(f"[DB] flushing {len(buffer)} records")
                        self.db_manager.insert_bulk(buffer)
                        buffer.clear()
                    except Exception as e:
                        # don't die; log and continue
                        print(f"[DB] insert_bulk failed on batch of {len(buffer)}: {e}")
                        buffer.clear()

            except queue.Empty:
                if buffer:
                    try:
                        print(f"[DB] timeout flush: writing {len(buffer)} records")
                        self.db_manager.insert_bulk(buffer)
                    except Exception as e:
                        print(f"[DB] timeout flush failed on {len(buffer)} records: {e}")
                    buffer.clear()


    def memory_monitor(self, threshold):
        while True:
            available_memory = psutil.virtual_memory().available
            if available_memory > threshold:
                break
            time.sleep(1)

    def process_files(self):
        threads = []
        completed_threads = 0
        lock = threading.Lock()

        def thread_callback():
            nonlocal completed_threads
            with lock:
                completed_threads += 1

        for ligand_file in self.FILES:
            thread = ProcessFileThread(
                [ligand_file],
                self.barrier,
                self.event,
                thread_callback,
                self.db_queue,
                self.vina_sem,
                gpu_id= next_gpu_id() 
            )
            threads.append(thread)

        # 1 KB threshold is effectively always true; keep but log once
        threshold = 1024
        if psutil.virtual_memory().available <= threshold:
            print("[WARN] extremely low available memory reported; waiting...")
            while psutil.virtual_memory().available <= threshold:
                time.sleep(1)

        for thread in threads:
            thread.start()
            print(f"[RUN] launched {len(threads)} worker threads")


        while completed_threads < len(threads):
            time.sleep(1)

        time.sleep(2)  # Give threads time to finish putting into the queue

        self.db_queue.put("INIT")
        self.db_queue.put("SHUTDOWN")
        self.db_queue.join()           # Wait for all DB queue tasks to complete
        self.db_thread.join()          # Wait for the DB thread to fully exit
        self.db_manager.close()        # Only now it's safe to close DB

        gc.collect()
        self.check_memory()

    def check_memory(self):
        bunch = []
        for i in range(0, len(self.FILES), 6):
            current_bunch = self.FILES[i:i+6]
            bunch.extend(current_bunch)
            memory_info = psutil.Process().memory_info()
            print(f"Memory usage: {memory_info.rss / (1024 * 1024)} MB")
            break  # ← prevent barrier from ever being reached


    def run(self):
        self.process_files()


'''
###this is a placeholder for the GPFGenerator class, which is not used in the current version of the script.
    def extract_atom_types_from_pdbqt(self, ligand_file):
        atom_types = set()
        try:
            with open(ligand_file, 'r') as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atom_type = line[77:79].strip()  # AutoDock-style atom type
                        if atom_type:
                            atom_types.add(atom_type)
        except Exception as e:
            print(f"[ERROR] Failed to extract atom types: {e}")
        return " ".join(sorted(atom_types))
'''

class GPFGenerator:
    def __init__(self, macromol_dir=None):
        # default to config’s MACRO_MOL_DIR if not provided
        self.mdir = macromol_dir or MACRO_MOL_DIR

    def calculate_grid_center_and_size(self, receptor_file, padding=10.0):
        ##in the past i was mixing macro molecules and ligands, which is a bad idea
        # This function calculates the geometric center and size of the grid based on the receptor file.
        x_coords, y_coords, z_coords = [], [], []
        with open(receptor_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        x_coords.append(x)
                        y_coords.append(y)
                        z_coords.append(z)
                    except ValueError:
                        continue

        center_x = np.mean(x_coords)
        center_y = np.mean(y_coords)
        center_z = np.mean(z_coords)

        size_x = np.max(x_coords) - np.min(x_coords) + padding
        size_y = np.max(y_coords) - np.min(y_coords) + padding
        size_z = np.max(z_coords) - np.min(z_coords) + padding
        # DEBUG
        print("[RECEPTOR] file:", receptor_file)
        print("[RECEPTOR] bbox min:", (f"{np.min(x_coords):.3f}", f"{np.min(y_coords):.3f}", f"{np.min(z_coords):.3f}"))
        print("[RECEPTOR] bbox max:", (f"{np.max(x_coords):.3f}", f"{np.max(y_coords):.3f}", f"{np.max(z_coords):.3f}"))
        print("[RECEPTOR] center (Å):", (f"{center_x:.3f}", f"{center_y:.3f}", f"{center_z:.3f}"))
        print("[RECEPTOR] raw size+pad (Å):", (f"{size_x:.3f}", f"{size_y:.3f}", f"{size_z:.3f}"))

        return (center_x, center_y, center_z), (size_x, size_y, size_z)


    def _clusters_3d(self, vol, thr):
        nx, ny, nz = vol.shape
        seen = np.zeros(vol.shape, dtype=bool)
        dirs = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
        clusters = []
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    if seen[x,y,z] or vol[x,y,z] < thr:
                        continue
                    # BFS
                    q = [(x,y,z)]
                    seen[x,y,z] = True
                    vox = []
                    s = 0.0
                    while q:
                        a,b,c = q.pop()
                        vox.append((a,b,c))
                        s += vol[a,b,c]
                        for dx,dy,dz in dirs:
                            u,v,w = a+dx, b+dy, c+dz
                            if 0<=u<nx and 0<=v<ny and 0<=w<nz and not seen[u,v,w] and vol[u,v,w] >= thr:
                                seen[u,v,w] = True
                                q.append((u,v,w))
                    clusters.append((s, vox))  # total score & voxels
        # sort by total score descending
        clusters.sort(key=lambda t: t[0], reverse=True)
        return clusters

    def hotspot_center_world_from_cluster(self, voxels, grid_meta):
        (nx, ny, nz), spacing, center = grid_meta
        cx, cy, cz = nx//2, ny//2, nz//2
        # centroid in voxel space
        vx = sum(v[0] for v in voxels)/len(voxels)
        vy = sum(v[1] for v in voxels)/len(voxels)
        vz = sum(v[2] for v in voxels)/len(voxels)
        dx = (vx - cx) * spacing
        dy = (vy - cy) * spacing
        dz = (vz - cz) * spacing
        return (center[0]+dx, center[1]+dy, center[2]+dz)


    def write_gpf_file(self, receptor_file, center, size, output_gpf, spacing=0.375, npts_cap=255):
        ###yup, i realised that generating all maps file in one go is more computationally efficient
        # than parsing all ligands for needed atom types
        # and then generating maps for each ligand separately.
        # This is because AutoDock can generate all maps in one go, and then use them
        # for all ligands separately, as needed. which is much faster than generating maps for each ligand separately.
        DEFAULT_AUTODOCK_ATOM_TYPES = [
            "A",  # aliphatic carbon
            "C",  # aromatic carbon
            "HD", # hydrogen donor
            "N",  # nitrogen
            "NA", # nonpolar nitrogen
            "OA", # oxygen acceptor
            "S",  # sulfur
            "SA", # sulfur acceptor (optional, rarely used)
            "Cl", "Br", "F", "I", # halogens
            "Zn", "Mg", "Ca", "Fe", "Mn"  # metals (optional based on use-case)
        ]
        try:
            # Ensure grid point counts are odd
            npts = [int(dim / spacing) | 1 for dim in size]  # force odd with bitwise OR

            base_name = os.path.splitext(os.path.basename(receptor_file))[0]
            fld_filename = f"{base_name}.maps.fld"
            output_gpf_path = os.path.join(MACRO_MOL_DIR, os.path.basename(output_gpf))

            # DEBUG: print raw inputs
            print("[GPF] receptor:", receptor_file)
            print("[GPF] center (Å):", tuple(f"{c:.3f}" for c in center))
            print("[GPF] requested size (Å):", tuple(f"{s:.3f}" for s in size))
            print("[GPF] initial spacing (Å/grid):", spacing, "  cap:", npts_cap)

            # fit spacing/npts so each dim ≤ 255 and odd
            spacing_fit, npts, box_fit = self._fit_spacing_and_npts(size, spacing, npts_cap=npts_cap)

            # DEBUG: print fitted parameters
            print("[GPF] fitted spacing (Å/grid):", f"{spacing_fit:.6f}")
            print("[GPF] npts (nx,ny,nz):", npts)
            print("[GPF] final box (Å):", tuple(f"{b:.3f}" for b in box_fit))
            print("[GPF] ~grid cells:", npts[0]*npts[1]*npts[2])

            with open(output_gpf_path, 'w') as f:
                f.write(f"npts {npts[0]} {npts[1]} {npts[2]}\n")
                f.write(f"gridfld {fld_filename}\n")
                f.write(f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
                f.write(f"spacing {spacing}\n")
                f.write(f"receptor {os.path.basename(receptor_file)}\n")
                f.write(f"ligand_types {' '.join(DEFAULT_AUTODOCK_ATOM_TYPES)}\n")

                for atom in DEFAULT_AUTODOCK_ATOM_TYPES:
                    f.write(f"map {base_name}.{atom}.map\n")

                f.write(f"elecmap {base_name}.e.map\n")
                f.write(f"dsolvmap {base_name}.d.map\n")
                f.write("dielectric -0.1465\n")


            print(f"[INFO] GPF file written: {output_gpf_path}")
        except Exception as e:
            print(f"[ERROR] Failed to write GPF file: {e}")
            raise
        return output_gpf_path, os.path.join(self.mdir, fld_filename)

    def run_autogrid(self, gpf_path):
        base = Path(gpf_path).stem
        fld_path = Path(self.mdir) / f"{base}.maps.fld"
        if fld_path.exists():
            print(f"[INFO] Skipping AutoGrid: maps already exist ({fld_path})")
            return str(fld_path)

        log_path = gpf_path.replace(".gpf", ".glg")
        if os.path.exists(log_path):
            print(f"[INFO] AutoGrid log already exists: {log_path}")
            return str(fld_path) if fld_path.exists() else None

        try:
            subprocess.run(
                [f"{AUTODOCK_GPU_DIR}/autogrid/autogrid4",
                "-p", os.path.basename(gpf_path),
                "-l", os.path.basename(log_path)],
                check=True, stderr=subprocess.PIPE, text=True,
                cwd=self.mdir, timeout=1200
            )
            print(f"[INFO] AutoGrid finished. Log written to: {log_path}")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] AutoGrid execution failed:\n{e.stderr}")
            raise

        return str(fld_path) if fld_path.exists() else None
    
    def _fit_spacing_and_npts(self, size_xyz, spacing, npts_cap=255):
        """
        Given requested physical box size (Å) and initial spacing (Å/grid),
        increase spacing if needed so that npts per dimension ≤ npts_cap
        and make each npts odd. Returns (spacing, npts_xyz, box_size_xyz).
        """
        import math
        sx, sy, sz = size_xyz
        # raw grid counts
        nx = max(1, int(math.ceil(sx / spacing)))
        ny = max(1, int(math.ceil(sy / spacing)))
        nz = max(1, int(math.ceil(sz / spacing)))

        # if any exceeds cap, scale spacing up
        max_n = max(nx, ny, nz)
        if max_n > npts_cap:
            scale = float(max_n) / float(npts_cap)
            spacing *= scale
            nx = max(1, int(math.ceil(sx / spacing)))
            ny = max(1, int(math.ceil(sy / spacing)))
            nz = max(1, int(math.ceil(sz / spacing)))

        # force odd npts
        if nx % 2 == 0: nx += 1
        if ny % 2 == 0: ny += 1
        if nz % 2 == 0: nz += 1

        # recompute physical box extents the grid will actually cover
        bx = nx * spacing
        by = ny * spacing
        bz = nz * spacing
        return spacing, (nx, ny, nz), (bx, by, bz)



    def ensure_maps(self, receptor_file: str, spacing: float = 0.5, npts_cap: int = 128, include_metals: bool = False):
        base = Path(receptor_file).stem
        fld = Path(self.mdir) / f"{base}.maps.fld"
        gpf = Path(self.mdir) / f"{base}.gpf"
        if fld.exists():
            return str(fld), str(gpf)
        center, size = self.calculate_grid_center_and_size(receptor_file)
        gpf_path, _ = self.write_gpf_file(
            receptor_file, center, size, output_gpf=f"{base}.gpf",
            spacing=spacing, npts_cap=npts_cap, include_metals=include_metals
        )
        self.run_autogrid(gpf_path)
        return str(fld), str(gpf_path)

    #tiny AutoGrid .map reader (ASCII)
    def _read_map_ascii(self, map_path: str):
        """
        Very tolerant reader for AutoGrid ASCII .map:
        - extracts npts (nx,ny,nz) and spacing from the paired .gpf if needed
        - reads all float tokens and reshapes to (nx, ny, nz)
        """
        p = Path(map_path)
        base = p.with_suffix("")  # e.g. 5i6x.OA
        gpf = p.parent / f"{base.name.split('.')[0]}.gpf"

        nx = ny = nz = None
        spacing = None
        origin = None

        # Try to infer grid size from the .map header; fall back to .gpf
        floats = []
        with open(map_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                # Many headers start with non-numeric lines; collect floats anyway
                # We’ll grab numbers; reshape later once we know npts
                for tok in s.split():
                    try:
                        floats.append(float(tok))
                    except ValueError:
                        pass

        # If we can’t detect npts from the map itself, parse the GPF (reliable)
        if gpf.exists():
            with open(gpf, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith("npts"):
                        _, a, b, c = line.split()
                        nx, ny, nz = int(a), int(b), int(c)
                    elif line.startswith("spacing"):
                        _, sp = line.split()
                        spacing = float(sp)
                    elif line.startswith("gridcenter"):
                        _, cx, cy, cz = line.split()
                        origin = (float(cx), float(cy), float(cz))

        if nx is None or ny is None or nz is None:
            raise ValueError(f"Could not determine grid dimensions for {map_path}")

        arr = np.array(floats, dtype=np.float32)
        if arr.size != nx * ny * nz:
            # Some map writers include header numbers; try to trim from the end
            arr = arr[-(nx*ny*nz):]
        vol = arr.reshape((nx, ny, nz), order="C")
        return vol, (nx, ny, nz), spacing, origin

    #simple hotspot score from a few maps 
    def load_hotspot_score(self, receptor_stem: str):
        """
        Build a composite hotspot score from just 3 maps:
        - e.map (electrostatic):   weight -1.0
        - d.map (desolvation):     weight -0.5
        - C.map (hydrophobicity):  weight -0.2
        We z-normalize each map first to balance scales.
        For C (hydrophobic) we clamp to negative only (favorable).
        Returns: score (nx,ny,nz), meta = ((nx,ny,nz), spacing, center)
        """
        base = Path(self.mdir) / receptor_stem
        components = [
            ("e",  -1.0),   # more negative electrostatics → better
            ("d",  -0.5),   # more negative desolvation → better
            ("C",  -0.2),   # negative carbon pockets → better
        ]

        score = None
        vol_ref = None
        for t, w in components:
            p = base.parent / f"{base.name}.{t}.map"
            if not p.exists():
                print(f"[WARN] missing map: {p}")
                continue
            vol, npts, spacing, origin = self._read_map_ascii(str(p))
            # Favorable contributions only where it makes sense
            if t in ("C",):
                vol = np.minimum(vol, 0.0)
            # z-normalize per map to equalize dynamic ranges
            volz = self._z(vol)
            contrib = w * volz
            score = contrib if score is None else (score + contrib)
            vol_ref = (npts, spacing, origin)

        if score is None:
            raise FileNotFoundError("No usable e/d/C maps found to build hotspot score.")
        return score, vol_ref
    # pick top hotspot voxel and convert to world coords
    def hotspot_center_world(self, score, grid_meta):
        (nx, ny, nz), spacing, center = grid_meta
        # score has shape (nx, ny, nz) with center at 'center' in world coords.
        # Convert voxel index of max score to world position:
        idx = np.unravel_index(np.argmax(score), score.shape)
        ix, iy, iz = [int(i) for i in idx]

        # Map index → offset from center (assuming index center at nx//2,…)
        cx, cy, cz = nx // 2, ny // 2, nz // 2
        dx = (ix - cx) * spacing
        dy = (iy - cy) * spacing
        dz = (iz - cz) * spacing
        wx = center[0] + dx
        wy = center[1] + dy
        wz = center[2] + dz
        return (wx, wy, wz)

    # ligand bbox (PDBQT) 
    def ligand_bbox(self, ligand_file: str, pad: float = 2.0):
        xs, ys, zs = [], [], []
        with open(ligand_file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except ValueError:
                        pass
        if not xs:
            # fallback small box
            return (6.0 + pad, 6.0 + pad, 6.0 + pad)
        sx = (max(xs) - min(xs)) + pad
        sy = (max(ys) - min(ys)) + pad
        sz = (max(zs) - min(zs)) + pad
        return (sx, sy, sz)

    # propose docking box from hotspot + ligand size
    def propose_box_for_ligand(self, receptor_file: str, ligand_file: str, fixed_size=(60.0, 60.0, 60.0)):
        """
        1) ensure maps, 2) build 3-map hotspot score, 3) pick peak, 4) use fixed box size.
        Returns: (center_xyz, size_xyz, fld_path)
        """
        fld, gpf = self.ensure_maps(receptor_file)
        stem = Path(receptor_file).stem
        score, meta = self.load_hotspot_score(stem)
        center = self.hotspot_center_world(score, meta)
        size = tuple(float(s) for s in fixed_size)
        # Debug prints so you can verify numbers easily
        (nx, ny, nz), spacing, grid_center = meta
        print(f"[HOTSPOT] grid dims: {nx}x{ny}x{nz}, spacing={spacing}, gridcenter={grid_center}")
        print(f"[HOTSPOT] chosen center: {center}  size: {size}")
        print(f"[HOTSPOT] score stats: min={score.min():.3f}  max={score.max():.3f}  mean={score.mean():.3f}  std={score.std():.3f}")
        print("[BOX] center (Å):", tuple(f"{c:.3f}" for c in center))
        print("[BOX] size   (Å):", tuple(f"{s:.3f}" for s in size))
        print("[BOX] fld:", fld)
        return center, size, fld
    
    def _z(self, arr: np.ndarray) -> np.ndarray:
        m = float(arr.mean())
        s = float(arr.std()) or 1.0
        return (arr - m) / s

    def create_gpf(self, receptor_file, output_gpf="grid_params.gpf", spacing=0.375):
        center, size = self.calculate_grid_center_and_size(receptor_file)
        self.write_gpf_file(receptor_file, center, size, output_gpf, spacing)
        self.run_autogrid(os.path.join(MACRO_MOL_DIR, os.path.basename(output_gpf)))


class ProcessFileThread(threading.Thread):
    
    def __init__(self, bunch, barrier, event, callback, db_queue, vina_sem, gpu_id=0):
        super().__init__()
        self.bunch = bunch
        self.barrier = barrier
        self.event = event
        self.callback = callback
        self.db_queue = db_queue
        self.vina_sem = vina_sem
        self.gpu_id = gpu_id


    # Function to parse ligand's .pdbqt file and extract atomic coordinates
    def calculate_grid_center_and_size(self, ligand_file):
        x_coords, y_coords, z_coords = [], [], []
        
        with open(ligand_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])  # Extract x, y, z coordinates
                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)
        
        # Calculate geometric center (mean of x, y, z coordinates)
        center_x = np.mean(x_coords)
        center_y = np.mean(y_coords)
        center_z = np.mean(z_coords)
        
        # Calculate grid size (based on the bounding box)
        size_x = np.max(x_coords) - np.min(x_coords) + 10  # Add padding to grid size
        size_y = np.max(y_coords) - np.min(y_coords) + 10
        size_z = np.max(z_coords) - np.min(z_coords) + 10
        
        return (center_x, center_y, center_z), (size_x, size_y, size_z)
    
    def parse_vina_output_file(self, filepath, receptor_name, ligand_file):
        # Wait briefly for filesystem
        for _ in range(3):
            if os.path.exists(filepath):
                break
            time.sleep(0.05)
        else:
            print(f"parse_vina_output_file: File not found after retries: {filepath}")
            return []

        results = []
        current_model = None
        current_affinity = None
        current_rmsd_lb = None
        current_rmsd_ub = None
        current_ligand_id = None
        inside_model_block = False

        def finalize():
            nonlocal current_model, current_affinity, current_rmsd_lb, current_rmsd_ub, current_ligand_id, inside_model_block
            if inside_model_block and (current_model is not None) and (current_affinity is not None):
                lid = current_ligand_id or Path(ligand_file).stem or "ligand"
                tag = f"{lid}-Model{current_model}-{receptor_name}"
                results.append((tag, current_affinity, current_rmsd_lb, current_rmsd_ub, ligand_file))
            # reset
            current_model = None
            current_affinity = None
            current_rmsd_lb = None
            current_rmsd_ub = None
            current_ligand_id = None
            inside_model_block = False

        with open(filepath, "r", encoding="utf-8", errors="ignore") as fp:
            for line in fp:
                line = line.strip()

                if line.startswith("MODEL"):
                    # start new model (do not finalize here; we finalize on ENDMDL)
                    inside_model_block = True
                    current_affinity = None
                    current_rmsd_lb = None
                    current_rmsd_ub = None
                    current_ligand_id = None
                    try:
                        current_model = int(line.split()[1])
                    except (IndexError, ValueError):
                        current_model = None

                elif line.startswith("REMARK VINA RESULT:") and inside_model_block:
                    parts = line.split()
                    if len(parts) >= 6 and parts[0:3] == ['REMARK', 'VINA', 'RESULT:']:
                        try:
                            current_affinity = float(parts[3])
                            current_rmsd_lb = float(parts[4])
                            current_rmsd_ub = float(parts[5])
                        except ValueError:
                            pass

                elif line.startswith("REMARK  Name") and inside_model_block:
                    # Some writers use: "REMARK  Name = LIGAND"
                    parts = line.split()
                    if len(parts) >= 4:
                        current_ligand_id = parts[-1]  # safer than [3]

                elif line.startswith("ENDMDL") and inside_model_block:
                    finalize()

        # EOF without ENDMDL: finalize any open block
        if inside_model_block:
            finalize()

        return results


    def parse_adgpu_xml(self, xml_path, receptor_name, ligand_file_fallback):
        """
        Parse AutoDock-GPU XML and return a list of tuples:
        (tag, affinity, rmsd_lb, rmsd_ub, ligand_file)
        tag format: <ligand_id>-Model<run_id>-<receptor_name>
        """
        results = []
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()

            # ligand id from <ligand> element (fallback to provided path)
            ligand_path = root.findtext("ligand") or ligand_file_fallback
            ligand_id = Path(ligand_path).stem

            runs = root.find("runs")
            if runs is None:
                return results

            for run in runs.findall("run"):
                # run id becomes the "Model" number in your schema
                run_id_text = run.get("id")
                try:
                    run_id = int(run_id_text) if run_id_text is not None else 1
                except ValueError:
                    run_id = 1

                # binding energy
                e_txt = run.findtext("free_NRG_binding")
                if e_txt is None:
                    # skip malformed entries
                    continue
                try:
                    affinity = float(e_txt)
                except ValueError:
                    continue

                tag = f"{ligand_id}-Model{run_id}-{receptor_name}"
                # AD-GPU XML doesn't include RMSD values
                results.append((tag, affinity, None, None, ligand_file_fallback))

        except Exception as e:
            print(f"[XML-parse] Failed on {xml_path}: {e}")

        return results


    def run(self):
        # Calculate grid center and size dynamically for each ligand
        macro_mols = glob.glob(f"{MACRO_MOL_DIR}/*.pdbqt")
        if not macro_mols:
           print("No macro molecule found in", MACRO_MOL_DIR)
           return
        macro_mol   = macro_mols[0]
        receptor_name = os.path.basename(macro_mol)

        # 2) ensure maps via GPFGenerator (no rerun if .maps.fld exists)
        gpf = GPFGenerator(MACRO_MOL_DIR)
        fld_file, _gpf_path = gpf.ensure_maps(macro_mol)
        try:
            buffer = []
            BATCH_SIZE = 500
            for ligand_file in self.bunch:
                # ensure output dir exists
                grid_center, grid_size, _ = gpf.propose_box_for_ligand(macro_mol, ligand_file)
                Path(DOCKING_DIR).mkdir(parents=True, exist_ok=True)
                lig_stem  = Path(ligand_file).stem
                out_stem  = Path(DOCKING_DIR) / f"{lig_stem}_{uuid.uuid4().hex[:8]}"
                gpu_pdbqt = f"{out_stem}_out.pdbqt"   # AD-GPU best pose (if --gbest 1)
                xml_out   = f"{out_stem}.xml"     # AD-GPU XML
                vina_out  = f"{out_stem}.pdbqt"       # Vina CPU output target

                # Per-GPU semaphore (one job per GPU at a time)
                sem = GPU_SEMAPHORES.get(self.gpu_id, next(iter(GPU_SEMAPHORES.values())))

                env = os.environ.copy()
                env["OMP_NUM_THREADS"] = "1"
                env["MKL_NUM_THREADS"] = "1"
                env["OPENBLAS_NUM_THREADS"] = "1"
                # If you need hard isolation per thread, uncomment:
                # env["CUDA_VISIBLE_DEVICES"] = str(self.gpu_id)
                acquired_gpu = False
                acquired_cpu = False
                print(ligand_file)

                try:
                    if GPU_TYPE == "NVIDIA" or GPU_TYPE == "CUDA":
                        sem.acquire()
                        acquired_gpu = True
                        result = subprocess.run(
                            [
                                f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_cuda_{NUMWI}wi",
                                "--lfile",   str(Path(ligand_file).resolve()),
                                "--ffile",   str(Path(fld_file).resolve()),
                                "--nrun",    "200",
                                "--gbest",   "5",            # write <resnam>_out.pdbqt (optional but nice)
                                "--xmloutput","1",           # ensure XML is produced
                                "--resnam",  str(out_stem),  # basename; outputs land in DOCKING_DIR
                                "--devnum",  str(int(self.gpu_id) + 1)  # AD-GPU is 1-indexed
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            cwd=MACRO_MOL_DIR,
                            env=env,
                            timeout=100
                        )
                        print(f"[AutoDock-GPU stdout]\n{result.stdout}")
                        print(f"[AutoDock-GPU stderr]\n{result.stderr}")

                        if result.returncode != 0:
                            print(f"Docking failed for: {ligand_file}")
                            print(f"STDOUT:\n{result.stdout.strip()}")
                            print(f"STDERR:\n{result.stderr.strip()}")
                            continue
                        # Ensure output files are created
                    
                    elif GPU_TYPE == "AMD" or GPU_TYPE == "OPENCL":
                        sem.acquire()
                        acquired_gpu = True
                        # --- AutoDock-GPU (XML-first) ---
                        result = subprocess.run(
                            [
                                f"{AUTODOCK_GPU_DIR}/bin/autodock_gpu_ocl_{NUMWI}wi",
                                "--lfile",   str(Path(ligand_file).resolve()),
                                "--ffile",   str(Path(fld_file).resolve()),
                                "--nrun",    "20",
                                "--gbest",   "5",            # write <resnam>_out.pdbqt (optional but nice)
                                "--xmloutput","1",           # ensure XML is produced
                                "--resnam",  str(out_stem),  # basename; outputs land in DOCKING_DIR
                                "--devnum",  str(int(self.gpu_id) + 1)  # AD-GPU is 1-indexed
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            cwd=MACRO_MOL_DIR,
                            env=env,
                            timeout=100
                        )
                        print(f"[AutoDock-GPU stdout]\n{result.stdout}")
                        print(f"[AutoDock-GPU stderr]\n{result.stderr}")

                        if result.returncode != 0:
                            print(f"Docking failed for: {ligand_file}")
                            print(f"STDOUT:\n{result.stdout.strip()}")
                            print(f"STDERR:\n{result.stderr.strip()}")
                            continue
                        # Ensure output files are created

                    elif GPU_TYPE == "CPU" or GPU_TYPE == "VINA": # CPU/Vina fallback, i mean, why not?
                        # --- Vina CPU fallback (PDBQT parsing) ---
                        # CPU/Vina path uses ONLY the CPU limiter
                        self.vina_sem.acquire()
                        acquired_cpu = True
                        print("[VINA] center_xyz:", tuple(f"{v:.3f}" for v in grid_center))
                        print("[VINA] size_xyz:",   tuple(f"{v:.3f}" for v in grid_size))
                        result = subprocess.run(
                            [
                                f"{VINA_DIR}/bin/vina",
                                "--receptor", str(macro_mol),
                                "--ligand",   str(ligand_file),
                                "--center_x", str(grid_center[0]),
                                "--center_y", str(grid_center[1]),
                                "--center_z", str(grid_center[2]),
                                "--size_x",   str(grid_size[0]),
                                "--size_y",   str(grid_size[1]),
                                "--size_z",   str(grid_size[2]),
                                "--cpu",      "2",
                                "--out",      str(vina_out)
                            ],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            timeout=1200
                        )
                        if result.returncode != 0:
                            print(f"Vina failed")   
                            print(f"STDOUT:\n{result.stdout.strip()}")
                            print(f"STDERR:\n{result.stderr.strip()}")
                            continue

                    # Ensure output files are created
                except subprocess.TimeoutExpired:
                    print(f"Docking timed out for: {ligand_file}")
                    sem.release()
                    self.vina_sem.release()
                    continue
                finally:
                    if acquired_gpu:
                        sem.release()
                    if acquired_cpu:
                        self.vina_sem.release()

                print(f"[{threading.current_thread().name}] Docking {ligand_file}")
                if GPU_TYPE == "NVIDIA" or GPU_TYPE == "CUDA":
                    parsed_results = self.parse_adgpu_xml(xml_out, receptor_name, ligand_file)
                elif GPU_TYPE == "AMD" or GPU_TYPE == "OPENCL":
                    parsed_results = self.parse_adgpu_xml(xml_out, receptor_name, ligand_file)
                else:
                    parsed_results = self.parse_vina_output_file(vina_out, receptor_name, ligand_file)

                # parsed_results now unified shape: (tag, affinity, rmsd_lb, rmsd_ub, ligand_file)
                if not parsed_results:
                    print(f"No valid docking data for: {ligand_file}")
                    continue
                buffer.extend(parsed_results)


                if len(buffer) >= BATCH_SIZE:
                    self.db_queue.put(buffer.copy())
                    buffer.clear()

                # memory log
                memory_info = psutil.Process().memory_info()
                print(f"Memory usage: {memory_info.rss / (1024 * 1024):.2f} MB")


                # Set event if memory usage is above the limit 
                if memory_info.rss > 2252800000000: #Set your own limit
                    print(f"Memory usage exceeded the limit.")
                    #self.event.set()
                    #self.barrier.wait()
            if buffer:
                self.db_queue.put(buffer)

        except Exception as e:
            print(e)

        gc.collect()
        # Delete the files
        del self.bunch
        # Notify the main thread that this thread has completed its work
        self.callback()


if __name__ == "__main__":
    start_time = time.time()
    if (GPU_TYPE == "NVIDIA" or GPU_TYPE == "AMD" or GPU_TYPE == "OPENCL" or GPU_TYPE == "CUDA") and (not os.path.exists(os.path.join(MACRO_MOL_DIR, "*.gpf"))):
        gpf_gen = GPFGenerator()
        for receptor in sorted(glob.glob(os.path.join(MACRO_MOL_DIR, "*.pdbqt"))):
            gpf_gen.create_gpf(receptor, output_gpf=f"{Path(receptor).stem}.gpf")
    else:
        print(f"[WARN] Docking engine is set to {GPU_TYPE}, skipping GPF generation.")
    processor = DockingProcessor()
    processor.run()
    print("Process finished --- %s seconds ---" % (time.time() - start_time))
###i can be stupid 
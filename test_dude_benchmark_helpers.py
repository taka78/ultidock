import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

DOCKING_ROOT = Path(__file__).resolve().parent / "docking"
if str(DOCKING_ROOT) not in sys.path:
    sys.path.insert(0, str(DOCKING_ROOT))

import make_grids
import profile_receptors

from benchmarks.dude_docking_benchmark import (
    collect_scores_from_db,
    collect_scores_from_db_by_site,
    reset_sqlite_outputs,
)
from benchmarks.evaluate_cavities import evaluate_sites
from benchmarks.summary_tools import write_method_outputs
from molguard.io.receptor_prep import _parse_excess_bond_residues, sanitize_pdb_for_meeko


def test_sanitize_pdb_for_meeko_preserves_ions_and_remaps_residues(tmp_path: Path) -> None:
    input_pdb = tmp_path / "receptor.pdb"
    output_pdb = tmp_path / "receptor.sanitized.pdb"
    input_pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   DIC A 170      60.292  28.711  -6.506",
                "ATOM      2  HN  DIC A 170      60.114  27.777  -6.166",
                "ATOM      3  CA  ALA A   1      38.390  22.607  -2.054",
                "ATOM      4  ZN   ZN  A 900      -4.452  56.769  81.496",
                "ATOM      5  OW  WAM A 901      -3.709  57.206  77.758",
                "TER",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    sanitize_pdb_for_meeko(input_pdb, output_pdb)
    lines = output_pdb.read_text(encoding="utf-8").splitlines()

    asp_line = next(line for line in lines if line.startswith("ATOM") and " ASP " in line)
    ala_ca_line = next(line for line in lines if line.startswith("ATOM") and " ALA " in line and " CA " in line)
    zn_line = next(line for line in lines if line.startswith("HETATM") and " ZN " in line)

    assert " DIC " not in output_pdb.read_text(encoding="utf-8")
    assert " WAM " not in output_pdb.read_text(encoding="utf-8")
    assert " HN " not in output_pdb.read_text(encoding="utf-8")
    assert asp_line[17:20].strip() == "ASP"
    assert ala_ca_line[76:78].strip() == "C"
    assert zn_line[76:78].strip() == "Zn"


def test_sanitize_pdb_for_meeko_assigns_blank_chain_segments(tmp_path: Path) -> None:
    input_pdb = tmp_path / "receptor.pdb"
    output_pdb = tmp_path / "receptor.sanitized.pdb"
    input_pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA     1       0.000   0.000   0.000",
                "ATOM      2  CA  ALA     1       1.000   0.000   0.000",
                "TER",
                "ATOM      3  N   LYS   443       1.200   0.000   0.000",
                "ATOM      4  CA  LYS   443       2.200   0.000   0.000",
                "TER",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    sanitize_pdb_for_meeko(input_pdb, output_pdb)
    lines = output_pdb.read_text(encoding="utf-8").splitlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    ter_lines = [line for line in lines if line.startswith("TER")]

    assert atom_lines[0][21] == "A"
    assert atom_lines[1][21] == "A"
    assert ter_lines[0][21] == "A"
    assert atom_lines[2][21] == "B"
    assert atom_lines[3][21] == "B"
    assert ter_lines[1][21] == "B"


def test_parse_excess_bond_residues_from_meeko_padding_error() -> None:
    stderr = (
        "matched with excess inter-residue bond(s): A:23\n"
        "matched with excess inter-residue bond(s): d:443\n"
        "RuntimeError: Expected 2 paddings for (A:23, d:443) "
        "with bonds [(8, 18)], but got 0"
    )

    assert _parse_excess_bond_residues(stderr) == ["A:23", "d:443"]


def test_evaluate_sites_handles_legacy_centers_tsv(tmp_path: Path) -> None:
    ligand = tmp_path / "crystal_ligand.mol2"
    centers = tmp_path / "centers.tsv"

    ligand.write_text(
        "\n".join(
            [
                "@<TRIPOS>MOLECULE",
                "ligand",
                "2 1 0 0 0",
                "SMALL",
                "NO_CHARGES",
                "@<TRIPOS>ATOM",
                "1 C1 0.0000 0.0000 0.0000 C.3 1 LIG 0.0",
                "2 C2 2.0000 2.0000 2.0000 C.3 1 LIG 0.0",
                "@<TRIPOS>BOND",
                "1 1 2 1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    centers.write_text(
        "\n".join(
            [
                "receptor\tS1\t1.0000\t1.0000\t1.0000\t20.0\t20.0\t20.0",
                "receptor\tS2\t5.0000\t5.0000\t5.0000\t20.0\t20.0\t20.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    evaluation = evaluate_sites(reference_ligand_path=ligand, centers_tsv_path=centers, hit_threshold_a=4.0)

    assert evaluation["n_sites"] == 2
    assert evaluation["best_site_id"] == "S1"
    assert evaluation["best_site_order"] == 1
    assert evaluation["any_site_success"] is True


def test_collect_scores_from_db_by_site_keeps_site_metrics_separate(tmp_path: Path) -> None:
    db_path = tmp_path / "dock.db"
    ligand_a = (tmp_path / "active__ligA.pdbqt").resolve()
    ligand_b = (tmp_path / "decoy__ligB.pdbqt").resolve()
    ligand_a.write_text("", encoding="utf-8")
    ligand_b.write_text("", encoding="utf-8")

    with sqlite3.connect(db_path) as conn:
        conn.execute(
            """
            CREATE TABLE docking_results (
                docking_file TEXT,
                "binding_affinity (kcal/mol)" REAL,
                binding_site TEXT
            )
            """
        )
        conn.executemany(
            'INSERT INTO docking_results (docking_file, "binding_affinity (kcal/mol)", binding_site) VALUES (?, ?, ?)',
            [
                (str(ligand_a), -7.0, "S1"),
                (str(ligand_a), -8.5, "S2"),
                (str(ligand_b), -5.0, "S1"),
                (str(ligand_b), -4.0, "S2"),
            ],
        )

    label_map = {ligand_a: 1, ligand_b: 0}
    combined = collect_scores_from_db(db_path, label_map)
    per_site = collect_scores_from_db_by_site(db_path, label_map)

    assert len(combined) == 2
    assert {row.best_binding_site for row in combined} == {"S1", "S2"}
    assert sorted(per_site) == ["S1", "S2"]
    assert [row.best_affinity_kcal_mol for row in per_site["S1"]] == [-7.0, -5.0]
    assert [row.best_affinity_kcal_mol for row in per_site["S2"]] == [-8.5, -4.0]


def test_write_method_outputs_tracks_completed_failed_and_pending_targets(tmp_path: Path) -> None:
    dataset_root = tmp_path / "datasets"
    dataset_root.mkdir()
    output_dir = tmp_path / "results"
    method_label = "ultidock_cpu"

    ace_dir = output_dir / "ace" / method_label
    ace_dir.mkdir(parents=True)
    (ace_dir / "metrics.json").write_text(
        json.dumps(
            {
                "target": "ace",
                "roc_auc": 0.75,
                "ef1": 4.0,
                "ef5": 2.0,
                "runner_elapsed_s": 12.0,
                "site_evaluation": {
                    "best_distance_a": 1.8,
                    "best_site_order": 1,
                    "any_site_success": True,
                },
            }
        ),
        encoding="utf-8",
    )

    ampc_dir = output_dir / "ampc" / method_label
    ampc_dir.mkdir(parents=True)
    (ampc_dir / "failure.json").write_text(
        json.dumps(
            {
                "target": "ampc",
                "status": "timeout",
                "reason": "timeout",
                "elapsed_s": 3.0,
            }
        ),
        encoding="utf-8",
    )

    refreshed = write_method_outputs(
        dataset_root=dataset_root,
        output_dir=output_dir,
        method_label=method_label,
        target_names=["ace", "ampc", "hivpr"],
        mode="cpu",
        center_source="auto",
        seed=42,
        invocation_wall_time_s=15.0,
    )

    summary = refreshed["summary"]
    status = refreshed["status"]
    runtime_summary = refreshed["runtime_summary"]
    failure_summary = refreshed["failure_summary"]

    assert summary["n_targets"] == 3
    assert summary["n_completed"] == 1
    assert summary["n_failed"] == 1
    assert summary["n_pending"] == 1
    assert summary["aggregate"]["roc_auc_mean"] == 0.75
    assert summary["site_validation"]["mean_best_distance_a"] == 1.8
    assert summary["last_invocation_wall_time_s"] == 15.0

    assert status["completed_targets"] == ["ace"]
    assert status["failed_targets"] == ["ampc"]
    assert status["pending_targets"] == ["hivpr"]

    assert runtime_summary["n_targets_with_runtime"] == 2
    assert runtime_summary["total_observed_runtime_s"] == 15.0
    assert runtime_summary["total_completed_runtime_s"] == 12.0
    assert runtime_summary["total_failed_runtime_s"] == 3.0

    assert failure_summary["n_failed"] == 1
    assert failure_summary["by_status"] == {"timeout": 1}
    assert failure_summary["failures"][0]["target"] == "ampc"

    for path in refreshed["paths"].values():
        assert path.exists()


def test_reset_sqlite_outputs_removes_db_and_sidecars(tmp_path: Path) -> None:
    db_path = tmp_path / "ultidock_results.db"
    wal_path = tmp_path / "ultidock_results.db-wal"
    shm_path = tmp_path / "ultidock_results.db-shm"

    for path in (db_path, wal_path, shm_path):
        path.write_text("x", encoding="utf-8")

    reset_sqlite_outputs(db_path)

    assert not db_path.exists()
    assert not wal_path.exists()
    assert not shm_path.exists()


def test_hotspot_box_sizing_uses_side_length_semantics() -> None:
    assert make_grids._npts_for_box_side(31.2, 0.375) == 85
    assert make_grids._npts_for_size(31.2, 0.375) == 85
    assert make_grids._npts_for_extent(31.2, 0.375, 129.0) == 85


def test_maps_pocket_shell_bounds_follow_profiled_r_min() -> None:
    assert make_grids._maps_pocket_shell_bounds(3.0) == (1.5, 12.0)
    assert make_grids._maps_pocket_shell_bounds(5.0) == (2.5, 12.0)
    assert make_grids._maps_pocket_shell_bounds(1.0) == (0.75, 12.0)


def test_detect_maps_hotspots_does_not_double_box_side() -> None:
    shape = (31, 31, 31)
    hotspot_slice = slice(13, 18)

    c_map = np.zeros(shape, dtype=np.float32)
    e_map = np.zeros(shape, dtype=np.float32)
    d_map = np.zeros(shape, dtype=np.float32)

    c_map[hotspot_slice, hotspot_slice, hotspot_slice] = -5.0
    e_map[hotspot_slice, hotspot_slice, hotspot_slice] = -20.0
    d_map[6:25, 6:25, 6:25] = 1.0

    sites = make_grids.detect_maps_hotspots(
        c_map,
        e_map,
        d_map,
        origin=(0.0, 0.0, 0.0),
        spacing=0.375,
        max_sites=1,
        box_side_A=31.2,
    )

    assert len(sites) == 1
    assert (sites[0]["nx"], sites[0]["ny"], sites[0]["nz"]) == (85, 85, 85)


def test_sparse_whole_protein_maps_are_not_marked_degenerate() -> None:
    c_map = np.zeros((20, 20, 20), dtype=np.float32)
    e_axis = np.linspace(-0.08, -0.02, 20, dtype=np.float32)
    e_map = np.broadcast_to(e_axis[None, None, :], (20, 20, 20)).copy()
    d_map = np.zeros((20, 20, 20), dtype=np.float32)

    c_map[8:12, 8:12, 8:12] = -0.2
    d_map[7:13, 7:13, 7:13] = 0.5

    assert make_grids._maps_are_degenerate(c_map, e_map, d_map) is False


def test_flat_electrostatics_are_marked_degenerate() -> None:
    c_map = np.zeros((10, 10, 10), dtype=np.float32)
    e_map = np.zeros((10, 10, 10), dtype=np.float32)
    d_map = np.zeros((10, 10, 10), dtype=np.float32)

    assert make_grids._maps_are_degenerate(c_map, e_map, d_map) is True


def test_select_r_min_cavity_uses_significant_peak_lower_tail() -> None:
    stats = {
        "p25": 2.99,
        "significant_peak_p10": 3.76,
    }

    assert profile_receptors._select_r_min_cavity(stats) == 3.76


def test_select_r_min_cavity_falls_back_to_p25_when_needed() -> None:
    stats = {
        "p25": 2.45,
    }

    assert profile_receptors._select_r_min_cavity(stats) == 2.45


def test_build_hybrid_sites_combines_internal_and_surface_candidates() -> None:
    internal_sites = [
        {"site_id": "S1", "cx": 0.0, "cy": 0.0, "cz": 0.0, "nx": 81, "ny": 81, "nz": 81, "spacing": 0.375, "r_peak": 4.0, "F": 0.9},
        {"site_id": "S2", "cx": 20.0, "cy": 0.0, "cz": 0.0, "nx": 81, "ny": 81, "nz": 81, "spacing": 0.375, "r_peak": 3.5, "F": 0.4},
    ]
    surface_sites = [
        {"site_id": "S1", "cx": 1.5, "cy": 0.0, "cz": 0.0, "nx": 85, "ny": 85, "nz": 85, "spacing": 0.375, "r_peak": 2.2, "F": 0.8},
        {"site_id": "S2", "cx": 40.0, "cy": 0.0, "cz": 0.0, "nx": 85, "ny": 85, "nz": 85, "spacing": 0.375, "r_peak": 2.0, "F": 0.2},
    ]

    hybrid_sites = make_grids._build_hybrid_sites(
        internal_sites,
        surface_sites,
        min_sep_A=4.0,
        max_sites=3,
    )

    assert hybrid_sites
    assert hybrid_sites[0]["family"] == "hybrid"
    assert hybrid_sites[0]["selection_score"] > 0.0
    assert abs(hybrid_sites[0]["cx"]) < 1.5


def test_receptor_search_picks_one_site_per_family_when_available() -> None:
    internal = [
        {"site_id": "I1", "cx": 0.0, "cy": 0.0, "cz": 0.0, "selection_score": 0.9, "F": 0.9, "family": "internal"},
    ]
    surface = [
        {"site_id": "S1", "cx": 15.0, "cy": 0.0, "cz": 0.0, "selection_score": 0.8, "F": 0.8, "family": "surface"},
    ]
    hybrid = [
        {"site_id": "H1", "cx": 30.0, "cy": 0.0, "cz": 0.0, "selection_score": 0.7, "F": 0.7, "family": "hybrid"},
    ]

    selected = make_grids._assemble_receptor_search_sites(
        internal,
        surface,
        hybrid,
        min_sep_A=10.0,
        target_count=3,
    )

    assert [site["site_id"] for site in selected] == ["S1", "S2", "S3"]
    assert [site["family"] for site in selected] == ["internal", "surface", "hybrid"]


def test_receptor_search_preserves_central_surface_candidate() -> None:
    internal = [
        {"site_id": "I1", "cx": 0.0, "cy": 0.0, "cz": 0.0, "selection_score": 0.9, "F": 0.9, "family": "internal"},
    ]
    surface = [
        {"site_id": "S1", "cx": 20.0, "cy": 0.0, "cz": 0.0, "selection_score": 1.0, "F": 1.0, "family": "surface", "center_closeness": 0.1},
        {"site_id": "S2", "cx": 0.0, "cy": 20.0, "cz": 0.0, "selection_score": 0.2, "F": 0.2, "family": "surface", "center_closeness": 0.95},
    ]
    hybrid = [
        {"site_id": "H1", "cx": 40.0, "cy": 0.0, "cz": 0.0, "selection_score": 0.8, "F": 0.8, "family": "hybrid"},
    ]

    selected = make_grids._assemble_receptor_search_sites(
        internal,
        surface,
        hybrid,
        min_sep_A=10.0,
        target_count=4,
    )

    assert [site["family"] for site in selected[:2]] == ["surface", "surface"]
    assert any(site.get("center_closeness") == 0.95 for site in selected)


def test_receptor_search_preserves_low_score_central_hybrid_candidate() -> None:
    surface = [
        {"site_id": "S1", "cx": 0.0, "cy": 0.0, "cz": 0.0, "selection_score": 1.0, "F": 1.0, "family": "surface"},
    ]
    hybrid = [
        {"site_id": "H1", "cx": 20.0, "cy": 0.0, "cz": 0.0, "selection_score": 1.0, "F": 1.0, "family": "hybrid", "center_closeness": 0.1},
        {"site_id": "H2", "cx": 40.0, "cy": 0.0, "cz": 0.0, "selection_score": 0.8, "F": 0.8, "family": "hybrid", "center_closeness": 0.2},
        {"site_id": "H3", "cx": 0.0, "cy": 20.0, "cz": 0.0, "selection_score": 0.1, "F": 0.1, "family": "hybrid", "center_closeness": 0.95},
    ]

    selected = make_grids._assemble_receptor_search_sites(
        [],
        surface,
        hybrid,
        min_sep_A=10.0,
        target_count=4,
    )

    assert any(site["family"] == "hybrid" and site.get("center_closeness") == 0.95 for site in selected)


def test_select_internal_search_r_min_relaxes_when_profiled_threshold_is_too_strict() -> None:
    relaxed = make_grids._select_internal_search_r_min(
        requested_r_min_A=3.76,
        peak_radii_A=[1.05, 1.18, 1.24, 1.31, 1.35],
    )

    assert relaxed == 1.33


def test_select_internal_search_r_min_keeps_requested_threshold_when_supported() -> None:
    relaxed = make_grids._select_internal_search_r_min(
        requested_r_min_A=3.20,
        peak_radii_A=[2.9, 3.4, 3.8, 4.1, 4.4],
    )

    assert relaxed == 3.2

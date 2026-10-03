from __future__ import annotations

import csv
import importlib
import sys
import types
from pathlib import Path

import numpy as np
import pytest


def _load_make_grids(monkeypatch):
    config = types.ModuleType("config")
    values = {
        "AUTODOCK_GPU_DIR": "",
        "DOCKING_DIR": "",
        "RESULTS_DIR": "",
        "NUMWI": 1,
        "LIGANDS_DIR": "",
        "CENTERS_TSV": "centers.tsv",
        "GRID_MODE": "blind",
        "GRID_SPACING": 0.375,
        "GRID_MARGIN": 5.0,
        "GRID_CAP": 150.0,
        "R_MIN_CAVITY_A": None,
        "HOTSPOT_BOX_ANGLE": 35.0,
        "HOTSPOT_NMS_MINSEP_A": 14.0,
        "SURFACE_SHELL__MIN_A": 2.0,
        "SURFACE_SHELL__MAX_A": 20.0,
        "SURFACE_NMS_MINSEP_A": 5.0,
        "MAX_CENTER_DIST_A": 10.0,
        "CONTACT_SHELL_A": 4.0,
        "MIN_SURFACE_FRAC": 0.01,
        "AUTOSITES": 6,
    }
    for key, value in values.items():
        setattr(config, key, value)
    monkeypatch.setitem(sys.modules, "config", config)
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "docking"))
    sys.modules.pop("make_grids", None)
    return importlib.import_module("make_grids")


def _write_map(path: Path, values: np.ndarray) -> None:
    path.write_text(
        "\n".join(
            [
                "GRID_PARAMETER_FILE grid.gpf",
                "GRID_DATA_FILE receptor.maps.fld",
                "MACROMOLECULE receptor.pdbqt",
                "SPACING 0.500",
                "NELEMENTS 4 6 8",
                "CENTER 10.000 20.000 30.000",
                *(f"{float(value):.6f}" for value in values.ravel(order="F")),
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_known_site_box_applies_to_the_selected_receptor(monkeypatch, tmp_path):
    from cli.ultidock import _write_known_site_tsv

    grids = _load_make_grids(monkeypatch)
    centers = tmp_path / "sites.tsv"
    _write_known_site_tsv(centers, (-36.106, -20.758, 4.897), 20.0)

    rows = grids._parse_centers_tsv(centers, receptor_key="5i6x_edited")

    assert len(rows) == 1
    assert rows[0]["site_id"] == "S1"
    assert rows[0]["center"] == (-36.106, -20.758, 4.897)


def test_predicted_boxes_still_require_the_matching_receptor(monkeypatch, tmp_path):
    grids = _load_make_grids(monkeypatch)
    centers = tmp_path / "sites.tsv"
    centers.write_text(
        "# meta policy=provided\n"
        "known_site\tS1\t1\t2\t3\t95\t95\t95\t0.375\n"
    )

    assert grids._parse_centers_tsv(centers, receptor_key="5i6x_edited") == []


def _pdbqt_atom_line(serial: int, atom_name: str, atom_type: str, x: float = 0.0) -> str:
    line = (
        f"ATOM  {serial:5d} {atom_name:<4s} ALA A{serial:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}"
    )
    return line.ljust(77) + f"{atom_type:>2s}"


def _write_receptor_pdbqt(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                _pdbqt_atom_line(1, "C1", "C", 0.0),
                _pdbqt_atom_line(2, "O1", "OA", 1.0),
                _pdbqt_atom_line(3, "H1", "HD", 2.0),
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_fld(tmp_path: Path, nelements=(4, 6, 8), dims=(5, 7, 9)) -> Path:
    fld = tmp_path / "receptor.maps.fld"
    fld.write_text(
        "\n".join(
            [
                "# AVS field file",
                "#SPACING 0.500",
                f"#NELEMENTS {nelements[0]} {nelements[1]} {nelements[2]}",
                "#CENTER 10.000 20.000 30.000",
                "ndim=3",
                f"dim1={dims[0]}",
                f"dim2={dims[1]}",
                f"dim3={dims[2]}",
                "veclen=3",
                "label=C-affinity",
                "label=Electrostatics",
                "label=Desolvation",
                "variable 1 file=receptor.C.map filetype=ascii skip=6",
                "variable 2 file=receptor.e.map filetype=ascii skip=6",
                "variable 3 file=receptor.d.map filetype=ascii skip=6",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return fld


def _write_fld_without_desolvation_label(tmp_path: Path) -> Path:
    fld = tmp_path / "receptor.maps.fld"
    fld.write_text(
        "\n".join(
            [
                "# AVS field file",
                "#SPACING 0.500",
                "#NELEMENTS 4 6 8",
                "#CENTER 10.000 20.000 30.000",
                "ndim=3",
                "dim1=5",
                "dim2=7",
                "dim3=9",
                "veclen=3",
                "label=C-affinity",
                "label=Electrostatics",
                "label=HD-affinity",
                "variable 1 file=receptor.C.map filetype=ascii skip=6",
                "variable 2 file=receptor.e.map filetype=ascii skip=6",
                "variable 3 file=receptor.HD.map filetype=ascii skip=6",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return fld


def _write_full_ad4_fld(tmp_path: Path, make_grids) -> Path:
    fld = tmp_path / "receptor.maps.fld"
    lines = [
        "# AVS field file",
        "#SPACING 0.500",
        "#NELEMENTS 4 6 8",
        "#CENTER 10.000 20.000 30.000",
        "ndim=3",
        "dim1=5",
        "dim2=7",
        "dim3=9",
        f"veclen={len(make_grids._AD4_TYPES) + 2}",
    ]
    lines.extend(f"label={atom_type}-affinity" for atom_type in make_grids._AD4_TYPES)
    lines.extend(["label=Electrostatics", "label=Desolvation"])
    for index, atom_type in enumerate(make_grids._AD4_TYPES, start=1):
        lines.append(f"variable {index} file=receptor.{atom_type}.map filetype=ascii skip=6")
    lines.append(f"variable {len(make_grids._AD4_TYPES) + 1} file=receptor.e.map filetype=ascii skip=6")
    lines.append(f"variable {len(make_grids._AD4_TYPES) + 2} file=receptor.d.map filetype=ascii skip=6")
    fld.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return fld


def test_fld_nelements_shape_origin_and_round_trip(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    shape = (5, 7, 9)
    values = np.zeros(shape, dtype=np.float32)
    _write_map(tmp_path / "receptor.C.map", values)
    _write_map(tmp_path / "receptor.e.map", values)
    _write_map(tmp_path / "receptor.d.map", values)
    fld = _write_fld(tmp_path)

    meta = make_grids.load_fld_or_map_meta(fld)

    assert meta["shape"] == shape
    assert meta["origin"] == pytest.approx((9.0, 18.5, 28.0))
    assert meta["affinity_types"] == ("C",)
    ijk = (4, 6, 8)
    xyz = make_grids.voxel_to_world(ijk, meta["origin"], meta["spacing"])
    round_trip = tuple(round((xyz[i] - meta["origin"][i]) / meta["spacing"]) for i in range(3))
    assert round_trip == ijk


def test_load_map_ascii_requires_exact_value_count(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    shape = (5, 7, 9)
    values = np.zeros(shape, dtype=np.float32).ravel(order="F")[:-1]
    _write_map(tmp_path / "short.C.map", values)

    with pytest.raises(ValueError, match="expected 315 values"):
        make_grids.load_map_ascii(tmp_path / "short.C.map", shape)


def test_known_synthetic_maximum_world_coordinate(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    shape = (5, 7, 9)
    values = np.zeros(shape, dtype=np.float32)
    expected_ijk = (2, 3, 4)
    values[expected_ijk] = 42.0
    _write_map(tmp_path / "receptor.C.map", values)
    _write_map(tmp_path / "receptor.e.map", np.zeros(shape, dtype=np.float32))
    _write_map(tmp_path / "receptor.d.map", np.zeros(shape, dtype=np.float32))
    meta = make_grids.load_fld_or_map_meta(_write_fld(tmp_path))

    loaded = make_grids.load_map_ascii(meta["map_paths"]["C"], meta["shape"])
    observed_ijk = tuple(int(i) for i in np.unravel_index(np.argmax(loaded), loaded.shape))
    observed_xyz = make_grids.voxel_to_world(observed_ijk, meta["origin"], meta["spacing"])

    assert observed_ijk == expected_ijk
    assert observed_xyz == pytest.approx((10.0, 20.0, 30.0))


def test_common_physics_ranking_overrides_legacy_role_prior(monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    score_field = np.zeros((5, 5, 5), dtype=np.float32)
    edt_grid = np.zeros_like(score_field)
    score_field[1, 1, 1] = 0.2
    score_field[2, 2, 2] = 0.8
    edt_grid[1, 1, 1] = 2.0
    edt_grid[2, 2, 2] = 3.0
    sites = [
        {
            "site_id": "S1",
            "cx": 1.0,
            "cy": 1.0,
            "cz": 1.0,
            "F": 0.5,
            "r_peak": 2.0,
            "portfolio_role": "surface_focus",
        },
        {
            "site_id": "S2",
            "cx": 2.0,
            "cy": 2.0,
            "cz": 2.0,
            "F": 0.5,
            "r_peak": 3.0,
            "portfolio_role": "core_reserve",
        },
    ]

    legacy = make_grids._assign_receptor_search_fitness_scores(sites)
    assert legacy[0]["F"] > legacy[1]["F"]

    ranked = make_grids._assign_common_physics_ranking_scores(
        legacy,
        score_context={
            "origin": (0.0, 0.0, 0.0),
            "spacing": 1.0,
            "score_field": score_field,
            "score_max": 0.8,
            "edt_grid": edt_grid,
            "pocket_min_A": 1.0,
            "pocket_max_A": 4.0,
        },
    )

    assert ranked[0]["F"] == pytest.approx(0.25)
    assert ranked[1]["F"] == pytest.approx(1.0)
    assert ranked[0]["legacy_F"] == pytest.approx(legacy[0]["F"])
    assert ranked[1]["legacy_F"] == pytest.approx(legacy[1]["F"])
    assert ranked[0]["portfolio_role"] == "surface_focus"
    assert ranked[1]["portfolio_role"] == "core_reserve"
    assert (ranked[0]["cx"], ranked[0]["cy"], ranked[0]["cz"]) == (1.0, 1.0, 1.0)
    assert (ranked[1]["cx"], ranked[1]["cy"], ranked[1]["cz"]) == (2.0, 2.0, 2.0)


def test_common_physics_ranking_uses_uniform_edt_fallback(monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    edt_grid = np.zeros((5, 5, 5), dtype=np.float32)
    edt_grid[1, 1, 1] = 2.0
    edt_grid[2, 2, 2] = 4.0
    sites = [
        {"site_id": "S1", "cx": 1.0, "cy": 1.0, "cz": 1.0, "F": 0.9},
        {"site_id": "S2", "cx": 2.0, "cy": 2.0, "cz": 2.0, "F": 0.1},
    ]

    ranked = make_grids._assign_common_physics_ranking_scores(
        sites,
        score_context={
            "origin": (0.0, 0.0, 0.0),
            "spacing": 1.0,
            "score_field": None,
            "score_max": 0.0,
            "edt_grid": edt_grid,
            "pocket_min_A": 1.0,
            "pocket_max_A": 5.0,
        },
    )

    assert ranked[0]["F"] == pytest.approx(0.5)
    assert ranked[1]["F"] == pytest.approx(1.0)
    assert {site["ranking_basis"] for site in ranked} == {"EDT_fallback"}


def test_hotspot_detection_can_return_common_ranking_context(monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    shape = (31, 31, 31)
    c_map = np.zeros(shape, dtype=np.float32)
    e_map = np.zeros(shape, dtype=np.float32)
    d_map = np.zeros(shape, dtype=np.float32)
    c_map[13:18, 13:18, 13:18] = -5.0
    e_map[13:18, 13:18, 13:18] = -20.0
    d_map[6:25, 6:25, 6:25] = 1.0

    sites, context = make_grids.detect_maps_hotspots(
        c_map,
        e_map,
        d_map,
        origin=(0.0, 0.0, 0.0),
        spacing=0.375,
        max_sites=1,
        box_side_A=31.2,
        return_score_context=True,
    )

    assert len(sites) == 1
    assert context["score_field"].shape == shape
    assert context["score_max"] > 0.0
    assert context["origin"] == (0.0, 0.0, 0.0)
    assert context["spacing"] == pytest.approx(0.375)


def test_desolvation_fallback_does_not_select_hd_affinity_map(
    tmp_path: Path, monkeypatch
) -> None:
    make_grids = _load_make_grids(monkeypatch)
    shape = (5, 7, 9)
    values = np.zeros(shape, dtype=np.float32)
    _write_map(tmp_path / "receptor.C.map", values)
    _write_map(tmp_path / "receptor.e.map", values)
    _write_map(tmp_path / "receptor.HD.map", np.ones(shape, dtype=np.float32))
    _write_map(tmp_path / "receptor.d.map", np.full(shape, 2.0, dtype=np.float32))

    meta = make_grids.load_fld_or_map_meta(_write_fld_without_desolvation_label(tmp_path))

    assert meta["map_paths"]["D"].name == "receptor.d.map"


def test_fld_dimension_mismatch_fails(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    _write_fld(tmp_path, nelements=(4, 6, 8), dims=(4, 6, 8))

    with pytest.raises(ValueError, match="NELEMENTS\\+1"):
        make_grids.load_fld_or_map_meta(tmp_path / "receptor.maps.fld")


def test_cav_emps_rejects_c_only_map_generation(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)

    with pytest.raises(ValueError, match="CaV-EMPS requires full AD4 ligand_types"):
        make_grids.autogenerate_centers_tsv(
            receptor_pdbqt=str(tmp_path / "receptor.pdbqt"),
            out_root=str(tmp_path / "maps"),
            centers_tsv_path=str(tmp_path / "centers.tsv"),
            mode="receptor_search",
            map_types=("C",),
        )


def test_internal_mode_rejects_c_only_map_generation(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)

    with pytest.raises(ValueError, match="CaV-EMPS requires full AD4 ligand_types"):
        make_grids.autogenerate_centers_tsv(
            receptor_pdbqt=str(tmp_path / "receptor.pdbqt"),
            out_root=str(tmp_path / "maps"),
            centers_tsv_path=str(tmp_path / "centers.tsv"),
            mode="internal",
            map_types=("C",),
        )


def test_cached_fld_affinity_type_mismatch_fails_before_reuse(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    maps_dir = tmp_path / "maps"
    maps_dir.mkdir()
    shape = (5, 7, 9)
    values = np.zeros(shape, dtype=np.float32)
    _write_map(maps_dir / "receptor.C.map", values)
    _write_map(maps_dir / "receptor.e.map", values)
    _write_map(maps_dir / "receptor.d.map", values)
    _write_fld(maps_dir)

    with pytest.raises(RuntimeError, match="Incompatible cached FLD"):
        make_grids.autogenerate_centers_tsv(
            receptor_pdbqt=str(tmp_path / "receptor.pdbqt"),
            out_root=str(maps_dir),
            centers_tsv_path=str(tmp_path / "centers.tsv"),
            mode="receptor_search",
        )


def test_full_ad4_fld_affinity_types_are_parsed(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    shape = (5, 7, 9)
    values = np.zeros(shape, dtype=np.float32)
    for atom_type in make_grids._AD4_TYPES:
        _write_map(tmp_path / f"receptor.{atom_type}.map", values)
    _write_map(tmp_path / "receptor.e.map", values)
    _write_map(tmp_path / "receptor.d.map", values)
    meta = make_grids.load_fld_or_map_meta(_write_full_ad4_fld(tmp_path, make_grids))

    assert meta["affinity_types"] == tuple(make_grids._AD4_TYPES)


def test_receptor_types_are_extracted_but_ligand_types_remain_full_ad4(
    tmp_path: Path, monkeypatch
) -> None:
    make_grids = _load_make_grids(monkeypatch)
    receptor = tmp_path / "receptor.pdbqt"
    _write_receptor_pdbqt(receptor)
    gpf = tmp_path / "grid.gpf"

    make_grids._write_site_gpf(gpf, receptor, center=(0.0, 0.0, 0.0), npts=(5, 5, 5), spacing=0.5)
    text = gpf.read_text(encoding="utf-8")

    assert "receptor_types C OA HD\n" in text
    assert "ligand_types " + " ".join(make_grids._AD4_TYPES) + "\n" in text


def test_unknown_receptor_atom_type_fails(tmp_path: Path, monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    receptor = tmp_path / "bad_receptor.pdbqt"
    receptor.write_text(_pdbqt_atom_line(1, "X1", "XX") + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="unsupported receptor atom type"):
        make_grids.receptor_atom_types(receptor)


def test_surface_refinement_stays_in_seed_component(monkeypatch) -> None:
    make_grids = _load_make_grids(monkeypatch)
    shape = (21, 21, 21)
    score = np.zeros(shape, dtype=np.float32)
    edt = np.ones(shape, dtype=np.float32)
    pocket = np.zeros(shape, dtype=bool)
    border = np.ones(shape, dtype=bool)

    pocket[4:7, 9:12, 9:12] = True
    score[4:7, 9:12, 9:12] = 1.0
    pocket[14:17, 9:12, 9:12] = True
    score[14:17, 9:12, 9:12] = 1.0

    refined, diagnostics = make_grids._refine_surface_hotspot_ijk(
        np.array([5, 10, 10], dtype=np.int32),
        score=score,
        edt_grid=edt,
        pocket_mask=pocket,
        border_mask=border,
        spacing=1.0,
        raw_score=1.0,
        radius_A=20.0,
        return_diagnostics=True,
    )

    assert diagnostics["refinement_component_count"] == 2
    assert diagnostics["refinement_component_size"] == 27
    assert refined[0] < 8
    assert pocket[refined]


def test_forensic_probe_full_maps_match_june_and_c_only_matches_july() -> None:
    probe = Path("/tmp/ultidock_holo4k_divergence/probe_center_match_summary.csv")
    if not probe.exists():
        pytest.skip("Holo4K forensic probe evidence is not present")

    rows = list(csv.DictReader(probe.open(newline="", encoding="utf-8")))
    assert rows
    for row in rows:
        assert int(row["probe_all_vs_june_set_same"]) == 6
        assert int(row["probe_c_vs_july_set_same"]) == 6


def test_forensic_probe_c_only_changes_dsolv_map() -> None:
    probe = Path("/tmp/ultidock_holo4k_divergence/probe_map_numeric_comparison.csv")
    if not probe.exists():
        pytest.skip("Holo4K forensic probe evidence is not present")

    rows = list(csv.DictReader(probe.open(newline="", encoding="utf-8")))
    assert rows
    by_map = {row["map"]: [] for row in rows}
    for row in rows:
        by_map[row["map"]].append(row)
    assert all(float(row["max_abs_diff"]) == 0.0 for row in by_map["C"])
    assert all(float(row["max_abs_diff"]) == 0.0 for row in by_map["e"])
    assert all(float(row["max_abs_diff"]) > 0.0 for row in by_map["d"])

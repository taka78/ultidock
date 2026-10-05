"""Associate reported scores with the coordinates actually saved by the engine."""

from pathlib import Path


def _atom_signature(lines):
    return tuple((line[6:16], line[30:54]) for line in lines
                 if line.startswith(("ATOM  ", "HETATM")))


def adgpu_best_pose(xml_path):
    """Find the run represented by -best.pdbqt, using the DLG coordinates.

    AutoDock-GPU selects this file by total energy, whereas XML ranks runs by
    binding energy. The first XML run therefore need not be the saved best pose.
    """
    xml_path = Path(xml_path).resolve()
    best_path = xml_path.with_name(xml_path.stem + "-best.pdbqt")
    signature = _atom_signature(best_path.read_text().splitlines())
    if not signature:
        raise ValueError(f"No pose atoms in {best_path}")
    poses = {}
    model = None
    for line in xml_path.with_suffix(".dlg").read_text().splitlines():
        if not line.startswith("DOCKED: "):
            continue
        line = line[8:]
        if line.startswith("MODEL"):
            model = int(line.split()[1])
            poses[model] = []
        elif model is not None:
            poses[model].append(line)
    matches = [run for run, lines in poses.items() if _atom_signature(lines) == signature]
    if len(matches) != 1:
        raise ValueError(f"Cannot uniquely match {best_path.name} to a DLG run: {matches}")
    return best_path, matches[0]


def write_vina_best_pose(output_path, model):
    """Extract the scored Vina model without changing its coordinates or ID."""
    output_path = Path(output_path).resolve()
    selected = []
    inside = False
    for line in output_path.read_text().splitlines():
        if line.startswith("MODEL"):
            inside = int(line.split()[1]) == model
        if inside:
            selected.append(line)
        if line.startswith("ENDMDL") and inside:
            break
    if not _atom_signature(selected):
        raise ValueError(f"No pose atoms for model {model} in {output_path}")
    best_path = output_path.with_name(output_path.stem + "-best.pdbqt")
    best_path.write_text("\n".join(selected) + "\n")
    return best_path

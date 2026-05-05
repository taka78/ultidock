"""Small README pointers shared by the repository CLI commands."""

from __future__ import annotations

from pathlib import Path

README_TOPICS = {
    "requirements": ("Requirements", "requirements"),
    "quick-start": ("Quick Start", "quick-start-end-to-end-run"),
    "commands": ("Command Reference", "command-reference"),
    "receptor-inputs": (
        "Deterministic Receptor Input Handling",
        "deterministic-receptor-input-handling",
    ),
    "benchmarks": ("Benchmarking", "benchmarking"),
    "examples": ("Working with the Example Pipelines", "working-with-the-example-pipelines"),
    "troubleshooting": ("Troubleshooting", "troubleshooting"),
    "site-finder": (
        "Grid Boxing & Cavity Finder Algorithm",
        "spotlight-grid-boxing--cavity-finder-algorithm",
    ),
}


def _repo_root() -> Path | None:
    for candidate in [Path(__file__).resolve().parent.parent, Path.cwd(), *Path.cwd().parents]:
        if (candidate / "README.md").is_file() and (candidate / "docking").is_dir():
            return candidate
    return None


def readme_hint(topic: str) -> str:
    """Return a concise terminal hint pointing at a relevant README section."""
    title, anchor = README_TOPICS.get(topic, README_TOPICS["commands"])
    root = _repo_root()
    readme = root / "README.md" if root else Path("README.md")
    return f"See README: {title} ({readme}#{anchor})"

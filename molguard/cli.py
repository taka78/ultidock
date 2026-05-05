"""Compatibility shim for the MolGuard command-line entry point.

The real CLI definitions live in the repository-level ``cli`` package so the
project can expose separate ``molguard`` and ``ultidock`` commands.
"""

from cli.molguard import cli

__all__ = ["cli"]

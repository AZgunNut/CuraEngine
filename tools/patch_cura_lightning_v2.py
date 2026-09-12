#!/usr/bin/env python3
"""Patch a Cura source/package tree with CuraLightning numeric controls.

Adds two real-valued Lightning controls:
- lightning_smoothing: 0.0..100.0
- lightning_offset_widths: 0.0..100.0 extrusion widths

Usage:
    python tools/patch_cura_lightning_v2.py PATH_TO_CURA_ROOT
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

SMOOTHING_KEY = "lightning_smoothing"
OFFSET_KEY = "lightning_offset_widths"
OLD_SELECTOR_KEY = "lightning_experimental_variant"


def find_resource(root: Path, relative: str) -> Path:
    candidates = [
        root / relative,
        root / "share" / "cura" / "resources" / Path(relative).relative_to("resources"),
        root / "resources" / Path(relative).relative_to("resources"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"Could not find {relative} below {root}")


def patch_definition(path: Path) -> None:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    children = data["settings"]["infill"].setdefault("children", {})
    children.pop(OLD_SELECTOR_KEY, None)

    children[SMOOTHING_KEY] = {
        "label": "Lightning Smoothing",
        "description": "Rounds Lightning paths progressively. 0 keeps the raw Lightning shape; 100 applies intentionally extreme broad smoothing. Smoothing backs off locally to avoid crossing another Lightning path.",
        "type": "float",
        "default_value": 0.0,
        "minimum_value": 0.0,
        "maximum_value": 100.0,
        "enabled": "infill_pattern == 'lightning'",
        "settable_per_mesh": True,
    }

    children[OFFSET_KEY] = {
        "label": "Lightning Offset Distance",
        "description": "Maximum companion-path offset measured literally in extrusion line widths. Decimal values are allowed: 5.6 means 5.6 line widths; 75 means 75 line widths. The offset backs off locally to avoid crossing another Lightning path.",
        "type": "float",
        "unit": "x line width",
        "default_value": 0.0,
        "minimum_value": 0.0,
        "maximum_value": 100.0,
        "enabled": "infill_pattern == 'lightning'",
        "settable_per_mesh": True,
    }

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, indent=4, ensure_ascii=False)
        handle.write("\n")


def patch_visibility(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    lines = [line for line in text.splitlines() if line != OLD_SELECTOR_KEY]

    try:
        infill_index = lines.index("[infill]")
    except ValueError as exc:
        raise RuntimeError(f"No [infill] section found in {path}") from exc

    insert_at = infill_index + 1
    while insert_at < len(lines) and not lines[insert_at].startswith("["):
        insert_at += 1

    for key in (SMOOTHING_KEY, OFFSET_KEY):
        if key not in lines:
            lines.insert(insert_at, key)
            insert_at += 1

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: patch_cura_lightning_v2.py PATH_TO_CURA_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1]).resolve()
    definition = find_resource(root, "resources/definitions/fdmprinter.def.json")
    visibility = find_resource(root, "resources/setting_visibility/expert.cfg")

    patch_definition(definition)
    patch_visibility(visibility)

    with definition.open("r", encoding="utf-8") as handle:
        patched = json.load(handle)
    children = patched["settings"]["infill"]["children"]
    for key in (SMOOTHING_KEY, OFFSET_KEY):
        setting = children.get(key)
        if not setting:
            raise RuntimeError(f"Missing CuraLightning control: {key}")
        if setting.get("minimum_value") != 0.0 or setting.get("maximum_value") != 100.0:
            raise RuntimeError(f"Bad range for CuraLightning control: {key}")

    print(f"Patched {definition}")
    print(f"Patched {visibility}")
    print("Verified Lightning controls: smoothing 0..100; offset 0..100 line widths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

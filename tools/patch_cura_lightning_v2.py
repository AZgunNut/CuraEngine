#!/usr/bin/env python3
"""Patch a Cura source/package tree with the CuraLightning V2 selector.

The custom CuraEngine understands the integer setting
``lightning_experimental_variant`` with values 0..22. This script adds a
visible Cura enum so all variants can be selected in one running Cura session.

Usage:
    python tools/patch_cura_lightning_v2.py PATH_TO_CURA_ROOT
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

SETTING_KEY = "lightning_experimental_variant"


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

    settings = data["settings"]
    infill = settings["infill"]
    children = infill.setdefault("children", {})

    options = {"0": "Stock Lightning"}
    options.update({str(i): f"Experimental {i}" for i in range(1, 23)})

    children[SETTING_KEY] = {
        "label": "Lightning Experimental Variant",
        "description": "Selects the CuraLightning V2 Lightning algorithm. Stock uses normal Cura Lightning; Experimental 1 through 22 select the test variants.",
        "type": "enum",
        "options": options,
        "default_value": "0",
        "enabled": "infill_pattern == 'lightning'",
        "settable_per_mesh": True,
    }

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, indent=4, ensure_ascii=False)
        handle.write("\n")


def patch_visibility(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    if SETTING_KEY in text:
        return

    lines = text.splitlines()
    try:
        infill_index = lines.index("[infill]")
    except ValueError as exc:
        raise RuntimeError(f"No [infill] section found in {path}") from exc

    insert_at = infill_index + 1
    while insert_at < len(lines) and not lines[insert_at].startswith("["):
        insert_at += 1
    lines.insert(insert_at, SETTING_KEY)
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

    # Re-read to make failures explicit in CI.
    with definition.open("r", encoding="utf-8") as handle:
        patched = json.load(handle)
    setting = patched["settings"]["infill"]["children"].get(SETTING_KEY)
    if not setting or len(setting.get("options", {})) != 23:
        raise RuntimeError("CuraLightning V2 selector verification failed")

    print(f"Patched {definition}")
    print(f"Patched {visibility}")
    print("Verified selector: Stock + Experimental 1..22")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

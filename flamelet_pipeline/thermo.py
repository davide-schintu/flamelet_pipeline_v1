from __future__ import annotations

from pathlib import Path
import sys


def load_rftermo_mixture():
    package_dir = Path(__file__).resolve().parent
    local_paths = [
        package_dir / "rfthermo" / "build",
        package_dir / "rfthermo",
    ]
    for path in local_paths:
        if path.exists() and str(path) not in sys.path:
            sys.path.insert(0, str(path))

    try:
        from RFThermo import Mixture
    except ImportError as exc:
        raise ImportError(
            "RFThermo non e' importabile. Build locale suggerito:\n"
            "  cd flamelet_pipeline/rfthermo\n"
            "  python setup.py build_ext --inplace\n"
            "oppure copia il modulo compilato in flamelet_pipeline/rfthermo/build/."
        ) from exc
    return Mixture

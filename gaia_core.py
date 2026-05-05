"""
GAIA_CORE: Módulo fundacional del pipeline Gaia
================================================

Patrones de ingeniería de software científico extraídos del ecosistema
vacuum_energy.py y train_gp.py, adaptados al dominio astronómico.

Provee:
- Custom exception hierarchy (GaiaPipelineError + subtipos)
- Validadores físicos (parallax, distancia, magnitud, RUWE, velocidad)
- Decorador plot_wrapper y context manager plot_context
- Flags de dependencias opcionales (HAS_*)
- Framework de selftest con asserts físicos
- Ejecución paralela con bypass de overhead para N pequeño

Filosofía: fail loud, fail early. Un bug silencioso es infinitamente peor
que un crash explícito.

Autor: Juan Galaz + Claude (Abr 2026)
"""
from __future__ import annotations

import logging
import os
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from functools import wraps
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import numpy as np


# ═════════════════════════════════════════════════════════════════════════════
# 1. DEPENDENCIAS OPCIONALES — GRACEFUL DEGRADATION
# ═════════════════════════════════════════════════════════════════════════════

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

try:
    from joblib import Parallel, delayed
    HAS_JOBLIB = True
except ImportError:
    HAS_JOBLIB = False

try:
    import sympy as sp
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False

try:
    from rich.console import Console
    from rich.table import Table
    HAS_RICH = True
except ImportError:
    HAS_RICH = False


# ═════════════════════════════════════════════════════════════════════════════
# 2. CUSTOM EXCEPTION HIERARCHY
# ═════════════════════════════════════════════════════════════════════════════

class GaiaPipelineError(Exception):
    """Raíz de todas las excepciones del pipeline Gaia."""
    pass


class DataQualityError(GaiaPipelineError):
    """Datos observacionales corruptos, faltantes o fuera de rango físico."""
    pass


class ContractViolationError(GaiaPipelineError):
    """Un script upstream no provee las columnas/tipos que downstream espera.

    Este error existe para matar al bug silencioso que sufrimos en v1 (Script
    06 caía a 'Desconocido' porque las columnas del 03 estaban renombradas).
    """
    pass


class KinematicsError(GaiaPipelineError):
    """Fallo en transformaciones astrométricas o cálculos cinemáticos."""
    pass


class IntegratorError(GaiaPipelineError):
    """Fallo en el integrador N-cuerpos (deriva energética, NaN, divergencia)."""
    pass


class MissingDependencyError(GaiaPipelineError):
    """Dependencia opcional requerida para la feature solicitada."""
    pass


# ═════════════════════════════════════════════════════════════════════════════
# 3. VALIDADORES FÍSICOS
# ═════════════════════════════════════════════════════════════════════════════
#
# Estos son guards. Se llaman al inicio de cada función pública y FALLAN
# ruidosamente si los inputs violan rangos físicos. Nunca silenciosos.

def validate_parallax(parallax_mas: float | np.ndarray, name: str = "parallax") -> None:
    """Valida paralaje en milli-arcsec. Gaia DR3: rango real ~[-5, 500] mas."""
    arr = np.atleast_1d(parallax_mas)
    if np.any(arr < -10):
        raise DataQualityError(
            f"{name} < -10 mas es físicamente imposible (error astrométrico > 10 mas)"
        )
    if np.any(arr > 1000):
        raise DataQualityError(
            f"{name} > 1000 mas: objeto estaría a menos de 1 pc (solo el Sol)"
        )


def validate_distance_pc(d_pc: float | np.ndarray, name: str = "distance") -> None:
    """Valida distancia en parsecs."""
    arr = np.atleast_1d(d_pc)
    if np.any(arr <= 0):
        raise DataQualityError(f"{name} debe ser positiva, got min={np.min(arr):.3e} pc")
    if np.any(arr > 1e6):
        raise DataQualityError(
            f"{name} > 1 Mpc: fuera de galaxia. Verifica unidades (¿parsecs vs kpc?)"
        )


def validate_magnitude(mag: float | np.ndarray, name: str = "magnitude",
                       abs_mag: bool = False) -> None:
    """Valida magnitud aparente o absoluta en rango físicamente plausible."""
    arr = np.atleast_1d(mag)
    if abs_mag:
        # Magnitudes absolutas: M_G ∈ [-30, +25] cubre desde quasars a enanas M
        if np.any(arr < -30) or np.any(arr > 25):
            raise DataQualityError(
                f"{name} absoluta fuera de [-30, +25]: min={np.min(arr):.2f}, max={np.max(arr):.2f}"
            )
    else:
        # Magnitudes aparentes Gaia: G ∈ [3, 22] típicamente
        if np.any(arr < 0) or np.any(arr > 25):
            raise DataQualityError(
                f"{name} aparente fuera de [0, 25]: min={np.min(arr):.2f}, max={np.max(arr):.2f}"
            )


def validate_ruwe(ruwe: float | np.ndarray) -> None:
    """Valida RUWE (Renormalized Unit Weight Error). Valor físico: ~[0.5, 10]."""
    arr = np.atleast_1d(ruwe)
    if np.any(arr < 0):
        raise DataQualityError(f"RUWE no puede ser negativo, got min={np.min(arr):.2f}")
    if np.any(arr > 100):
        logging.warning(
            f"RUWE > 100 detectado (max={np.max(arr):.1f}): posible error de cálculo"
        )


def validate_velocity_kms(v: float | np.ndarray, name: str = "velocity") -> None:
    """Valida velocidad en km/s. Límite físico: escape velocity galáctica ~550 km/s."""
    arr = np.atleast_1d(v)
    finite = arr[np.isfinite(arr)]
    if len(finite) == 0:
        return  # todo NaN, no validamos
    if np.any(np.abs(finite) > 3000):
        raise DataQualityError(
            f"{name} > 3000 km/s detectada: posible error de unidades "
            f"(max|v|={np.max(np.abs(finite)):.1f} km/s)"
        )


def validate_columns_contract(df, required_cols: list[str], producer: str, consumer: str) -> None:
    """Verifica contrato de columnas entre scripts. Falla RUIDOSAMENTE si falta alguna.

    Este validador existe específicamente para matar el bug silencioso de v1
    donde el Script 06 esperaba `luminosidad_absoluta_g` pero el Script 03
    producía `luminosidad_absoluta_g_corregida`, y la clasificación caía a
    NaN sin avisar.
    """
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ContractViolationError(
            f"{consumer} requiere columnas que {producer} no produjo: {missing}. "
            f"Columnas disponibles: {list(df.columns)[:10]}..."
        )


# ═════════════════════════════════════════════════════════════════════════════
# 4. PLOTTING INFRASTRUCTURE
# ═════════════════════════════════════════════════════════════════════════════

def _require_matplotlib() -> None:
    if not HAS_MPL:
        raise MissingDependencyError(
            "matplotlib requerido para plotting. Install: pip install matplotlib"
        )


@contextmanager
def plot_context(style: str = "scientific"):
    """Context manager que aplica estilo y restaura rcParams al salir.

    Evita contaminar el estado global de matplotlib — útil cuando varios
    scripts del pipeline corren en la misma sesión Python (Jupyter, tests).
    """
    _require_matplotlib()
    original_params = plt.rcParams.copy()

    if style == "scientific":
        plt.rcParams.update({
            "font.family": "DejaVu Serif",
            "mathtext.fontset": "stix",
            "font.size": 11,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "lines.linewidth": 1.4,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "grid.linestyle": ":",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
        })
    try:
        yield
    finally:
        plt.rcParams.update(original_params)


def plot_wrapper(filename_stem: str, formats: tuple[str, ...] = ("png", "pdf")):
    """Decorador que estandariza save + cleanup de plots.

    Uso:
        @plot_wrapper("cmd_diagram")
        def plot_cmd(data, dpi=300, out_dir="figures"):
            fig, ax = plt.subplots()
            ...
            return fig   # IMPORTANTE: debe retornar la Figure
    """
    def decorator(plot_func: Callable) -> Callable:
        @wraps(plot_func)
        def wrapper(*args, **kwargs):
            _require_matplotlib()
            dpi = kwargs.get('dpi', 300)
            out_dir = kwargs.get('out_dir', 'figures')

            out_path = Path(out_dir)
            out_path.mkdir(parents=True, exist_ok=True)

            fig = plot_func(*args, **kwargs)
            if fig is None:
                logging.warning(f"plot_func {plot_func.__name__} no retornó Figure")
                return None

            saved = []
            for fmt in formats:
                fname = out_path / f"{filename_stem}.{fmt}"
                fig.savefig(fname, dpi=dpi)
                saved.append(fname)

            logging.info(f"Figura guardada: {saved}")
            plt.close(fig)
            return saved

        return wrapper
    return decorator


# ═════════════════════════════════════════════════════════════════════════════
# 5. PARALLEL EXECUTION
# ═════════════════════════════════════════════════════════════════════════════

def parallel_map(func: Callable, xs: Iterable, threshold: int = 50,
                 prefer: str = "threads") -> np.ndarray:
    """Ejecuta `func` sobre `xs` en paralelo, con bypass para datasets chicos.

    Para N < threshold se ejecuta secuencial porque el overhead de
    serialization + spawning de joblib supera el beneficio del paralelismo.
    Patrón extraído de vacuum_energy.py._parallel_map.
    """
    data_list = list(xs)

    if len(data_list) < threshold or not HAS_JOBLIB:
        return np.array([func(x) for x in data_list])

    max_workers = int(os.getenv('MAX_WORKERS', max(1, (os.cpu_count() or 2) // 2)))
    return np.array(
        Parallel(n_jobs=max_workers, prefer=prefer)(
            delayed(func)(x) for x in data_list
        )
    )


# ═════════════════════════════════════════════════════════════════════════════
# 6. SELFTEST FRAMEWORK
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class TestResult:
    """Resultado de un test individual."""
    name: str
    passed: bool
    value: Any = None
    expected: Any = None
    tolerance: Optional[float] = None
    message: str = ""

    def __str__(self) -> str:
        status = "✔ PASS" if self.passed else "✘ FAIL"
        s = f"[{status}] {self.name}"
        if self.value is not None:
            s += f" (val={self.value}"
            if self.expected is not None:
                s += f", expected={self.expected}"
            s += ")"
        if self.message:
            s += f" — {self.message}"
        return s


class SelfTestSuite:
    """Colección de tests con asserts físicos para validar el pipeline.

    Patrón inspirado en vacuum_energy._selftest: cada test afirma una
    propiedad física que DEBE cumplirse. Si falla, algo se rompió en el
    pipeline.
    """

    def __init__(self, name: str):
        self.name = name
        self.results: list[TestResult] = []

    def assert_close(self, name: str, value: float, expected: float,
                     tolerance: float, relative: bool = True) -> TestResult:
        """Afirma que value ≈ expected dentro de tolerancia."""
        if relative and expected != 0:
            error = abs((value - expected) / expected)
        else:
            error = abs(value - expected)
        passed = error <= tolerance
        result = TestResult(
            name=name, passed=passed, value=value, expected=expected,
            tolerance=tolerance,
            message=f"error_{'rel' if relative else 'abs'}={error:.3e}"
        )
        self.results.append(result)
        return result

    def assert_bounds(self, name: str, value: float,
                      lower: float = -np.inf, upper: float = np.inf) -> TestResult:
        """Afirma que value está en [lower, upper]."""
        passed = lower <= value <= upper
        result = TestResult(
            name=name, passed=passed, value=value,
            expected=f"[{lower}, {upper}]"
        )
        self.results.append(result)
        return result

    def assert_true(self, name: str, condition: bool, message: str = "") -> TestResult:
        """Afirma que `condition` es True."""
        result = TestResult(name=name, passed=bool(condition), message=message)
        self.results.append(result)
        return result

    def summary(self) -> tuple[int, int]:
        """Retorna (n_passed, n_total)."""
        n_passed = sum(1 for r in self.results if r.passed)
        return n_passed, len(self.results)

    def all_passed(self) -> bool:
        return all(r.passed for r in self.results)

    def print_report(self, use_rich: bool = True) -> None:
        """Imprime el reporte. Usa Rich si disponible, fallback a print."""
        n_passed, n_total = self.summary()
        header = f"═══ {self.name} [{n_passed}/{n_total} passed] ═══"

        if HAS_RICH and use_rich:
            from rich.console import Console
            console = Console()
            color = "green" if self.all_passed() else "red"
            console.print(f"[bold {color}]{header}[/bold {color}]")
            for r in self.results:
                c = "green" if r.passed else "red"
                console.print(f"  [{c}]{r}[/{c}]")
        else:
            print(header)
            for r in self.results:
                print(f"  {r}")


# ═════════════════════════════════════════════════════════════════════════════
# 7. DETECCIÓN AUTOMÁTICA DE COLUMNAS (utilidad común)
# ═════════════════════════════════════════════════════════════════════════════

def autodetect_column(df, candidates: list[str], semantic_name: str) -> str:
    """Retorna el primer nombre de columna de `candidates` que existe en df.

    Útil para manejar retrocompatibilidad entre versiones del pipeline
    (p.ej. `V_radial_kms` vs `radial_velocity` vs `dr2_radial_velocity`).
    """
    for c in candidates:
        if c in df.columns:
            return c
    raise ContractViolationError(
        f"No se encontró ninguna columna para '{semantic_name}'. "
        f"Candidatos buscados: {candidates}. "
        f"Columnas disponibles: {list(df.columns)[:15]}..."
    )


# ═════════════════════════════════════════════════════════════════════════════
# 8. CONFIGURACIÓN DE LOGGING ESTÁNDAR
# ═════════════════════════════════════════════════════════════════════════════

def setup_logging(verbosity: int = 0, log_file: Optional[str] = None) -> None:
    """Logging estándar del pipeline con niveles escalonados.

    verbosity: 0=WARNING, 1=INFO, 2+=DEBUG
    """
    level = max(logging.WARNING - 10 * verbosity, logging.DEBUG)

    logger = logging.getLogger()
    logger.setLevel(level)
    if logger.hasHandlers():
        logger.handlers.clear()

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
    logger.addHandler(ch)

    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setFormatter(logging.Formatter(
            "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
        ))
        logger.addHandler(fh)


# ═════════════════════════════════════════════════════════════════════════════
# EXPORTS
# ═════════════════════════════════════════════════════════════════════════════

__all__ = [
    # Exceptions
    "GaiaPipelineError", "DataQualityError", "ContractViolationError",
    "KinematicsError", "IntegratorError", "MissingDependencyError",
    # Validators
    "validate_parallax", "validate_distance_pc", "validate_magnitude",
    "validate_ruwe", "validate_velocity_kms", "validate_columns_contract",
    # Plotting
    "plot_context", "plot_wrapper",
    # Parallel
    "parallel_map",
    # Selftest
    "TestResult", "SelfTestSuite",
    # Utils
    "autodetect_column", "setup_logging",
    # Flags
    "HAS_MPL", "HAS_H5PY", "HAS_JOBLIB", "HAS_SYMPY", "HAS_RICH",
]

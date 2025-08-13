#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refactored COMETS simulations with detailed English comments.

Requirements:
- cometspy
- cobra (COBRApy)
- numpy, pandas, matplotlib

What it does:
- Runs COMETS simulations of the same SBML model under several media/carbon
  conditions (ATCC vs MBM; glucose/lactate/pyruvate; and MBM no-carbon).
- Saves each scenario as a pickle (.pkl) with the full COMETS experiment.
- Plots unified-time growth curves (gDW and estimated OD600) and exports CSVs.
"""

import os
import pickle
from pathlib import Path
from typing import Dict, Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cobra
import cometspy as c


# =========================
# Paths and global settings
# =========================

# SBML model path (edit if needed)
SBML_PATH = "refs/kbase/Hot5F3RastDraftModelATCCMediaGapFillModel.xml"

# Output directory
OUT_DIR = Path("./outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# COMETS install path
DEFAULT_COMETS_HOME =  "/Users/huangyu/Desktop/internship/comets_macos/comets_2.12.3/"

# Simulation settings
TIME_STEP_H = 0.01             # hours per step
SIM_HOURS   = 20.0             # total simulated hours
MAX_CYCLES  = int(SIM_HOURS / TIME_STEP_H)
DEFAULT_VMAX = 18.5            # default Vmax for uptake
OXYGEN_VMAX  = 10.0            # Vmax for oxygen exchange
DEFAULT_KM   = 1.5e-5          # default Km
INITIAL_BIOMASS_GDW = 5e-6     # initial total biomass (gDW)
SPACE_WIDTH = 1                # 1x1 grid
MAX_SPACE_BIOMASS = 10
MIN_SPACE_BIOMASS = 1e-11

# Unified plotting time step
UNIFIED_DT = 0.05  # hours

# OD600 estimation parameters (rough conversion)
VOLUME_L = 0.001                # COMETS test tube ~ 1 cm^3 = 1 mL
GDW_PER_OD_PER_L = 0.30         # 1 OD600 ≈ 0.30 gDW/L (replace with your value)


# =========================
# Media and carbon sources
# =========================

# Oxygen (always added as a limited amount)
OXYGEN_ID = 'cpd00007_e0'
OXYGEN_AMOUNT = 2.0  # mmol (in the 1 cm^3 volume)

# Carbon sources (use one at a time; None for "no carbon")
CARBON_SOURCES = {
    "glucose":  ('cpd00027_e0', 0.012),
    "lactate":  ('cpd00221_e0', 0.012),
    "pyruvate": ('cpd00020_e0', 0.012),
}

# Exchange reactions to open (avoid COBRA bounds limiting COMETS uptake)
EX_TO_OPEN = [
    'EX_cpd00027_e0',  # glucose
    'EX_cpd00221_e0',  # lactate
    'EX_cpd00020_e0',  # pyruvate
    'EX_cpd00029_e0',  # acetate (often secreted/consumed)
]

# ATCC-like base media (excluding oxygen/carbon source)
ATCC_BASE_MEDIA = {
    'cpd00058_e0': 1000,  # Cu2+
    'cpd00971_e0': 1000,  # Na+
    'cpd00063_e0': 1000,  # Ca2+
    'cpd00048_e0': 1000,  # Sulfate
    'cpd10516_e0': 1000,  # Fe3+
    'cpd00254_e0': 1000,  # Mg2+
    'cpd00009_e0': 1000,  # Phosphate
    'cpd00205_e0': 1000,  # K+
    'cpd00013_e0': 1000,  # NH3/NH4+
    'cpd00099_e0': 1000,  # Cl-
    'cpd00030_e0': 1000,  # Mn2+
    'cpd00001_e0': 1000,  # H2O
    'cpd00034_e0': 1000,  # Zn2+
    'cpd00209_e0': 1000,  # NO3-
    'cpd00011_e0': 1000,  # CO2
    'cpd00067_e0': 1000,  # H+
    'cpd17005_e0': 1000,  # CO3(2-)
    'cpd00240_e0': 1000,  # EDTA
    'cpd01843_e0': 1000,  # Citric acid
    'cpd00149_e0': 1000,  # Co2+
    'cpd00423_e0': 1000,  # Vitamin B12r
}

# MBM-like (defined/minimal) base media
MBM_BASE_MEDIA = {
    'cpd00971_e0': 1000,  # Na+
    'cpd00099_e0': 1000,  # Cl-
    'cpd00048_e0': 1000,  # Sulfate
    'cpd00205_e0': 1000,  # K+
    'cpd00966_e0': 1000,  # Br-
    'cpd09225_e0': 1000,  # Boric acid (H3BO3)
    'cpd00552_e0': 1000,  # F-
    'cpd00242_e0': 1000,  # HCO3-
    'cpd00254_e0': 1000,  # Mg2+
    'cpd00063_e0': 1000,  # Ca2+
    'cpd09695_e0': 1000,  # Sr2+
    'cpd10516_e0': 1000,  # Fe3+ (FeEDTA)
    'cpd00240_e0': 1000,  # EDTA
    'cpd28238_e0': 1000,  # tris-HCl
    'cpd00013_e0': 1000,  # NH3/NH4+
    'cpd00009_e0': 1000,  # Phosphate
    'cpd01301_e0': 1000,  # Biotin
    'cpd00393_e0': 1000,  # Folate
    'cpd00263_e0': 1000,  # Pyridoxine
    'cpd00220_e0': 1000,  # Riboflavin
    'cpd00793_e0': 1000,  # Thiamine
    'cpd00133_e0': 1000,  # Nicotinic acid
    'cpd00644_e0': 1000,  # Pantothenic acid
    'cpd01826_e0': 1000,  # Cyanocobalamin
    'cpd00443_e0': 1000,  # p-Aminobenzoic acid
    # Extra metals needed for growth
    'cpd00058_e0': 1000,  # Cu2+
    'cpd00030_e0': 1000,  # Mn2+
    'cpd00034_e0': 1000,  # Zn2+
    'cpd00149_e0': 1000,  # Co2+
}


# =========================
# Utility functions
# =========================

def ensure_comets_home() -> str:
    """Ensure COMETS_HOME is set and exists; return the path."""
    path = os.environ.get("COMETS_HOME", "").strip()
    if not path:
        path = DEFAULT_COMETS_HOME
        os.environ["COMETS_HOME"] = path
    if not Path(path).exists():
        raise FileNotFoundError(
            f"COMETS_HOME does not exist: {path}\n"
            f"Set COMETS_HOME or edit DEFAULT_COMETS_HOME."
        )
    return path


def load_cobra_model(sbml_path: str) -> cobra.Model:
    """Load the SBML model once and reuse it."""
    print(f"[INFO] Loading SBML model: {sbml_path}")
    return cobra.io.read_sbml_model(sbml_path)


def build_comets_model(cobra_model: cobra.Model) -> c.model:
    """
    Wrap the COBRA model into a COMETS model, set pFBA objective style,
    open selected exchange bounds, and set initial biomass.
    """
    m = c.model(cobra_model)
    m.obj_style = 'MAX_OBJECTIVE_MIN_TOTAL'
    for ex_id in EX_TO_OPEN:
        try:
            m.change_bounds(ex_id, -1000, 1000)
        except Exception as e:
            # Some exchanges may not exist in a given model; that's fine.
            print(f"[WARN] Exchange not found: {ex_id} ({e})")
    # Initial biomass at grid (0, 0)
    m.initial_pop = [0, 0, INITIAL_BIOMASS_GDW]
    return m


def make_layout(base_media: Dict[str, float],
                carbon: Optional[Tuple[str, float]] = None,
                oxygen_amount: float = OXYGEN_AMOUNT) -> c.layout:
    """
    Create a 1x1 layout and set media.
    - base_media: background nutrients (excluding O2 and carbon)
    - carbon: (compound_id, amount) or None for no-carbon
    - oxygen_amount: total amount of oxygen to add
    """
    L = c.layout()
    for met, val in base_media.items():
        L.set_specific_metabolite(met, float(val))
    L.set_specific_metabolite(OXYGEN_ID, float(oxygen_amount))
    if carbon is not None:
        met, amt = carbon
        L.set_specific_metabolite(met, float(amt))
    return L


def make_params(time_step: float = TIME_STEP_H,
                max_cycles: int = MAX_CYCLES,
                default_vmax: float = DEFAULT_VMAX,
                default_km: float = DEFAULT_KM) -> c.params:
    """Create a COMETS params object with commonly used settings."""
    P = c.params()
    P.set_param('timeStep', time_step)
    P.set_param('maxCycles', max_cycles)
    P.set_param('defaultVmax', default_vmax)
    P.set_param('defaultKm', default_km)
    P.set_param('spaceWidth', SPACE_WIDTH)
    P.set_param('maxSpaceBiomass', MAX_SPACE_BIOMASS)
    P.set_param('minSpaceBiomass', MIN_SPACE_BIOMASS)
    P.set_param('writeMediaLog', True)
    P.set_param('writeFluxLog', True)
    P.set_param('FluxLogRate', 1)
    return P


def run_experiment(tag: str,
                   base_media: Dict[str, float],
                   carbon: Optional[Tuple[str, float]],
                   cobra_model: cobra.Model,
                   out_dir: Path = OUT_DIR) -> Path:
    """
    Run one scenario and save the full COMETS experiment as a pickle.
    - tag: name used for output (e.g., 'ATCC_Glucose')
    - base_media: ATCC or MBM base media (dict)
    - carbon: (compound_id, amount) or None
    - cobra_model: preloaded COBRA model
    Returns: path to the saved .pkl file
    """
    print(f"\n[INFO] Running scenario: {tag}")

    ensure_comets_home()

    # 1) Build media/layout
    layout = make_layout(base_media=base_media, carbon=carbon)

    # 2) Build COMETS model
    model = build_comets_model(cobra_model)

    # 3) Set a specific Vmax for the oxygen exchange reaction
    #    The exchange id is typically 'EX_cpd00007_e0'
    model.change_vmax('EX_cpd00007_e0', OXYGEN_VMAX)

    # 4) Add model to the layout
    layout.add_model(model)

    # 5) Simulation parameters
    params = make_params()

    # 6) Create and run the experiment
    experiment = c.comets(layout, params)
    experiment.run()

    # 7) Save the experiment object (.pkl)
    out_pkl = out_dir / f"{tag}.pkl"
    with open(out_pkl, "wb") as f:
        pickle.dump(experiment, f)
    print(f"[INFO] Saved: {out_pkl}")

    return out_pkl


def extract_time_biomass(pkl_path: Path) -> Tuple[np.ndarray, np.ndarray, float]:
    """Load a saved experiment and return (time[h], biomass[gDW], dt[h])."""
    with open(pkl_path, "rb") as f:
        sim = pickle.load(f)
    biomass = sim.total_biomass.iloc[:, 1].to_numpy()  # column 1 = total gDW
    dt = float(sim.parameters.get_param('timeStep'))
    time = np.arange(len(biomass)) * dt
    return time, biomass, dt


def biomass_to_od600(biomass_gdw: np.ndarray,
                     volume_L: float = VOLUME_L,
                     alpha_gdw_per_OD_per_L: float = GDW_PER_OD_PER_L) -> np.ndarray:
    """Convert total gDW to a rough OD600 estimate."""
    return (biomass_gdw / volume_L) / alpha_gdw_per_OD_per_L


def plot_and_export(series: Dict[str, Tuple[np.ndarray, np.ndarray]],
                    unified_dt: float = UNIFIED_DT,
                    out_dir: Path = OUT_DIR) -> None:
    """
    Plot all scenarios on a unified time axis.
    - gDW curves + PNG
    - Estimated OD600 curves + PNG
    - Export both as CSV
    """
    max_time = max(t[-1] for t, _ in series.values())
    unified_time = np.arange(0, max_time + unified_dt, unified_dt)

    # Plot gDW
    plt.figure(figsize=(10, 6))
    csv_gdw = {"Time (h)": unified_time}

    for name, (time, biomass) in series.items():
        aligned = np.interp(unified_time, time, biomass)
        plt.plot(unified_time, aligned, label=name, linewidth=1.6)
        csv_gdw[name] = aligned

    plt.xlabel("Time (hours)")
    plt.ylabel("Total Biomass (gDW)")
    plt.title("Growth curves under different media/carbon sources")
    plt.legend()
    plt.tight_layout()
    fig_path = out_dir / "growth_curves_gDW.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[INFO] Saved: {fig_path}")

    # Plot estimated OD600
    plt.figure(figsize=(10, 6))
    csv_od = {"Time (h)": unified_time}

    for name, (time, biomass) in series.items():
        aligned = np.interp(unified_time, time, biomass)
        od = biomass_to_od600(aligned, VOLUME_L, GDW_PER_OD_PER_L)
        plt.plot(unified_time, od, label=f"{name} (OD600)", linewidth=1.6)
        csv_od[name] = od

    plt.xlabel("Time (hours)")
    plt.ylabel("OD600 (arb.)")
    plt.title("Estimated OD600 (requires volume and conversion factor)")
    plt.legend()
    plt.tight_layout()
    fig_path = out_dir / "growth_curves_OD600.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[INFO] Saved: {fig_path}")

    # CSV exports
    pd.DataFrame(csv_gdw).to_csv(out_dir / "all_growth_curves_gDW.csv", index=False)
    pd.DataFrame(csv_od ).to_csv(out_dir / "all_growth_curves_OD600.csv", index=False)
    print(f"[INFO] Saved CSVs to: {out_dir.resolve()}")


# =========================
# Main pipeline
# =========================

def main():
    # 1) Load the SBML model once
    cobra_model = load_cobra_model(SBML_PATH)

    # 2) Define scenarios: (output_name, base_media, carbon_source_or_None)
    scenarios: List[Tuple[str, Dict[str, float], Optional[Tuple[str, float]]]] = [
        # ATCC media
        ("ATCC_Glucose",  ATCC_BASE_MEDIA, CARBON_SOURCES["glucose"]),
        ("ATCC_Lactate",  ATCC_BASE_MEDIA, CARBON_SOURCES["lactate"]),
        ("ATCC_Pyruvate", ATCC_BASE_MEDIA, CARBON_SOURCES["pyruvate"]),
        # MBM media
        ("MBM_NoCarbon",  MBM_BASE_MEDIA,  None),
        ("MBM_Glucose",   MBM_BASE_MEDIA,  CARBON_SOURCES["glucose"]),
        ("MBM_Lactate",   MBM_BASE_MEDIA,  CARBON_SOURCES["lactate"]),
        ("MBM_Pyruvate",  MBM_BASE_MEDIA,  CARBON_SOURCES["pyruvate"]),
    ]

    # 3) Run all scenarios and store paths
    result_paths: Dict[str, Path] = {}
    for name, base_media, carbon in scenarios:
        pkl_path = run_experiment(tag=name,
                                  base_media=base_media,
                                  carbon=carbon,
                                  cobra_model=cobra_model,
                                  out_dir=OUT_DIR)
        result_paths[name] = pkl_path

    # 4) Load results for plotting/export
    series: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for name, p in result_paths.items():
        time, biomass, dt = extract_time_biomass(p)
        series[name] = (time, biomass)

    # 5) Plot and export
    plot_and_export(series, unified_dt=UNIFIED_DT, out_dir=OUT_DIR)

    print("\n[ALL DONE] Results saved under:", OUT_DIR.resolve())


if __name__ == "__main__":
    main()
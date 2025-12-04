"""
Core simulation code for the DNA replication asymmetry / proteostasis model.

This module provides:
- step_mutations: symmetric vs asymmetric mutation generation
- simulate_load: single-cell load trajectories with optional aging decline
- simulate_population: tissue-scale collapse and survival
- survival_curve: fraction alive over time
- collapse_stats: median, quartiles
- run_sensitivity: parameter sweep over d, B, and mutation rates
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# MODEL PARAMETERS
# ---------------------------------------------------------------------------

@dataclass
class ModelParams:
    """Encapsulates all parameters for replication & proteostasis model."""
    # Genomic / mutation parameters
    G: float = 12.0
    f_E: float = 0.7
    f_L: float = 0.3

    MU_BASE: float = 0.003
    MU_LE_E: float = 0.002
    MU_LA_E: float = 0.006
    MU_LE_L: float = 0.0024
    MU_LA_L: float = 0.0072

    # Proteostasis parameters
    p_mis: float = 0.5
    Y: float = 1000.0
    B: float = 2000.0


# ---------------------------------------------------------------------------
# MUTATION STEP
# ---------------------------------------------------------------------------

def step_mutations(rng: np.random.Generator, asym: bool, size: int, p: ModelParams):
    """Return mutation counts per lineage for symmetric (H0) or asymmetric (H1)."""
    if not asym:
        return rng.poisson(p.G * p.MU_BASE, size=size)

    # Asymmetric case: leading/lagging × early/late
    mut_le_e = rng.poisson(p.G * p.f_E * p.MU_LE_E, size=size)
    mut_la_e = rng.poisson(p.G * p.f_E * p.MU_LA_E, size=size)
    mut_le_l = rng.poisson(p.G * p.f_L * p.MU_LE_L, size=size)
    mut_la_l = rng.poisson(p.G * p.f_L * p.MU_LA_L, size=size)

    return mut_le_e + mut_la_e + mut_le_l + mut_la_l


# ---------------------------------------------------------------------------
# SINGLE-CELL SIMULATION
# ---------------------------------------------------------------------------

def simulate_load(
    asym: bool,
    reps: int = 1000,
    T: int = 200,
    d: float = 0.8,
    gamma: float = 0.0,
    p: Optional[ModelParams] = None,
    Lcrit: Optional[float] = None,
    seed: int = 123
):
    """Simulate proteostasis load in parallel for many lineages."""
    if p is None:
        p = ModelParams()

    rng = np.random.default_rng(seed)
    L = np.zeros(reps)
    collapsed = np.zeros(reps, dtype=bool)
    collapse_time = np.full(reps, np.nan)

    for t in range(T):
        d_t = max(d - gamma * t, 0.1) if gamma > 0 else d

        M = step_mutations(rng, asym, reps, p)
        misfold = rng.binomial(M, p.p_mis) * p.Y
        total_new = misfold + p.B

        L = (L + total_new) * (1 - d_t)

        if Lcrit is not None:
            newly = (~collapsed) & (L >= Lcrit)
            collapse_time[newly] = t + 1
            collapsed[newly] = True

    if Lcrit is None:
        return L
    return L, collapsed, collapse_time


# ---------------------------------------------------------------------------
# TISSUE-SCALE SIMULATION
# ---------------------------------------------------------------------------

def simulate_population(
    d: float,
    gamma: float,
    B: float,
    Lcrit: float,
    T: int = 2000,
    reps: int = 3000,
    seed: int = 999,
    p: Optional[ModelParams] = None,
):
    """Simulate many stem-cell–like lineages with aging decline in clearance."""
    if p is None:
        p = ModelParams(B=B)
    else:
        p = ModelParams(**{**p.__dict__, "B": B})

    rng = np.random.default_rng(seed)

    def run(asym: bool):
        L = np.zeros(reps)
        collapsed = np.zeros(reps, dtype=bool)
        collapse_time = np.full(reps, np.nan)

        for t in range(T):
            d_t = max(d - gamma * t, 0.1) if gamma > 0 else d

            M = step_mutations(rng, asym, reps, p)
            misfold = rng.binomial(M, p.p_mis) * p.Y
            total_new = misfold + p.B

            L = (L + total_new) * (1 - d_t)

            newly = (~collapsed) & (L >= Lcrit)
            collapse_time[newly] = t + 1
            collapsed[newly] = True

        return collapsed, collapse_time

    coll0, t0 = run(False)
    coll1, t1 = run(True)

    return coll0, t0, coll1, t1


# ---------------------------------------------------------------------------
# SURVIVAL CURVE
# ---------------------------------------------------------------------------

def survival_curve(collapse_times: np.ndarray, T: int):
    """Return fraction alive (not collapsed) at each time step."""
    alive = []
    for t in range(1, T + 1):
        alive.append(np.mean((np.isnan(collapse_times)) | (collapse_times > t)))
    return np.array(alive)


# ---------------------------------------------------------------------------
# COLLAPSE STATS
# ---------------------------------------------------------------------------

def collapse_stats(times: np.ndarray, label: str):
    """Return median collapse time and quartiles for a set of lineages."""
    valid = times[~np.isnan(times)]
    if len(valid) == 0:
        return {"label": label, "n": 0, "median": np.nan, "q25": np.nan, "q75": np.nan}

    return {
        "label": label,
        "n": int(len(valid)),
        "median": float(np.median(valid)),
        "q25": float(np.percentile(valid, 25)),
        "q75": float(np.percentile(valid, 75)),
    }


# ---------------------------------------------------------------------------
# SENSITIVITY ANALYSIS
# ---------------------------------------------------------------------------

def run_sensitivity(
    d_vals=(0.6, 0.7, 0.8, 0.9),
    B_vals=(1000.0, 2000.0, 4000.0),
    mu_factors=(0.5, 1.0, 2.0),
    reps: int = 600,
    T: int = 200
):
    """
    Sweep over clearance (d), basal misfolding (B), and mutation-rate scaling.

    Returns DataFrame with Δ% load between H0 and H1.
    """
    rows = []
    base = ModelParams()

    for d_val in d_vals:
        for B_val in B_vals:
            for mf in mu_factors:

                p_scaled = ModelParams(
                    G=base.G,
                    f_E=base.f_E,
                    f_L=base.f_L,
                    MU_BASE=base.MU_BASE * mf,
                    MU_LE_E=base.MU_LE_E * mf,
                    MU_LA_E=base.MU_LA_E * mf,
                    MU_LE_L=base.MU_LE_L * mf,
                    MU_LA_L=base.MU_LA_L * mf,
                    p_mis=base.p_mis,
                    Y=base.Y,
                    B=B_val
                )

                L0 = simulate_load(False, reps=reps, T=T, d=d_val, p=p_scaled, seed=111)
                L1 = simulate_load(True,  reps=reps, T=T, d=d_val, p=p_scaled, seed=222)

                mean0 = L0.mean()
                mean1 = L1.mean()
                delta = (mean1 - mean0) / mean0 * 100

                rows.append({
                    "d": d_val,
                    "B": B_val,
                    "mu_factor": mf,
                    "mean_L_H0": mean0,
                    "mean_L_H1": mean1,
                    "delta_pct": delta
                })

    return pd.DataFrame(rows)

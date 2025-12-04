# DNA Replication Asymmetry and Proteostasis Collapse Dynamics
*A stochastic modeling framework for mutation-driven proteotoxic stress and failure thresholds*

---

## Overview

This repository contains a computational framework for quantifying how **DNA replication asymmetry**—differences in mutation input between leading and lagging strands—affects long-term **proteostasis load**, **steady-state misfolding**, and **collapse probability** in long-lived proliferative tissues.

The model is relevant for biological systems where cells undergo **hundreds to thousands of divisions**, experience **strong proteostasis early in life**, and show **gradual age-related clearance decline**. Examples include:

**hematopoietic stem and progenitor cells, intestinal stem cells, and epidermal basal progenitors**,  
where experimental studies have documented:

- replication-dependent mutation asymmetry  
- high basal protein-turnover burden (large \(B\))  
- slow age-dependent loss of proteostasis capacity (γ ≈ 10⁻³)  
- long-term maintenance under stochastic misfolding stress  

To avoid confusion:  
**This project models replication-induced asymmetry in DNA synthesis, not cell-division or partitioning asymmetry.**

---

## Scientific components implemented

The model follows the standard **load–capacity** framework widely used in theoretical proteostasis and aging research. It incorporates:

- Poisson mutation processes with strand-specific asymmetry  
- Binomial misfolding amplification  
- Basal proteotoxic input \(B\)  
- Clearance decline \(d_t = \max(d - \gamma t, d_{\min})\)  
- Analytical derivation of steady-state load \(μ^*\)  
- Stochastic first-passage collapse at threshold \(L_{\text{crit}}\)  
- Tissue-level survival curves across thousands of simulated lineages  
- Sensitivity analyses (varying \(B\), \(d\), \(\gamma\), mutation-rate scaling)  
- Random-seed robustness tests  

Deterministic expressions are used **only for deriving equilibrium expectations**.  
All collapse/survival outcomes arise exclusively from **full stochastic simulations**.

Sensitivity analyses vary **parameters**, not structural model assumptions.

---

## Repository structure

```
project_root/
│
├── src/
│   └── model_proteostasis.py          
│
├── notebooks/
│   ├── 01_mutation_sanity.ipynb       
│   ├── 02_load_dynamics.ipynb        
│   ├── 03_sensitivity.ipynb           
│   ├── 04_population_near_failure.ipynb
│   ├── 05_population_stress_test.ipynb
│   └── 06_seed_robustness.ipynb
│
├── data/
│   ├── README.md                      
│
├── results/
│   ├── mutation_sanity_stats.csv
│   ├── load_dynamics_summary.csv
│   ├── B_sensitivity.csv
│   ├── Lcrit_sweep.csv
│   ├── stress_test_stats.csv
│   └── seed_robustness.csv
│
└── figures/
    ├── mutation_histogram.png
    ├── load_distribution.png
    ├── B_sensitivity_plot.png
    ├── survival_near_failure.png
    └── survival_stress_test.png
```

Executing the notebooks will recreate all figures and CSV files.

---

## Installation

Install dependencies:

```
pip install -r requirements.txt
```

Minimal `requirements.txt`:

```
numpy
pandas
matplotlib
```



---

## Running the model

1. Launch JupyterLab  
2. Open the `notebooks/` directory  
3. Run notebooks **01 → 06** sequentially  
4. Figures and CSVs will be generated in `figures/` and `results/`

---

## Reproducibility

- All simulations use explicit NumPy random seeds.  
- All mathematical definitions and parameters are centralized in `src/model_proteostasis.py`.  
- No hidden state; no notebook overrides.  
- Results reproduce exactly across systems and runs.

This structure meets reproducibility standards for computational biology and theoretical biophysics publications.

---

## License

Released under the **MIT License**, allowing reuse, adaptation, and extension while preserving author attribution.

---

## Citation

If you use this code or its analytical derivations, please cite:

```
Regi Avdiaj, 2025. "DNA Replication Asymmetry and Proteostasis Collapse Dynamics."
```


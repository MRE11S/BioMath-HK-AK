# Land-based spread of the bubonic plague across Europe (1346–1353): reaction–diffusion model with documented routes

MATLAB code for the Biomathematics 374 report on the land-based spread of the Second Pandemic of plague.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=YOUR-USER/YOUR-REPO&file=bubonic_plague_PDE_model.m)

> Replace `YOUR-USER/YOUR-REPO` in the badge link with your GitHub user and repository name.

## Files

| File | Purpose |
|---|---|
| `bubonic_plague_PDE_model.m` | the 3D model with multiple entry points; produces every result in the article |
| `documented_arrival_year_grid.csv` | documented arrival year (1346–1353) on the model grid, digitised from Map 1 of Benedictow (2021); must stay in the same folder as the script |

## What the code does

It solves a susceptible–infected–removed (SIR) reaction–diffusion system with births and deaths on a 0.25° grid of Europe (12°W–50°E, 30–62°N) built from real topography. **Transmission happens over connected land only**: no sea or port spread and no flea or rat compartment.

Data taken from Benedictow (2021), *The Complete History of the Black Death* (Boydell Press), Map 1, "The advance of the Black Death year by year":

1. **Corridors.** Mobility is raised along 25 land corridors that follow the land arrows of Map 1 (city chains are listed in the script). Sea arrows are not used.
2. **Entry points.** The documented 1346 zone is infected on day 0; the documented 1347 zone and 21 listed ports are infected in the middle of their documented year (`prm.entryMode`, `prm.entryDay`), only if the land front has not already arrived. Infection never crosses the sea.
3. **Barrier.** The border between the Golden Horde and the Russian principalities, marked on the map as "passage prohibited for Christians", is a zero-flux wall (`prm.barrierFactor = 0`).
4. **Comparison.** The year-by-year colour zones were digitised onto the model grid (`documented_arrival_year_grid.csv`) and are used only to compare with the model, never imposed on it. The map was georeferenced on 29 cities (root-mean-square error about 35 km). Blank ("nan") cells are "no data" or outside the map.

The script then runs the 2D baseline (1 January 1346 – 31 December 1353), maps arrival times and safe zones, compares model and documented arrival years, verifies the front speed against the Fisher–KPP formula in 1D, tests the R0 = 1 threshold,.

Sea-borne arrivals documented on the map are taken as entry points (inputs), so the comparison with the map excludes all entry cells and is not an independent validation of the entry years.

## Requirements

* MATLAB R2019b or newer (the script uses `readmatrix`, `exportgraphics`, `movmax`, `smooth3` and `isosurface`).
* No toolboxes. No downloads: the topography is MATLAB's built-in `topo` data set.
* A few minutes of run time and modest memory.

## How to run

1. Clone or download this repository (keep both files together).
2. Open `bubonic_plague_PDE_model.m` in MATLAB (or click the *Open in MATLAB Online* badge).
3. Press **Run**.

All parameters are in one annotated block at the top (section `1. PARAMETERS`); an assumptions register with the reason for each assumption is in the header comment. Switches: `runExtras` (1D verification) and `makeVideo` (MP4 animation, skipped automatically if unavailable).

## Outputs (folder `BubonicPlague_PDE_Outputs/`)

| Content | File |
|---|---|
| Domain, corridors, local R0 | `Fig1_domain_routes_R0.png` |
| Wave snapshots | `Fig2_3D_wave_snapshots.png` |
| Arrival time | `Fig3_arrival_time.png` |
| Space–time isosurface (extra) | `Fig4_spacetime_isosurface.png` |
| Safe zones | `Fig5_safe_zones.png` |
| Documented vs model arrival year, difference map, city curves | `Fig6_documented_vs_model.png` |
| 1D verification | `Fig7_1D_validation_threshold.png` |
| Animation | `plague_PDE_wave_3D.mp4` |
| Tables | `city_table.csv`, `documented_comparison.csv`, `documented_zone_reach.csv`, `safe_zone_regions.csv`, `wave_speed_1D.csv` |
| Summary | `summary.txt`, `results.mat` |

## Main parameters

| Symbol | Value | Meaning |
|---|---|---|
| β | 0.30 day⁻¹ | transmission rate |
| γ | 0.15 day⁻¹ | removal rate |
| μ | 5 × 10⁻⁵ day⁻¹ | birth/death rate |
| D0 | 5 km² day⁻¹ | background diffusion over all land |
| b | 1 | corridor boost (corridor D = 10 km² day⁻¹) |
| φ | 0.45 | fraction of removed individuals who die |

Reference R0 = β/(γ+μ) ≈ 2.0.

## References

* Benedictow, O. J. *The Complete History of the Black Death* (Boydell Press, 2021): Map 1.
* Davydovych, V., Dutka, V. A. & Cherniha, R. Reaction–diffusion equations in mathematical models arising in epidemiology. *Symmetry* **15**, 2025 (2023).

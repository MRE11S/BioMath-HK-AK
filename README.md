# Land-based spread of the bubonic plague across Europe (1346–1353): reaction–diffusion model with documented routes

MATLAB code for the Biomathematics 374 report on the land-based spread of the Second Pandemic of plague.

## Files

| File | Purpose |
|---|---|
| `bubonic_plague_AK_HK.m` | the 3D model with multiple entry points; produces every result in the article |
| `documented_arrival_year_grid.csv` | documented arrival year (1346–1353) on the model grid, digitised from Map 1 of Benedictow (2021); must stay in the same folder as the script |

## Requirements

* MATLAB R2019b or newer (the script uses `readmatrix`, `exportgraphics`, `movmax`, `smooth3` and `isosurface`).
* No toolboxes. No downloads: the topography is MATLAB's built-in `topo` data set.
* A few minutes of run time and modest memory.

## How to run

1. Clone or download this repository (keep both files together).
2. Open `bubonic_plague_AK_HK.m` in MATLAB (or click the *Open in MATLAB Online* badge).
3. Press **Run**.

All parameters are in one annotated block at the top (section `1. PARAMETERS`); an assumptions register with the reason for each assumption is in the header comment. Switches: `runExtras` (1D verification) and `makeVideo` (MP4 animation, skipped automatically if unavailable).


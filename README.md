# BlueCrabs

Bioeconomic model and data for *"Timely Policy Responses Can Save Fisheries from (Some) Temperature Shocks"*.

## Repository contents

```text
functions/                 Model functions
inputs/                    Model inputs and parameter files
results/                   Simulation outputs used to reproduce manuscript figures

Variable_setup.R           Initializes model parameters
run_optimization_on_MPI.R  Executes large-scale parallel simulations

fig1_c_d.R                 Reproduces Fig. 1c,d
fig2.R                     Reproduces Fig. 2
fig3.R                     Reproduces Fig. 3
fig4.R                     Reproduces Fig. 4
```

## Software requirements

* R version 4.6.0 or later
* RStudio, optional
* Required R packages are listed in `Variable_setup.R`

Large-scale simulations were executed on a high-performance computing cluster using MPI. Figure-generation scripts can be run on a standard desktop computer using the processed simulation outputs included in this repository.

## Installation

Clone or download this repository. Typical installation time is less than 5 minutes.

## Reproducing the manuscript

The repository contains the processed simulation outputs required to reproduce all manuscript figures.

To regenerate the figures, run:

```text
fig1_c_d.R
fig2.R
fig3.R
fig4.R
```

The script `run_optimization_on_MPI.R` executes the full simulation framework using MPI-based parallel processing on a high-performance computing cluster. It is provided to document the complete simulation workflow and requires an MPI environment. It is not intended to be executed on a standard desktop computer.

## Data availability

The repository contains the processed data required to reproduce the figures presented in the manuscript.

Fishery-independent trawl survey and commercial landings data used for model parameterization and validation were provided by the Delaware Department of Natural Resources and Environmental Control (DNREC) and are not publicly available.

# survSched
Code for the [surveillance scheduling paper](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13341).

To install the package:
```
devtools::install_github("DRuiCHEN/survSched")
```

To reproduce the simulation results and the figures:
1. Download the `paper/` folder.
2. Install the package and run `simu_run.R` (or `simu_run_supp.R` for the experiments in the supplementary materials). 
The results from this step are already stored in `paper/rds/`.
3. `simu_plot.R` (or `simu_plot_supp.R`) loads the results in `paper/rds/` and produces the figures in the paper.


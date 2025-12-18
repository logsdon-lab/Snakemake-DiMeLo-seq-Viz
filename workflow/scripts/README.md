# Scripts

## `plot_enrichment.py`
Original script (`analysis_workflow_unfilter30.py`) contributed by Shenghan Gao. Normalizes the treatment signal by the control and generates `cenplot` images.

Briefly, it:
* Generates 5 kbp non-overlapping windows of the control and treatment pileup.
* Subtracts the control pileup from the treatment pileup.
* Generate a cenplot config toml file.
* Plots the CENP-A signal.

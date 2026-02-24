# Experimental modules (optional)

These scripts are **experimental** and may require additional dependencies and dataset-specific endpoint mapping.

They are provided as a starting point only and are **not** required to run the core reproducible pipeline.

## Contents
- `survival_101/`: template for running a larger survival-model comparison (multiple algorithms / combinations)

## Important
- You must define a survival endpoint (`time`, `event`) consistently across cohorts.
- Some models require extra R packages or external resources.

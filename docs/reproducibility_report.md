# Computational Reproducibility Report

## Project Information

- Student: Yashwant Nama
- Institution: Atria University, Bengaluru
- Programme: MS by Research in Life Sciences
- Course: LS-604 Computational Biology
- Project: P0
- Date: 27 November 2026

## Objective

To assess computational reproducibility by independently reconstructing a documented Conda environment on a second computer and comparing analytical outputs.

## Experimental Setup

| Parameter | Original System | Reconstructed System |
|---|---|---|
| Computer | HP Omen | Dell |
| Environment | compbio-p0 | compbio-p0-test |
| Operating System | Windows | Windows |
| Python | 3.11.18 | 3.11.18 |
| NumPy | 2.4.6 | 2.4.6 |
| pandas | 3.1.0 | 3.1.0 |
| Matplotlib | 3.11.2 | 3.11.2 |
| JupyterLab | 4.6.4 | 4.6.4 |

## Methodology

The original computational environment was exported into environment.yml and version-controlled using GitHub. The environment was independently reconstructed on the Dell laptop.

An identical numerical analysis was performed on both systems using synthetic gene expression values: 10, 20, 30, 40 and 50.

## Results

| Metric | HP Omen | Dell |
|---|---|---|
| Mean | 30.0 | 30.0 |
| Median | 30.0 | 30.0 |
| Population standard deviation | 14.142135623730951 | 14.142135623730951 |
| Population variance | 200.0 | 200.0 |
| Sum | 150 | 150 |

All five numerical outputs matched exactly.

## Conclusion

The successful reconstruction of the computational environment and identical numerical outputs demonstrate computational reproducibility for the selected analytical test across two independent Windows computers.

## Limitations

The experiment uses synthetic data and a simple deterministic numerical analysis. It does not establish reproducibility for all biological datasets, operating systems or computational workflows.

## Final Status

Computational reproducibility test: PASSED.

## Clean Environment Deletion and Reconstruction

### Date
27 November 2026

### System
HP Omen (Windows)

### Procedure

1. Deactivated the original Conda environment.
2. Removed the compbio-p0 environment completely.
3. Confirmed its removal using conda env list.
4. Reconstructed the environment using conda env create -f environment.yml.
5. Activated the reconstructed compbio-p0 environment.
6. Successfully imported NumPy, pandas and Matplotlib.
7. Verified all five essential software versions.
8. Repeated the original numerical analysis.

### Results

| Parameter | Original Environment | Rebuilt Environment |
|---|---|---|
| Python | 3.11.18 | 3.11.18 |
| NumPy | 2.4.6 | 2.4.6 |
| pandas | 3.1.0 | 3.1.0 |
| Matplotlib | 3.11.2 | 3.11.2 |
| JupyterLab | 4.6.4 | 4.6.4 |
| Mean | 30.0 | 30.0 |
| Median | 30.0 | 30.0 |
| Population standard deviation | 14.142135623730951 | 14.142135623730951 |
| Population variance | 200.0 | 200.0 |
| Sum | 150 | 150 |

All five software versions and all five numerical outputs matched exactly.

### Conclusion

The successful clean reconstruction using environment.yml demonstrates environment reproducibility on HP Omen. Together with the independent reconstruction on Dell, the results support cross-machine computational reproducibility for the selected numerical test on Windows.

### Final Validation Status

- Clean environment deletion: PASSED
- Environment reconstruction: PASSED
- Software version verification: PASSED
- Scientific library import verification: PASSED
- Numerical reproducibility: PASSED
- Slide Step 6: COMPLETED

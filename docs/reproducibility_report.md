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
# A mesh-free hybrid Chebyshev-Tucker tensor format with applications to multi-particle modelling

[![][doi-badge]][doi-link]
[![][arxiv-badge]][arxiv-link]
[![][code-badge]][code-link]

[doi-badge]: https://img.shields.io/badge/DOI-10.48550/arXiv.2503.01696-blue
[doi-link]: https://doi.org/10.48550/arXiv.2503.01696

[arxiv-badge]: https://img.shields.io/badge/arxiv-2503.01696-red
[arxiv-link]: https://arxiv.org/abs/2503.01696

[code-badge]: https://img.shields.io/badge/MATLAB-≥R2022a-blue.svg
[code-link]: https://www.mathworks.com/products/matlab.html

This repository contains the code to reproduce the results of the following paper:

Peter Benner, Boris N. Khoromskij, Venera Khoromskaia and Bonan Sun  
*A mesh-free hybrid Chebyshev-Tucker tensor format with applications to multi-particle modelling*  
arXiv prerint, [arXiv:2503.01696]([arxiv-link]), 2025, https://doi.org/10.48550/arXiv.2503.01696

## Requirements
This code requires MATLAB version R2022a or later and the [Chebfun](https://www.chebfun.org/) package being installed. See the instructions there for details.

## Reproducing the results
After installing the Chebfun, run `reproduce_all.m` to reproduce the figures and tables in the paper. It should take less than 5 minutes on modern laptops. See the table below for a description of the figures. The tables will be printed in the command window.

| File name | Description |
|--------|-------------|
| [`cp_vectors`](figures/cp_vectors.png) | Fig. 1 left |
| [`Pn_cp_vectors`](figures/Pn_cp_vectors.png) | Fig. 1 right |
| [`newton_err_vs_deg_l`](figures/newton_err_vs_deg_l.png) | Fig. 2 left |
| [`newton_err_vs_deg_s`](figures/newton_err_vs_deg_s.png) | Fig. 2 middle |
| [`newton_err_vs_k`](figures/newton_err_vs_k.png) | Fig. 2 right |
| [`ls_coefficients`](figures/ls_coefficients.png) | Fig. 3 |
| [`Pn_original`](figures/Pn_original.png) | Fig. 5 left |
| [`Pn_compress_err_n256`](figures/Pn_compress_err_n256.png) | Fig. 5 middle |
| [`Pn_total_err_n256`](figures/Pn_total_err_n256.png) | Fig. 5 right |
| [`Pn_total_err_n512`](figures/Pn_total_err_n512.png) | Fig. 6 left |
| [`Pn_total_err_n1024`](figures/Pn_total_err_n1024.png) | Fig. 6 middle |
| [`Pn_total_err_n2048`](figures/Pn_total_err_n2048.png) | Fig. 6 right & Fig. 7 left |
| [`Pn_total_err_n2048_lin`](figures/Pn_total_err_n2048_lin.png) | Fig. 7 middle |
| [`Pn_total_err_n2048_nn`](figures/Pn_total_err_n2048_nn.png) | Fig. 7 right |
| [`svals_CU`](figures/svals_CU.png) | Fig. 8 left |
| [`svals_CU_zoom`](figures/svals_CU_zoom.png) | Fig. 8 right |
| [`La_original_Rs11`](figures/La_original_Rs11.png) | Fig. 9 left |
| [`La_compress_err_Rs11`](figures/La_compress_err_Rs11.png) | Fig. 9 middle |
| [`La_total_err_Rs11`](figures/La_total_err_Rs11.png) | Fig. 9 right |
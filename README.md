# Measuring General Associations in Time Series: An Adaptation and Empirical Evaluation of the CODEC Coefficient in Determining Autoregressive Dynamics

> [!IMPORTANT]
>
> The construction of this README is still in progress and we may do some updates.

On this GitHub we provide all the code and files used in the paper *Measuring General Associations in Time Series: An Adaptation and Empirical Evaluation of the CODEC Coefficient in Determining Autoregressive Dynamics (preprint)* by Juan Pablo Montaño & Mario Arrieta-Prieto

## Abstract

<blockquote><sub> Identifying the maximum informative lag to include in an autoregressive model remains an open research problem due to the computational burden of treating it as a hyperparameter, especially in complex models. This study explores model-agnostic association measures, including Pearson, Spearman, and an adaptation of the recently proposed conditional dependence coefficient (CODEC), for guiding lag selection in time series. We adapt and implement the CODEC-based Feature Ordering by Conditional Independence (CODEC-FOCI) algorithm and evaluate its performance through extensive simulations across linear, nonlinear, stationary, nonstationary, seasonal, and heteroskedastic processes. Results show that CODEC outperforms classical correlation-based measures in nonlinear and nonstationary settings, especially for large sample sizes, whereas Pearson performs better in purely linear models. Moreover, the proposed methodology consistently identifies informative lag structures directly from the original series, with differencing and seasonal decomposition generally providing no improvement and often degrading performance. Applications to benchmark datasets confirm that the detected lag structures are consistent with those reported in the literature. These findings highlight CODEC's potential as a practical, model-free tool for exploratory autoregressive lag identification in time series analysis. </sub></blockquote>

## About this repository

This repository contains the following components:

### 📁 Folders

1.  **results:** Contains csv files with the estimations made by each auto-association method.
2.  **plots:** Stores all the plots and tables used in the article.
3.  **Miscellaneous Material**: Includes supplementary materials, mainly from earlier versions of the paper.

### 📄 R Scripts

4.  **Functions.R:** This file defines some functions that are used in the *`Simulation.R`* file.
5.  **Simulation.R**: Monte Carlo simulation. This code constructs each of the time series generating processes and implements the auto-association measures considered for estimating the autocorrelation parameter $p$. The results are saved in the *`results`* folder.
6.  **plots.R**: After running *Simulation.R*, this script generates the plots shown in the paper (output is saved in the *`plots`* folder).

## Next steps:

- If necessary, we may add a link to a Google Drive folder with some csv files or RData to give the option of avoid intensive computation parts.

- In case you wanna cite this GitHub, please use the same citation of the paper.

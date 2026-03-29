# 🧬 Medicago-EMS-Phenomics

![R](https://img.shields.io/badge/R-Language-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Publication Ready](https://img.shields.io/badge/Status-Publication_Ready-orange.svg)

An automated, robust R-based pipeline designed for the analysis of EMS (Ethyl methanesulfonate) mutagenesis experiments in *Medicago polymorpha*. This toolkit provides end-to-end analytical solutions ranging from toxicological dose-response modeling (LD50 estimation) to rigorous statistical enrichment of phenotypic mutation spectra.

This repository contains two core analytical modules that ensure statistical rigor, false-discovery rate (FDR) control, and publication-quality (SCI standard) data visualization.

## 📦 Dependencies

Ensure you have the following R packages installed before running the scripts:

```R
install.packages(c("dplyr", "tidyr", "ggplot2", "drc", "readr", "broom", "tools", "forcats", "ggsci", "scales", "stringr", "conflicted"))

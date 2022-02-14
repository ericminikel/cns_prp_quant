### PrP concentration in the central nervous system: regional variability, genotypic effects, and pharmacodynamic impact

![](/docs/graphical-abstract.png)

This repository holds the data and source code for the following manuscript:

[Mortberg MA, Zhao HT, Reidenbach AG, Gentile JE, Kuhn E, O'Moore J, Dooley PM, Connors TR, Mazur C, Allen SW, Trombetta BA, McManus AJ, Moore MR, Liu J, Cabin DE, Kordasiewicz HB, Mathews J, Arnold SE, Vallabh SM, Minikel EV. **PrP concentration in the central nervous system: regional variability, genotypic effects, and pharmacodynamic impact.** _JCI Insight._ 2022 Feb 8:e156532. doi: 10.1172/jci.insight.156532. PMID: 35133987.](https://doi.org/10.1172/jci.insight.156532).

This manuscript describes the development and validation of a cross-species PrP ELISA, and its application to dissecting prion disease risk factors, mapping regional variability in PrP expression, and demonstrating that pharmacologic PrP lowering in the brain is reflected in cerebrospinal fluid. 

![](/docs/elisa-setup.png)

_All the reagents and materials set up to run the assay._

#### What's here

In this repository, you can:

+ Read the [full text](/docs/mortberg-2021-medrxiv-v2.pdf) of the paper.
+ Download the [full assay protocol](/docs/elisa-protocol.pdf) and [1-page working checklist](/docs/elisa-checklist.pdf) to implement the assay yourself.
+ Browse the study's [raw data](/data). (Note: skyline files from mass spectrometry are deposited in [Panorama](https://panoramaweb.org/mortprpjci2202fz.url) under ProteomeXchange submission [PXD031432](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD031432).)
+ Run our source code to reproduce the figures in the manuscript from the raw data. Just sync the repository and run `Rscript src/cns_prp_manuscript_figures.R`. Dependencies are `tidyverse`, `sqldf`, `reshape2`, `minpack.lm`, and `magick`. You can also browse the other [source files](/src) to see how we fit curves, process data, and so on.


<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.


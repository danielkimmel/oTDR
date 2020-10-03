# Optimal Targeted Dimensionality Reduction

Optimal Targeted Dimensionality Reduction (oTDR) is a software package in MATLAB offering a systematic, principled and statistically rigorous road-map for neural population analysis. oTDR assimilates several aspects of population coding—dimensionality reduction and the separability and stability of the resulting low-dimensional representations—into a single statistical framework that can be applied either cohesively or modularly, including:

* De-mixing low-dimensional representations of task-relevant variables by synthesizing several assumptions—linear regression, orthogonalization, weighting by observation count, and (optional) noise reduction—in a single objective function. (oTDR is based on Targeted Dimensionality Reduction by Mante et al., *Nature*, 2013, which applied subsets of these assumptions, but did so serially, and thus approximated a solution to the original objective.
* Measuring the sensitivity and specificity of any arbitrary representation (whether computed via oTDR or any other de-mixing technique), as well as the relationships between representations, whether of different variables at the same time point (i.e., separability) or of the same variable at different times (i.e., stability).
* Combining the above tools with a rigorous statistical framework to test the significance of these metrics (again, applicable to any de-mixing technique) and rule-out epiphenomenal findings attributable to the population’s intrinsic, low-level features—i.e., dimensionality and temporal smoothness—to which high-dimensional systems are particularly susceptible. oTDR estimates the contribution of these features more accurately—and estimates statistical significance more conservatively—than prior approaches, such as distributing random dimensions isotropically or shuffling across time or conditions, which can account for either the data’s dimensionality or temporal smoothness, respectively, but not both.

For a full description of the methods and their application see:
Kimmel, D. L., Elsayed, G. F., Cunningham, J. P. & Newsome, W. T. [Value and choice as separable and stable representations in orbitofrontal cortex](https://www.nature.com/articles/s41467-020-17058-y). *Nat. Commun.* 11, 1–19 (2020).
*Please cite this paper when using all or part of the oTDR software package*

## Demo
oTDR's features are demonstrated in the provided [script](test/test_oTDR.m), as applied to the [data](https://github.com/danielkimmel/Kimmel_NatComm_2020) from Kimmel et al. (2020).

Note that the separate [data repository](https://github.com/danielkimmel/Kimmel_NatComm_2020) includes both .mat files with the raw data (directory "raw_data") as used for the paper, as well as, for convenience, .mat files containing the processed results (directory "results"). However, the results files have been generated since publication of the paper, and while the values are nearly identical (e.g., we confirmed that all summary statistics are within 0.1 tolerance), any random distributions will inherently be slightly different between the provided files and the results used for the paper, and consequently any values based on those distributions (e.g., p-values and confidence intervals) may differ slightly.

Once obtaining the processed results (either by running oTDR against the raw data via test/test_oTDR.m or by loading the preprocessed results), one can recreate the figures used in Kimmel et al. (2020) by running the provided [plotting code](test/CB_plot_all_figures.m). The plotting code contains an optional code block that will automatically load the preprocessed results from the separate [data repository](https://github.com/danielkimmel/Kimmel_NatComm_2020).

## Notes
* The term "oTDR" is used in two ways, which may lead to some confusion. (1) "oTDR" refers to the overarching set of dimensionality reduction and statistical techniques provided in this package. (2) "oTDR" also refers to the specific implementation of our dimensionality reduction technique that computes *static* regression axes (sRAs), as in Eq 6, Kimmel et al. (2020). For instance, the functions oTDR.m, runOTDR.m, optimizeOTDR.m, etc. are part of this algorithm to compute sRAs. In contrast, "TDR" refers to the implementation that computes *dynamic* regression axes (dRAs), as in Eq 7, Kimmel et al. (2020), and is executed by functions such as TDR.m, runTDR.m, etc.

* This package has been developed so that the vast majority of the tools are written in a general form applicable to any dataset. However, there are a few exceptions that remain in a form specific to the original study. In particular, the stability analysis (e.g., Figure 5, Kimmel et al. (2020)), as implemented in TDRstabilityAnalysis.m, is hard-coded for three task-relevant variables. It should work for *any* three variables (no fewer, no more), but this has not been tested. In future releases, we will update this package so that all functions are in a general form.

## Dependencies
* MATLAB version 9.6.0.1335978 (R2019a) with the following toolboxes (this list is not necessarily comprehensive): 
	* Statistics and Machine Learning Toolbox
	* Parallel Computing Toolbox
* The following publicly available MATLAB packages are bundled with the current package—including manoptToolBox and tensorMaxEntropy (see [packages](packages/))—and are not required. However, problems could arise if you have alternative versions of these packages installed and listed in your MATLAB path before the current package.

## Credits
oTDR was developed by [Daniel L. Kimmel](https://github.com/danielkimmel) and [Gamaleldin F. Elsayed](https://github.com/gamaleldin), who would like to thank John Cunningham (Columbia University), Bill Newsome (Stanford University) and Valerio Mante (University of Zurich) for guidance and support. This work was supported by the Howard Hughes Medical Institute, Stanford Bio-X, Leon Levy Foundation, and NIH Grants T32 GM007365, R25 MH086466, and T32 MH015144 (D.L.K.) and the Columbia Neurobiology and Behavior Program (G.F.E.).

## License
Copyright (c) 2020, Daniel L. Kimmel and Gamaleldin F. Elsayed.

See [license](LICENSE) for permitted use.

*When using all or part of the oTDR software package, please cite:*  
Kimmel, D. L., Elsayed, G. F., Cunningham, J. P. & Newsome, W. T. [Value and choice as separable and stable representations in orbitofrontal cortex](https://www.nature.com/articles/s41467-020-17058-y). *Nat. Commun.* 11, 1–19 (2020).


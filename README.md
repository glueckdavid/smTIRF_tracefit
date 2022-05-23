# smTIRF_tracefit
Shared repository for upload, discussion and troubleshooting of smTIRF scripts.

This ReadMe-file serves as a notepad for project goals, problems and points of discussion.

Project goals:
- Ingestion of single molecule trace files created by the Oxford Nanoimager
- Conversion of traces according to global parameters
- Automated estimation of local parameters (blinking- and bleaching-events)
- Application of local parameters (categorization according to groups, truncation of traces before blinking/bleaching)
- Visualization of datasets (trace-wise with/without applied parameters, global histograms)
- Fitting of step function to individual traces (see github.com/LandesLab for implementation in MatLab)
- Creation of FRET-Histograms with variable bin sizes (see matplotlib documentation for matplotlib.pyplot.hist)
- Fitting of 1-5gaussian functions to histograms
- Saving of output-data for converted traces, local parameters, Trace fitting levels, histogram data, fit data
- Saving of relevant plots with annotation

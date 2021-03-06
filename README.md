# smTIRF_tracefit (Michael Kirchmeyr, David Glück)
Shared repository for upload, discussion and troubleshooting of smTIRF scripts.

This ReadMe-file serves as a notepad for project goals, problems and points of discussion.
Edit and directly commit this file to post any comments or bulletpoints.


Short term goals:
-
- Finish the current script as "Trace import and data clean-up".
- Export traces as csv files for further handling.
- Apply Haar-wavelet transform to single molecule FRET traces.
- Display high frequency noise features as time-dependent curve.


Project goals:
- 
- Ingestion of single molecule trace files created by the Oxford Nanoimager.
- Conversion of traces according to global parameters.
- Automated estimation of local parameters (blinking- and bleaching-events).
- Application of local parameters (categorization according to groups, truncation of traces before blinking/bleaching).
- Visualization of datasets (trace-wise with/without applied parameters, global histograms).
- Fitting of step function to individual traces (see github.com/LandesLab/StaSi for implementation in MatLab).
- Creation of FRET-Histograms with variable bin sizes (see matplotlib documentation for matplotlib.pyplot.hist) both for raw and fitted data.
- Fitting of 1-5gaussian functions to histograms.
- Saving of output-data for converted traces, local parameters, Trace fitting levels, histogram data, fit data.
- Saving of relevant plots with annotation.

Tips/general notes:
-
- Break up code on smaller functions whenever possible. This helps to track down bugs and facilitates patching of bugs/improvement of routines.
- Break up the code in several files for better stability. Output-files in form of csv or json-dumps can be written and read easily and act as checkpoints between scripts.
- Try to rely on as few modules as possble. It is easy to introduce bugs and unexpected behaviour by foreign modules.
- Do not hesistate to delete old code. It is often faster to re-write a section than to fix clunky code. GitHub-branches make this easier, if you fear to lose progress.
- Let loops and conditionals communicate by frequent print-outputs.
- Always ask for help if you get stuck. We are no programmers and learn by sharing and discussing even the smallest details.

Current problems/pre-patch notes:
-
-  Indexing errors after trace truncation.

Points of discussion:
-
-  May be worth to consider breaking up the code at the following points: After application of global/local parameters, after histogram of raw data, after trace fit, after histogram of fit data.

# TreeGrowth

This document explains how the files in this repository should be run in order
to replicate the analyses of the manuscript "Regularization: a new tool for
investigating and predicting tree growth".

Some of the scripts in this repository were originally run on a high-performance
computing cluster because they require a lot of computation. As a result, we do
not recommend anyone try to rerun all these scripts on a single computer. We 
provide instructions below for how to test these scripts without running them in
their entirety. As the output of many scripts is used as input for other
scripts, the outputs of all scripts are already contained in the file. Please 
note that if you change the output of any script, other scripts that depend on
that output may no longer run without error.

### R packages required

To run all the scripts in this repository you will need to have the following
R packages installed on your system:

* tidyverse
* paletteer
* ForestPlotR

*Note: the ForestPlotR package is not available on CRAN and must be downloaded
from GitHub. Also, this analysis was conducted with v0.1.0 of ForestPlotR.
Install with "devtools::install_github("sgraham9319/ForestPlotR@v0.1.0")

### Contents

1. Construct tree neighborhoods
2. Define training and test sets
3. Define common neighbor species
4. Run regularized regression models
5. Run likelihood models
6. Create figures and tables

### 1. Construct tree neighborhoods

Run the script "Scripts/Constructing_neighborhoods.R" to construct neighborhoods
for all trees in the raw data. This script calls the ForestPlot R package, which
must be installed (see instructions above) for this script to run.

### 2. Define training and test sets

Run the script "Scripts/Training_test_specification.R" to define the four
training and four test sets.

### 3. Define common neighbor species

Run the script "Scripts/Specify_common_competitors.R" to create a .csv file 
containing a list of the common neighbor species of each focal species.

### 4. Run regularized regression models

Run the script "Scripts/Regularized_regression_model.R" to apply the
regularized regression model to each training set for each focal species. This
script calls the ForestPlot R package, which must be installed (see instructions
above) for this script to run. Please note that this script takes around 15 
minutes to run on a personal laptop. If you wish to check the script in < 15
minutes, please edit line2 16 and 17 of the script to reduce the number of 
focal species and training sets that are included.

### 5. Run likelihood models

The likelihood models were run on a UNIX-based high-performance computing
cluster because the optimizations require a lot of computation. One particular
function used in the modeling scripts (*mclapply*) will only work on a 
UNIX-based system so the scripts cannot be run in their current form on a 
Windows machine. However, because the optimizations take such a long time to 
run, we do not recommend that any of these scripts be run in their current form
on a single machine.

We instead advise anyone interested in testing or reviewing these scripts to run
a single optimization within one or more scripts. To do this, you will need to
uncomment a few lines of code (remove the # from the beginning of the line) and
then run one block of the code. Instructions for how to do this for each script
are provided in the table below, along with an estimated time for the script to
run on a personal laptop.

Script | Lines to uncomment | Lines to run | Estimated run time
:----- | :----------------- | :----------- | :-----------------
no_comp.R | 111-112 | 1-112 | 8 minutes
no_comp_cv.R | 132-133 | 1-133 | 11 minutes
eq_comp.R | 128-129 | 1-129 | 5 minutes
eq_comp_cv.R | 135-136 | 1-136 | 5 minutes
int_comp.R | 141-142 | 1-142 | 5 minutes
int_comp_cv.R | 138-139 | 1-139 | 4 minutes
ss_comp_ABAM.R | 175-176 | 1-176 | > 30 minutes
ss_comp_cv_ABAM.R | 158-159 | 1-159 | > 30 minutes
ss_comp_CANO.R | 163-164 | 1-164 | > 20 minutes
ss_comp_cv_CANO.R | 155-156 | 1-156 | > 20 minutes
ss_comp_PSME.R | 175-176 | 1-176 | > 20 minutes
ss_comp_cv_PSME.R | 158-159 | 1-159 | > 20 minutes
ss_comp_THPL.R | 163-164 | 1-164 | < 10 minutes
ss_comp_cv_THPL.R | 155-156 | 1-156 | < 10 minutes
ss_comp_TSHE.R | 178-179 | 1-179 | > 30 minutes
ss_comp_cv_TSHE.R | 159-160 | 1-160 | > 30 minutes
ss_comp_TSME.R | 157-158 | 1-158 | 4 minutes
ss_comp_cv_TSME.R | 153-154 | 1-154 | 4 minutes

### 6. Create figures and tables

To recreate the figures and tables included in the manuscript, please run the
scripts "Scripts/Tables_figures_species_summary.R",
"Scripts/Tables_figures_modeling_and_prediction.R" and 
"Scripts/Tables_figures_inference.R". These scripts call a number of function
files, which can be found in "Functions/"

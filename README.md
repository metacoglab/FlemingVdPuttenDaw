# FlemingVdPuttenDaw

This repository contains analysis code for the following paper:

Fleming, van der Putten & Daw “Neural mediators of changes of mind about perceptual decisions”

Anonymised behavioural data files and preprocessed ROI data are included in the repository to enable replication of data analyses and rapid generation of the figures in the paper. The group-level T-maps reported in the paper are available at NeuroVault.

FigureX.m files will reproduce the panels of the respective figure in the paper by loading in .mat and/or .csv files containing relevant behavioural data, regression outputs or preprocessed fMRI data. The paths in these scripts require altering to point to the relevant folder in your local version of the repo. E.g. to run Figure2.m, you would alter:

```
baseDir = '~/Dropbox/Research/Metacognition/stateactionexpt/github/stan/modelfits'; %% path to stan model fits
```

to 

```
baseDir = ‘~/pathToGithubRepo/stan/modelfits'; %% path to stan model fits
```

**Supporting code**

**mri**

Contains scripts for setting up the fMRI GLMs. fMRI data were preprocessed using standard pipelines available from our lab repository MetaLabCore (https://github.com/metacoglab/MetaLabCore), using SPM 12 v. 6225 (www.fil.ion.ucl.ac.uk/spm). 

**stan**

Contains code for fitting the computational models to behavioural data using STAN. See stan_metaConf_group.R.
Fitting the computational models to data requires a STAN installation (http://mc-stan.org).

**stats**

Contains R code for hierarchical regression analysis of behavioural and ROI data used to write out the .csv files containing the regression coefficients plotted in Figure 5, and in the Supplementary Tables. Also stores permutation test outputs for fMRI time course regressions calculated in Figure5.m, to avoid needing to recompute.

We make use of the following packages under R:

lme4 (http://cran.r-project.org/web/packages/lme4)

R.matlab (https://cran.r-project.org/web/packages/R.matlab/index.html)

car (https://cran.r-project.org/web/packages/car/index.html)

doBy (https://cran.r-project.org/web/packages/doBy/index.html)

optimx (https://cran.r-project.org/web/packages/optimx/index.html)

stargazer (https://cran.r-project.org/web/packages/stargazer/index.html)

Implementing the mediation analyses in MATLAB requires the Mediation Toolbox, available here:
https://canlabweb.colorado.edu/wiki/doku.php/help/mediation/m3_mediation_fmri_toolbox

**License**

This code is being released with a permissive open-source license. You should feel free to use or adapt the utility code as long as you follow the terms of the license, which are enumerated below. If you make use of or build on the computational models or behavioural/neuroimaging analyses, we would appreciate that you cite the paper.

Copyright (c) 2017, Stephen Fleming

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Degree Variability Project

Code to compare hub location and distribution across diffusion MRI preprocessing pipelines.
This work has been motivated by papers and code written by [Stuart Oldham](https://scholar.google.com.au/citations?hl=en&user=jj5dZe0AAAAJ&view_op=list_works&sortby=pubdate) and colleagues.

## Installation
1. Clone this repo and add to your MATLAB path
2. Clone the helper functions [here](https://github.com/magnesium2400/nihelp)
3. Add in the dependencies [BCT (the Brain Connectivity Toolbox)](https://sites.google.com/site/bctnet/) and [Stuart Oldham's Brain Surface Visualiser _plotSurfaceROIBoundary_](https://github.com/StuartJO/plotSurfaceROIBoundary)
4. Data can be downloaded from https://doi.org/10.26180/c.6352886.v1 and placed in the DATA folder

## Usage
0. `CalcSTRthr` takes in the raw data and generates the processed data
	1. This step is optional - the processed data can be downloaded
	2. Thresholds can be adjusted here, if desired
	3. New group reconstruction algorithms can be added here, if desired
1. `RunAnalysis` generates summary plots of all the pipelines including heatmaps and surface area correlations
	1. This step is optional
	2. Note that this takes a long time
2. The files in `./SCRIPTS/FigureFunctions` generate the figures from the paper

## Compatibility
The code has been tested on versions of MATLAB from R2019b to R2021b.

## See Also
- [Paper (OA)](https://direct.mit.edu/netn/article/7/4/1326/116174/Can-hubs-of-the-human-connectome-be-identified) available from _Network Neuroscience_
- [Data](https://doi.org/10.26180/c.6352886.v1) available at https://doi.org/10.26180/c.6352886.v1
- This work follows on from work by Stuart Oldham and colleagues - see their [paper](https://www.sciencedirect.com/science/article/pii/S1053811920307382), [code](https://github.com/BMHLab/MotionStructuralConnectivity), and [data](https://doi.org/10.26180/5e7313d012cee) (all openly available)

## Further Details
Please contact Mehul Gajwani by [email](mailto:mehul.gajwani1@monash.edu) or through GitHub (issue/pull request) for any further details.

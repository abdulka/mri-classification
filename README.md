# Classification of structural T1 weighted MRI images.
This repository contains code to replicate the classification procedure used in [1].

## Pre-requisits
The code in this repository depends on the proprietary and free software and data listed below.
  * Matlab R2015a
  * SPM8 (http://www.fil.ion.ucl.ac.uk/spm/software/spm8/)
  * VBM8 (http://www.neuro.uni-jena.de/vbm/download/)
  * libsvm-3.22 (https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
  * LONI Probabilistic Brain Atlas (LPBA) (http://www.loni.usc.edu/atlases/Atlas_Detail.php?atlas_id=12,[2])

## Usage
Classification of MRI images requires four main steps that are described below.
1) Pre-processing (CoregSegment.m)
This steps produces spatially normalized maps of modulated gray matter from the native T1 weighted image.
```matlab
CoregSegmentVBM8('/path/to/t1_image.nii');
```
Consult the help text to learn about further parameters.

2) Feature extraction (extractstructuralfeatures.m)
This step will extract a number of feature sets from a cell array of folders
and corresponding basenames, e.g.:
```
folders = {'./S1' './S2' './S3'};
basenames = {'T1.nii' 'T1.nii' 'T1.nii'};
voxelGM = extractstructuralfeatures(folders,basenames);
```
Note that the LPBA [2] must be installed on the system and the
path currently is hard coded on lines 10,12,14, and 178. This
will be updated in the future. Refer to the help text of the
function for more information.

3) Adjusting for covariates (gprkerneldetrend.m)
The function takes a kernel matrix, the covariates, and an index
vector that indicates which examples to use to estimate the
regression model. See the help text of the function for more details.

4) Classification
Classification can be performed with any kernel-based classifier.
If Q is the adjusted kernel matrix and TR an index vector of the training
examples, then the training and prediction is done as follows:
```
% train
model = svmtrain([(1:sum(TR))' Q(TR,TR)],'-s 0 -t 4 -b 1 -q');
% predict
yhat = svmtrain([(1:sum(~TR))' Q(~TR,TR)],model,'-b 1 -q');
```


## References
[1] Abdulkadir et al. Separating symptomatic Alzheimer’s Disease from depression based on structural MRI (in preparation/review)

[2] Shattuck DW, Mirza M, Adisetiyo V, Hojatkashani C, Salamon G, Narr KL, Poldrack RA, Bilder RM, Toga AW, Construction of a 3D Probabilistic Atlas of Human Cortical Structures, NeuroImage (2007), doi: 10.1016/j.neuroimage.2007.09.031 
# split_tissue_maps

split_tissue_maps estimates sepatated tissue measures given tissue probability maps and input measures using TV optimizations.

## Required Packages
* **NIfTI_20140122**: https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
* **TVAL3_v1.0**: http://www.caam.rice.edu/~optimization/L1/TVAL3/welcome-download.html
  * Download packages and extract them to the same directory.
  * To build matlab binaries, see ```build.sh``` file (tested with version R2012b).

## Command Line Parameters
```
==========================================================================
split_tissue_maps
  Version 1.0
  Developed by Dongjin Kwon
==========================================================================

Required:
  -i <file> : Input image.
  -p <file> : Tissue probabilities. More than one need to be specified, each one preceded by its own -p.
  -o <path> : Output prefix.
Optional:
  -mu <real> : Regualarization parameter (default: 10).
```

# YandaSoft

Yandasoft is a suite of applications and software developed by CSIRO for the calibration and imaging of Interferometric Radio Telescope data. 

## In this version

This release of the software is the first and consequently *beta* release. It contains the (at least) the following applications:

### Measurement Set creation and manipulation
* `csimulator`: simulation of visibilities. 
* `ccontsubtract`: continuum subtraction.
* `mslist`: measurement set interogation.
* `msconcat`: concatenation of measurement sets.
* `msmerge`: merging of measurement sets. 

### Calibration tools
* `cbpcalibrator`: bandpass calibrator.
* `ccalibrator`: for performing gain calibration 
* `ccalapply`: for the application of calibration solutions. 

### Imaging tasks
* `cimager`: Original ASKAP imager.
* `imager`: New imager - permits more parallisation options.
* `linmos`: Linear mosaicking of images


## About Yandasoft

These tasks were originally developed to form the real-time calibration and imaging pipeline for the ASKAP telescope. In order to distribute this software more widely we have extracted these tools from the main codebase of the telescope system and distributed it separately.

### Dependencies

The dependencies are listed in detail in the INSTALL.txt file. But it should be noted that there are internal "ASKAP" dependencies that are required. They are all public and are automatically pulled from their respective repositories by the included `build_all.sh` script.

### How to get it

Yadasoft and its required ASKAP dependencies are available from the CSIRO bitbucket server at:

* https://bitbucket.csiro.au/scm/askapsdp/lofar-common.git 
* https://bitbucket.csiro.au/scm/askapsdp/lofar-blob.git 
* https://bitbucket.csiro.au/scm/askapsdp/base-askap.git 
* https://bitbucket.csiro.au/scm/askapsdp/base-logfilters.git 
* https://bitbucket.csiro.au/scm/askapsdp/base-imagemath.git 
* https://bitbucket.csiro.au/scm/askapsdp/base-scimath.git 
* https://bitbucket.csiro.au/scm/askapsdp/base-askapparallel.git 
* https://bitbucket.csiro.au/scm/askapsdp/base-accessors.git 
* https://bitbucket.csiro.au/scm/askapsdp/yandasoft.git 



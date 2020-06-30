![Image](./common/Figure1.png)


### Demo code for kinematic preprocessing and behavioral analysis
This is demo code for preprocessing raw 3D kinematic data collected using DANNCE or CAPTURE, and then performing behavioral segmentation using t-SNE and watershed clustering. It also includes scripts for visualizing kinematic data.


Main demo files:
```
preprocessing/preprocess_dannce.m
CAPTURE_quickdemo.m
```

### Data Directories
You can download the DANNCE dataset and test CAPTURE dataset from:


### This code uses functions from the following open source software:
- Chronux: https://www.mathworks.com/matlabcentral/fileexchange/68537-chronux-analysis-software
- Pca Randomized: ()
- MTimesX: https://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support
- Motionmapper: https://github.com/gordonberman/MotionMapper
- Structvars: https://www.mathworks.com/matlabcentral/fileexchange/26216-structure-fields-to-variables?focused=8211475&tab=function

**Compatibility**.
This code has been tested in 64-bit MATLAB 2017b and 2019b running on Windows 10.

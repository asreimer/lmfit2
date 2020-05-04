## Levenburg-Marquardt fitting of SuperDARN auto-correlation functions (ACFs). **With NO ad hoc assumptions/conditions!**
Fitted SuperDARN data with accurate error bars. Filter data based on fitted error:

Python version!

# Summary
This repository contains code to that can fit the SuperDARN rawacf files using an error weighted non-linear least-squares fitting algorithm called the First-Principles Fitting Methodology (FPFM). The details of this algorithm as described in [Reimer et. al. (2018)](https://doi.org/10.1002/2017RS006450). The real and imaginary components are fitted against a decaying complex sinusoidal model of the ACF. The Levenburg-Marquardt technique is used to minimize the chi-squared sum. No ad hoc conditions are utilized. The errors in the real and imaginary components of the ACF lag estimates are based on the work by [Reimer et. al. (2016)](https://doi.org/10.1002/2016RS005975). Self-clutter due to the multiple-pulse technique is accounted for using the Mean Power-based Self-clutter Estimator (MPSE) detailed in [Reimer and Hussey (2015)](https://doi.org/10.1002/2015RS005706).

# Python Code
The python version of lmfit2 requires the numpy and lmfit python packages to be installed. The code also relies on the python version of the dmap library [(pydmap written by Keith Kotyk)](https://github.com/SuperDARNCanada/pydmap), which is included in this repository.

The python lmfit package can be installed using something like:
$ pip install lmfit

# Limitations and TO DO
The current version of the C code is not supported for usage with anything but RSTLite. A version compatible the new RST>=4.0 is currently under development.

Currently, this code only fits the ACF not the XCF, so no elevation angles are fit for yet. Also, only the exponential model of the ACF is fit. A version of the code that fits both Gaussian and exponential and ACF and XCF is under development.

The python code may not be updated to fit XCFs nor to fit the Gaussian model. It will depend on time/demand.

# Speed
The python code hasn't yet been speed tested.


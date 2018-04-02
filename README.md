## Levenburg-Marquardt fitting of SuperDARN auto-correlation functions (ACFs). **With NO ad hoc assumptions/conditions!**

# Summary
This repository contains code to that can fit the SuperDARN rawacf files using an error weighted non-linear least-squares fitting algorithm. The real and imaginary components are fitted against an exponential decaying complex sinusoidal model of the ACF. The Levenburg-Marquardt technique is used to minimize the chi2 sum. No ad hoc conditions are utilized. The errors in the real and imaginary components of the ACF lag estimates are based on the work by [Reimer and Hussey (2016)](http://onlinelibrary.wiley.com/doi/10.1002/2016RS005975/full).

# C Code
The C code requires either [RST](https://github.com/vtsuperdarn/VTRST3.5) or [RSTLite](https://github.com/vtsuperdarn/RSTLite) to be installed (need the C dmap library and the hdw.dat files). It also utilizes the [cmpfit library](https://www.physics.wisc.edu/~craigm/idl/cmpfit.html) written by Craig B. Marquardt. The codebase is adapted from the [fitacf3.0 code](https://github.com/SuperDARNCanada/fitacf.3.0) written by Keith Kotyk.

Once the RST dependency is installed, to compile the lmfit2 code, simply change into the C code directory and use the make command. The compiled binary will be in the bin directory.

# Python Code
The python version of lmfit2 requires the numpy and lmfit python packages to be installed. The code also relies on the python version of the dmap library [(pydmap written by Keith Kotyk)](https://github.com/SuperDARNCanada/pydmap), which is included in this repository.

# Limitations and TO DO
Currently, this code only fits the ACF, so no elevation angles are fit for yet. Also, only the exponential model of the ACF is fit.

# Speed
The C code takes approximately 9.5 minutes for 1 hour worth of raw data. Python hasn't yet been speed tested.


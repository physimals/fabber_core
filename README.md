Fabber
------

Bayesian nonlinear model fitting provides a more flexible alternative to traditional linear analysis of perfusion data. In multi-echo or multi-inversion-time data, this approach uses all time points to estimate all the parameters of interest simultaneously. This can be more accurate because it uses a detailed nonlinear model of the pulse sequence, rather than assuming that different parts of the signal are linearly related to each quantity of interest.

To make this analysis technique fast enough for routine use, we have adapted a fast approximate method known as Variational Bayes (VB) to work with non-linear generative models. Experimentally we have found that these calculations provide nearly identical results to sampling method and require only a fraction of the computation time (around a minute per slice).

fabber is a tool for approximate Bayesian inference using generative signal models. The design is modular and new nonlinear forward models can be incorporated into the source code easily. Dual-echo ASL models are documented here, fabber can also be used for multi inversion time resting ASL data, but this can most easily be achieved through the BASIL tool.

**Referencing Fabber**

If you use fabber in your research, please make sure that you reference at least the first of the articles listed below, and ideally the complete list.

*Chappell, M.A., Groves, A.R., Woolrich, M.W., “Variational Bayesian inference for a non-linear forward model”, IEEE Trans. Sig. Proc., 2009, 57(1), 223–236.*

*Woolrich, M., Chiarelli, P., Gallichan, D., Perthen, J., Liu, T. “Bayesian Inference of Haemodynamic Changes in Functional ASL Data”, Magnetic Resonance in Medicine, 56:891-906, 2006.*

*Groves, A. R., Chappell, M. A., & Woolrich, M. W. (2009). Combined spatial and non-spatial prior for inference on MRI time-series. NeuroImage, 45(3), 2009. <doi:10.1016/j.neuroimage.2008.12.027>.*

GIFTI support
=============

This branch has been created to add support for GIFTI input.

What is required for this?

 - [x] Add option --surface=<gii surface file> to specify reference surface file

 - [x] Add FabberRunDataGifti to load/save voxel data from GIFTI files

Currently we are using libfslsurface. This seems to be a bit buggy and under development. A better
option is probably to use libgiftiio directly until it has matured. In particular handling of
timeseries data seems poor and there is ugly debug messages being written to stdout. The main
advantage of libfslsurface is that we can support formats other than GIFTI.

 - [x] Add support in fabber_core.cc to detect --surface option and use FabberRunDataGifti 
       instead of FabberRunDataNewimage

This is ugly - we have to create a temporary FabberRunData object to parse the options initially
and if we find --surface we then re-parse with FabberRunDataGifti. Consequence of the way
we implemented run data which is not ideal for this use case.

 - [x] Move calculation of neighbours into the FabberRunData instance

This is required because knowledge of the connectivity of voxels/vertices is now 
dependent on the input data format - volume or surface.

 - [ ] Implement neighbour calculation for FabberRunDataGifti

Straighforward from vertex/triangle connectivity data

 - [ ] Evaluate whether spatial priors are valid on the assumption that a voxel may have 
   any number of neighbours at potentially varying distances. 
   
This is a difficult question!

 - [x] Evaluate use of FabberRunData.SetExtent which seems kind of useless but is used in the
   C API so it presumably has some purpose. Maybe can be limited to FabberRunDataArray


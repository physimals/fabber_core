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

Theory behind Fabber
====================

Fabber uses a technique called *Variational Bayes* to perform Bayesian inference in a computationally
efficient way. This section gives a brief overview of the theory behind Fabber, its 
advantages and limitations. For a fuller description see [1]_ .

Forward model
-------------

The forward model :math:`M` varies according to application, and produces a predicted *time series*
for a set of parameter values :math:`P_n`:

.. math::

    S(t) = M(t; P_0, P_1, P_2, ...)

A model may have any number of parameters, however we are principally interested in those whose
values we wish to estimate (infer) from the data. For example in a CASL model the bolus duration 
is a parameter of the model which predict the ASL signal but we may choose to regard this as fixed
by the acquisition sequence and not infer it. 

From here we will use the term *parameter* only for parameters of the model whose values we intend
to infer.

'Best Fit' modelling
--------------------

One conventional approach to this problem is to calculate the 'best fit' values of the parameters.
This is done relative to some *cost function* for example the squared difference between the 
model prediction and the actual data in non-linear least-squares (NLLS) fitting. The partial derivatives of the
model prediction for each parameter are calculated numerically and used to change their values
to reduce the cost function. When a minimum has been achieved to within some tolerance, the
resulting parameter values are returned.

Linear modelling
----------------

For some kinds of forward model, a pseudo-linear approach can be used where certain features of the
data such as the mean are assumed to be linearly dependent on a model parameter. For example if the
model predicts a peak, it may have a parameter which (approximately) determines the height of the peak,
another which is related to the peak width, another which affects the time position of the peak
and another which adds a constant signal offset. These parameters can then be related to
measurable features of the data and estimated independently.

The effectiveness of this kind of approach depends on the model and in general it is less reliable
with more complex nonlinear models.

Bayesian inference
------------------

In the Bayesian picture of inference, we always start with some *prior* distribution for 
the value of each parameter. This represents what we know about the parameter's value *before*
we have seen the data we are fitting it to.

The prior may be based on experimental evidence in which case it may have a mean (the accepted
experimental value) and a limited variance (reflecting the range of values obtained in 
previous experiments). This is an example of an *informative* prior.

Alternatively if the value of the parameter could vary by a very wide degree we our prior may
be given some default mean with very large variance that effectively allows it to take on any
possible value. This is an *uninformative* prior. 

Prior distributions should not be informed by the data we are fitting, instead Bayesian inference
calculates a *posterior* distribution for each parameter value which takes into account both
our prior knowledge (if any) and what the data says. Mathematically this is based on Bayes's theorem:

.. math::

    P(parameters \mid data) = \frac{P(data \mid parameters) \, P(parameters)}{P(data)}

Here :math:`P(parameters \mid data)`: is the posterior distribution of the parameter values given
the data we have. :math:`P(parameters)` is the prior distribution of the parameter values. 
:math:`P(data \mid parameters)` is the probability of getting the data we have given a set of
parameter values and is determined from the forward model together with some model of random
noise. This term is known as the *likelihood*. The final term :math:`P(data)` is known as the
*evidence*, and can often be neglected as it only provides a normalization of the
posterior distribution.

So, rather than determining the 'best fit' values of model parameters, the output is a
set of posterior distributions for each parameter which include information
about the predicted mean value, and also its variance. If the variance of 
a parameter's posterior distribution is high this means it's value is not precisely determined 
by the data (i.e. there are a range of probable value which are consistent with the data), whereas
if the variance of the posterior is low the value is well determined by the data (and/or the
prior).

Advantages of the Bayesian approach include:

 - The ability to incorporate prior information into our inference. For example this allows
   us to treat a tissue T1 value as a parameter in the model and constrain the inference by
   the known range of expected T1 values. In conventional model fitting we would either need
   to fix the value of T1 or allow it to vary in any way to fit the data, potentially 
   resulting in physically unlikely parameter values (which nevertheless fit the data well!)

 - Related to the above, we can potentially fit models which are formally 'over specified' 
   by parameters (i.e. where there are different combinations of parameter values which 
   produce identical output). Because of the priors, these different combinations of
   values will not be equally likely and hence we can choose between them.

 - The output includes information about how confident we can be in the parameter values 
   as well as the values themselves.

 - It is possible to incorporate other kinds of prior information into the inference, for
   example how rapidly a parameter may vary spatially.

The 'gold standard' approach to Bayesian inference is the Markov chain Monte Carlo (MCMC)
method where the posterior distribution is determined by sampling the prior distributions
and applying Bayes's theorem. However this is generally too computationally intensive to 
use routinely with problems involving :math:`10^5` - :math:`10^6` voxels as in
fMRI applications.

Variational Bayes
-----------------

The method used in fabber makes a variational approximation which nevertheless is able
to reproduce the results of MCMC closely. One consequence of this is that the prior and
posterior distributions must be of the same type. In the Fabber approach these are
modelled as a multi-variable Gaussian distributions which are characterised by parameter
means, variances and covariances.

The key measure used by Fabber in fitting the model is the ``Free energy`` which is
effectively the negative of the Kullback-Leibler divergence between the approximate multi-variate
Gaussian posterior and the ``true`` posterior distribution. In [1]_ an expression for the
free energy based on this form of the posterior is derived, along with a set of update equations
to maximise the free energy (and hence minimise the KL divergence) by varying 
the parameter means, variances and covariances.

References
----------

.. [1] *Chappell, M.A., Groves, A.R., Woolrich, M.W., "Variational Bayesian
   inference for a non-linear forward model", IEEE Trans. Sig. Proc., 2009,
   57(1), 223–236.*

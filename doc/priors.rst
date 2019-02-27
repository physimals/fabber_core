Priors in Fabber
================

Each parameter in a Fabber model has a *prior* which describes our existing knowledge
of the parameter's value before we see any data.

A model must provide a set of priors for all of it's parameters and in general it is
not good practice to modify them - especially in light of knowledge derived from 
the data as this undermines the Bayesian principles.

Nevertheless there are a number of options that can be set for priors which can
be used reasonably.

Spatial priors
--------------

A spatial prior applies spatial regularization to the parameter so that the spatial
variation in it's value is limited by the information present in the data. This has
the effect of smoothing the parameter map in areas where there is not enough 
information in the data to justify more detail. 

This can be beneficial since it produces smoother parameter maps with clearer 
structure but done in a principled way which treats each parameter independently
and applies a degree of smoothing related to the information in the data.

A spatial prior would be defined as follows::

    --PSP_byname1=myparam
    --PSP_byname1_type=M

The first options specifies which named parameter any additional ``--PSP_byname1_*`` options
refer to. The second option sets the prior type as ``M`` which is the most common 
type of spatial prior. Other supported types are ``m``, ``P`` and `p`.

Spatial priors are normally only applied to a single parameter which is representative
of the overall scale of the data. Since all the parameters are linked in the model, 
the result will generally be that all parameters are smoothed appropriately.

ARD priors
----------

Automatic Relevance Detection (ARD) is a type of prior in which a parameter's value
can 'collapse' to zero if there is not sufficient information in the data to justify
it having a nonzero value. This is useful for parameters which may be relevant only
in certain parts of the data, for example an arterial signal component which only
exists in large arteries.

An ARD prior would be defined as follows::

    --PSP_byname1=myparam
    --PSP_byname1_type=A
    
Image priors
------------

An image prior is a prior whose mean value may be different at each voxel. For example
if the tissue's local T1 value is a model parameter, it may be useful to use a 
T1 map calculated by some independent means (e.g. VFA or MOLLI sequences) to provide the
prior value at each voxel, while still allowing for the possibility of variation.

Image priors can be specified as follows::

    --PSP_byname1=myparam
    --PSP_byname1_type=I
    --PSP_byname1_prec=100

Note that the precision can be specified, this controls how free the model is to
vary the parameter. Choosing a high precision (e.g. 1e6) effectively makes the
image 'ground truth'. In this case we have given a precision of 100 which translates
into a standard deviation of 0.1, allowing *some* variation in the inferred value
but ensuring it will remain close to the image value.

Customizing priors
------------------

.. warning::

    Customizing priors, especially in response to information from the data is 
    opposed to the Bayesian methodology and should not be done unless you have
    good reason!

It is possible to override the model's built-in priors and specify their mean and
precision directly. This is done as follows::

    --PSP_byname1=myparam
    --PSP_byname1_mean=1.5
    --PSP_byname1_prec=0.1

This would set the prior for parameter 'myparam' to have a mean of 1.5 and a precision
of 0.1 (variance=10).

While this is normally discouraged, there are cases where it may be appropriate, for 
example when studying a population whose physiological parameters are known to differ
systematically from the average, or for similar reasons to allow a parameter to vary
more from the 'standard' prior value than the model normally allows.

Parameter transformations
-------------------------

Parameter transformations can be used when the default Gaussian distribution does 
not seem appropriate for a parameter. An example would be a parameter which for 
physical reasons cannot be negative. In this case we might guess that a log-normal
distribution would be more appropriate. This can be handled in Fabber by telling
the core inference engine to work with the log of the parameter value (which is
distributed as a Gaussian) and transform it to the actual value when evaluating
the model.

.. warning::

    Transformations are normally built into the model where they are appropriate.
    Inappropriate transformations can lead to numerical instability and poor
    fitting.

Since transformations are transparent to the model they can be modified as follows::

    --PSP_byname1=myparam
    --PSP_byname1_trans=L

This sets the parameter named ``myparam`` to have a log-transform.

Prior mean/precision and transformations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A natural question is how should the prior mean and variance be modified
when using a transformation. For example suppose we have a parameter representing
a transit time and it's normal prior has a mean of 1.3s and a precision of 5.
Unfortunately this defines a Gaussian which has a significant probability of 
being negative, which is probably not physically reasonable.

We might choose to apply a log-transform to this parameter to avoid this problem. 
But what should the mean and variance of the underlying Gaussian distribution
(i.e. the distribution of the log of the value) be.

We might naively assume that the same transform applies fir the mean, however this is not the
case. If we choose :math:`log(1.3)` as our mean we are modelling the prior as 
a log-normal distribution with a *geometric* mean of 1.3, which is subtly different.




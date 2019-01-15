Building new models
===================

For most applications, a model will need to be constructed. This will
include adjustable parameters which Fabber will then fit.

A complete example model is provided in the ``examples`` subdirectory of
the Fabber source code. This provides an easy template to implement a
new model. In the next section we will go through this example.

A simple example
----------------

To create a new Fabber model it is necessary to create an instance of
the class ``FwdModel``. As an example, we will create a model which fits the
data to a sine function with an amplitude, frequency and phase shift,
plus an optional constant offset:

.. math::
    A\sin(B(t+C)) [+D]

.. note::
    The source code and ``CMakeLists.txt`` file for this example are in the
    fabber source code, in the examples subdirectory. We will assume you
    have this to hand as we go through the process!

First we will create the interface .h file which shows the methods we
will need to implement:

.. code::

    // fwdmodel_sine.h - A simple sine curve fitting model
    #pragma once

    #include "fabber_core/fwdmodel.h"

    #include "newmat.h"

    #include <string>
    #include <vector>

    class SineFwdModel : public FwdModel {
    public:
        static FwdModel* NewInstance();

        SineFwdModel()
            : m_include_offset(false)
        {
        }

        std::string ModelVersion() const;
        std::string GetDescription() const;
        void GetOptions(std::vector<OptionSpec> &opts) const;

        void Initialize(FabberRunData &args);
        void EvaluateModel(const NEWMAT::ColumnVector &params, 
                        NEWMAT::ColumnVector &result, 
                        const std::string &key="") const;
        
    private:
        bool m_include_offset;
        static FactoryRegistration<FwdModelFactory, SineFwdModel> registration;
    };

We have not made our methods virtual, so nobody will be able to create a
subclass of our model. If we wanted this to be the case all the
non-static methods would need to be virtual, and we would need to add a
virtual destructor.

We have also declared a private variable to hold whether we are going to
use the optional offset or not.

We will now implement these methods one by one. Many of them are
straightforward. We start our implementation file as follows

::

    //  fwdmodel_sine.cc - Implements a simple sine curve fitting model
    #include "fwdmodel_sine.h"

    #include "fabber_core/fwdmodel.h"

    #include <math.h>

    using namespace std;
    using namespace NEWMAT;

Here are the first few methods:

::

   FactoryRegistration<FwdModelFactory, SineFwdModel> SineFwdModel::registration("sine");

   FwdModel* SineFwdModel::NewInstance()
   {
       return new SineFwdModel();
   }

The first line here registers our model so that it is known to Fabber by
the name *sine*. The second line is a *Factory method* used so that
Fabber can create a new instance of our model when its name appears on
the command line.

::

   string SineFwdModel::ModelVersion() const
   {
       return "1.0";
   }

    string SineFwdModel::GetDescription() const
    {
        return "Example model which uses a sine function";
    }

We’ve given our model a version number, if we update it at some later stage we should change the number
returned so anybody using the model will know it has changed and what version they have. There's also
a brief description which fabber will return when the user requests help on the model.

::

   static int NUM_OPTIONS = 1;
   static OptionSpec OPTIONS[] =
   {
     {"use-offset", OPT_BOOL, "If True, allow an additional constant offset parameter", OPT_NONREQ, "false"},
     {""}
   };

   void SineFwdModel::GetOptions(vector<OptionSpec> &opts) const
   {
       for (int i = 0; OPTIONS[i].name != ""; i++)
       {
               opts.push_back(OPTIONS[i]);
       }
   }

This is the suggested way to declare the options that your model can
take. It is a little cumbersome when there is only one option, but if
you have many options it will make it clear to see what they are.

The OptionSpec definition contains, in order, the option name, a type
indicator (OPT_STR, OPT_INT, OPT_BOOL, OPT_FILE) which indicates what
kind of data is expected, a human readable description, whether the
option must be specified (OPT_REQ) or not (OPT_NONREQ), and finally the
default value if any.

We have a single non-mandatory Boolean option - whether to allow the
extra constant offset. If not specified, it defaults to false, so not
including an offset.

::

    void SineFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
    {
        params.clear();

        int p=0;
        params.push_back(Parameter(p++, "a", DistParams(1, 1e6), DistParams(1, 1e6)));
        params.push_back(Parameter(p++, "b", DistParams(1, 1e6), DistParams(1, 1e6)));
        params.push_back(Parameter(p++, "c", DistParams(0, 1e6), DistParams(0, 1e6)));
        if (m_include_offset) {
            params.push_back(Parameter(p++, "d", DistParams(0, 1e6), DistParams(0, 1e6)));
        }
    }

The GetParameterDefaults method is quite important. It declares the parameters our
model takes, and their prior and initial posterior distributions. 

The code above declares three parameters ``a``, ``b``, ``c`` and a fourth ``d`` when
the optional offset is included. Each parameter has a name and is passes two ``DistParams``
instances defining the *prior* and *initial posterior* distribution for the parameter. 
The ``DistParams`` instances take two parameters - a mean and a variance.

Priors and Posteriors
~~~~~~~~~~~~~~~~~~~~~
*Priors* are central to Bayesian inference, and describe the extent of our belief about a parameter's
value before we have seen any data. 

For example if a parameter represents the T_1 value of
grey matter in the brain there is a well known range of plausible values. By declaring a
suitable prior we ensure that probabilities are calculated correctly and unlikely values 
of the parameter are avoided unless the data very strongly supports this. 

In our case we have no real prior information, so we are using an *uninformative* prior.
This has a large variance so the model has a lot of freedom in fitting the parameters and 
will try to get as close to matching the data as it can. This is reflected in the high
variance we are using (``1e6``). For the mean values, ``a`` and ``b`` are multiplicative so
it makes sense to give them defaults of ``1`` wherease ``c`` and ``d`` are additive so 
prior means of ``0`` seems more appropriate.

The second ``DistParams`` instance represents the initial *posterior*. This is the starting
point for the optimisation as it tries to find the best values for each parameter. Usually this
does not matter too much and can often be set to be identical to the prior. 

Sometimes it
is helpful to set this to a more restrictive value (lower variance) to avoid numerical 
instability. it is also possible to adjust this on a per-voxel basis - i.e. when the data
being fitted is available. We will not do that here, but it can be useful when fitting, for
example, a constant offset, where we can tell the optimisation to start with a value that 
is the mean of the data. This may help avoid instability and local minima.

In general it is against the spirit of the Bayesian approach to modify the priors on the
basis of the data, and no means are provided to do this. It is possible to modify the priors
on a global basis but this is not encouraged and in general a model should try to provide
good priors that will not need modification.

::

   void SineFwdModel::Initialize(FabberRunData& rundata)
   {
           m_include_offset = rundata.GetBool("use-offset");
   }

The initialize function is called before the model will be used. Its
purpose is to allow the model to set up any internal variables based on
the user-supplied options. Here we simply find out whether the user
asked us to include an offset or not.

::

   void SineFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
   {
       // Check we have been given the right number of parameters
       assert(params.Nrows() == NumParams());
       result.ReSize(data.Nrows());

       for (int i = 1; i <= data.Nrows(); i++)
       {
           double res = params(1) * sin(params(2) * (i - params(3)));
           if (m_include_offset) res += params(4);
           result(i) = res;
       }
   }

Finally the real work! We are given a list of parameter values (params)
and need to produce a time series of predicted data values (result). We
use our sine formula to do this. Note that since we have a discrete time
series, the variable i plays the role of x in the formula. We also only
include the offset parameter if it is configured, and the index for
vectors starts at 1, not zero as you might expect.

Building Fabber with our new model
----------------------------------

The example template comes with build files (``Makefile`` and
``CMakeLists.txt``) and convenience build scripts. To build the model we
simple do:

::

   scripts/build.sh release

or in an FSL development environment we can use the Makefile instead and
just do:

::

   make

Both of these commands produce a new executable fabber_sine which
contains our model. The ``cmake`` build puts this (and other build
files) in the ``build_release`` subdirectory.

The ``cmake`` build also produces a shared library (for example
``libfabber_models_sine.so`` on Linux) which can be loaded into the
generic ``fabber`` executable as follows:

::

   fabber --loadmodels=<path to library>/libfabber_models_sine.so

Both of these should contain the new model, as you can show by using the
``--listmodels`` option

Running our simple example
--------------------------

We will show how to run the model using the GUI tool as it makes the
results more visually obvious. However you can run from the command line
if you prefer.

On creating a new Fabber configuration file we need to ensure that we
are using the correct Fabber executable with our model built in to it.
Once this is set, you will see that the ‘sine’ model appears in the list
of forward models. Clicking on ‘Options’ then shows our single option.
We will leave it unset initially.

.. figure:: /uploads/e50fd8068e61fa5d778e582df8be0ec3/sine_options.PNG
   :alt: sine_options

   sine_options

We now click on ‘Run’ using the nlls inference method and a sample fMRI
data set. After some time, Fabber completes and we can see a list of our
output files. If we make the model fit visible and click on points in
the image data we can see how well our model fits the data.

.. figure:: /uploads/efe7d4f4f9e91e1af70da4b22a803135/modelfit_sine_nlls_no_offset.PNG
   :alt: modelfit_sine_nlls_no_offset

   modelfit_sine_nlls_no_offset

Not very well - because we allowed no offset, the best fit is nearly
linear. Let’s switch that option on, and try again.

.. figure:: /uploads/03e57bbd242d9d8de1ef3c47627d4dac/modelfit_sine_nlls_offset.PNG
   :alt: modelfit_sine_nlls_offset

   modelfit_sine_nlls_offset

That’s a bit better - it’s able to pick up the general variation in the
signal. Although it’s clear that this example is not giving us any
physical information we can see how the model has tried to reproduce the
general shape of the data by varying the allowed parameters.

Changing the example to your own model
--------------------------------------

To implement a single new model, it should be as simple as:

-  Edit the source files, ``Makefile`` and ``CMakeLists.txt`` to change
   references to ``sine`` to the name of your model
-  Rename source files, e.g. \ ``fwdmodel_sine.cc`` ->
   ``fwdmodel_<mymodel>.cc``
-  Add your model options to the options list in the ``.cc`` file
-  Implement the ``Evaluate`` and ``HardcodedInitialDists`` methods for
   your model


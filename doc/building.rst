Building Fabber
===============

Building Fabber using an installed FSL distribution
---------------------------------------------------

**You will need FSL to build Fabber** - it requires a number of
libraries distributed as part of FSL. In addition the Fabber
``Makefile`` is based around the FSL build system.

.. note::
    An additional ``cmake`` based build system also exists
    for use particularly on Windows. We will not describe this
    here.

Setting up an FSL development environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to have your system set up to compile FSL code. If you're already
building other FSL tools from source you've probably already done this,
and can skip this section. Otherwise, run the following commands (you
may want to put this in your ``.profile`` or ``.bash_profile`` script
if you are going to be doing a lot of recompiling)::

   source $FSLDIR/etc/fslconf/fsl-devel.sh
   export FSLDEVDIR=<prefix to install into>
   export PATH=$FSLDEVDIR/bin:$PATH

``FSLDEVDIR`` is an alternate prefix to ``FSLDIR`` which is used to 
store updated code separately from the official FSL release. Most
FSL-based scripts should use code installed in ``FSLDEVDIR`` in preference
to the main FSL release code.

Sometimes this is all you need to do, however often you will need to set
up a link so that FSL knows what compiler flags to use for your platform.
A couple of samples are given below:

OSX
###

.. code::
    
    sudo ln -s $FSLDIR/config/apple-darwin13-llvm6.0 $FSLDIR/config/$FSLMACHTYPE

Ubuntu
######

.. code::
    
    sudo ln -s $FSLDIR/config/linux_64-gcc4.8 $FSLDIR/config/$FSLMACHTYPE

Building ``fabber_core``
~~~~~~~~~~~~~~~~~~~~~~~~

You can probably skip this if you are just building an updated model
library. If you need to recompile the core, however, it should be a case of::

   cd fabber_core
   make install

This approach uses the same build tools as the rest of FSL which is
important on some platforms, notably OSX. It will install the updated
code into whatever prefix you selected as ``FSLDEVDIR``.

Building new or updated model libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Model libraries are distributed separately from the Fabber core.
If you need an updated version of a model library, for example
the ASL model library, you first need to get the source code
for the models library. A number of model libraries are
available in our `Github repositories <https://github.com/ibme-qubic/>`_
all named ``fabber_models_<name>``.

Then to build and install the updated model libraries you would then 
run, for example::

    cd fabber_models_asl
    make install

Adding your own models
----------------------

If you want to create your own model to use with the Fabber core
model fitting engine, see `Building a new model library`_. Once you've
designed and coded your model there are two ways to incorporate
it into the Fabber system:

Adding models directly to the core
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you wish, you can add your own models directly into the Fabber source
tree and build the executable as above. This is not generally
recommended because your model will be built into the core executable, however
it can be the quickest way to get an existing model built in. You will
need to follow these steps:

1. Add your model source code into the fabber_core directory, e.g. 

::

   fabber_core/fwdmodel_mine.cc
   fabber_core/fwdmodel_mine.h

2. Edit ``Makefile`` to add your model to the list of core objects, e.g. 

::

   COREOBJS =  fwdmodel_mine.o noisemodel.o fwdmodel.o inference.o fwdmodel_linear.o fwdmodel_poly.o convergence.o motioncorr.o priors.o transforms.o

3. Run ``make install`` again to build and install a new executable

Creating a new models library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the preferred approach if you want to distribute your new models. A template
for a new model library including a simple sine-function implementation is
included with the Fabber source code in ``fabber_core/examples``. See
`Building a new model library`_ for a full tutorial on this example which includes
how to set up the build scripts.

.. _Building a new model library: models.html



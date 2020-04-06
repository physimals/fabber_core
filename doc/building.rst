Building Fabber
===============

In most cases **you don't need to build Fabber** - executables covering a variety of models are available in FSL (v6.0.1 or later recommended). You might need to use the following instructions, however if you:

 - Need to install updated code which has not yet been released in FSL
 - Want to write your own model or otherwise modify the code

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

First you need to have your system set up to compile FSL code. If you're already
building other FSL tools from source you've probably already done this,
and can skip this section. Otherwise, run the following commands::

   source $FSLDIR/etc/fslconf/fsl-devel.sh
   export FSLDEVDIR=<prefix to install into>
   export PATH=$FSLDEVDIR/bin:$PATH

.. note::
    you may want to put this in your ``.profile`` or ``.bash_profile`` script
    if you are going to be doing a lot of recompiling

``FSLDEVDIR`` is an alternate prefix to ``FSLDIR`` which is used to 
store updated code separately from the official FSL release. You might want 
to set it to something in your home directory, e.g. ``$HOME/fsldev``. Most
FSL-based scripts should use code installed in ``FSLDEVDIR`` in preference
to the main FSL release code.

Setting up compiler/linker flags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    This is not always necessary - it depends on the system you're building on.
    So try skipping to the 'Building...' sections below and come back if you get
    compiler/linker errors.

Firstly, if you're getting errors trying to build you need to do ``make clean``
and (to be safe!) ``rm depend.mk`` before building again with new settings. Also
note that the FSL build system changed in v6.0.3 - here we provide only information
for the new build system.

The relevant settings are in ``$FSLDIR/config/buildSettings.mk``. In particular
you may need to modify ``ARCHFLAGS`` for your system - be careful as there are
separate definitions for Linux and Mac ('Darwin') systems, so make sure you
change the right one!

For recent versions of Ubuntu, you need to turn off use of the C++11 ABI as FSL libraries are not
compiled using this. To do this add the following to ``ARCHFLAGS``::

    -D_GLIBCXX_USE_CXX11_ABI=0 -no-pie

If you are having difficulty with other systems, please raise an issue and we will
investigate.

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



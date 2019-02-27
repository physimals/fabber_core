Building Fabber
===============

I have an official FSL distribution installed
---------------------------------------------

You can build Fabber using the FSL build system. First you need to set
up your development environment:

::

   source $FSLDIR/etc/fslconf/fsl-devel.sh
   export FSLDEVDIR=<prefix to install into>

Then building Fabber should be a case of:

::

   cd fabber_core
   make

This approach uses the same build tools as the rest of FSL which is
important on some platforms, notably OSX.

I don't have an FSL distribution which supports development
-----------------------------------------------------------

Fabber has an alternative build system using ``cmake`` as its build
tool. CMake is cross-platform and is designed for out-of-source builds,
so you create a separate build directory and all the compiled files end
up there. This helps to keep your source tree free of build artifacts.

Some convenience scripts are provided to build and install Fabber. You
need to ensure that ``FSLDIR`` is set to point to wherever the FSL
dependencies are installed, then for example to do a build which
includes debug symbols on Linux or OSX you run:

::

   scripts/build.sh debug

The executables and libraries will end up in ``build_debug``. Other
options are:

::

   scripts/build.sh release

   scripts/build.sh relwithdebinfo

The last example is quite useful - it produces a release build with full
optimization (typically about 2x faster than a debug build) but keeps
the debug symbols in the executable so debugging is possible.

You can also look at the scripts - they are very simple - to see the
actual CMake commands being run.

Pitfalls
~~~~~~~~

gcc vs Clang on OSX
^^^^^^^^^^^^^^^^^^^

OSX has Apple’s ``Clang`` compiler and the LLVM C++ standard library
``libc++`` as its default for compiling C++. However FSL uses ``gcc``
and the GNU ``libstdc++`` library for OSX builds. These are not
compatible. If you are building against a standard FSL installation you
need to force CMake to use ``gcc`` and ``libstdc++``. To do this:

1. Uncomment the following line in CMakeLists.txt

   set(CMAKE_CXX_FLAGS “${CMAKE_CXX_FLAGS} -std=c++11
   -stdlib=libstdc++”)

2. Add the following directives to the ``cmake`` command itself (edit
   the scripts/build.sh script if you are using it):

   -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++

Adding your own models
~~~~~~~~~~~~~~~~~~~~~~

If you wish, you can add your own models directly into the Fabber source
tree and build the executable as above. This is not generally
recommended because your model will be built into libfabbercore, however
it can be the quickest way to get an existing model built in. You will
need to follow these steps:

1. Add your model source code into the fabber_core directory,
   e.g. \ ``fabber_core/fwdmodel_mine.cc`` and
   ``fabber_core/fwdmodel_mine.h``

2. Edit CMakeLists.txt to add your model sources as follows::

   # Core objects - things that implement the framework for inference
   set(CORE_SRC noisemodel.cc fwdmodel.cc inference.cc)

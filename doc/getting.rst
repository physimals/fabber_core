Getting Fabber
==============

From FSL
--------

Fabber is distributed as part of  `FSL <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki>`_,
and this is the easiest way to get Fabber if you want to use existing
models for fMRI data . This documentation describes the 
version of Fabber included with *FSL v6.0.1 and above*.

Addition tools that can use Fabber will work with a correctly installed
FSL distribution although currently not all models we have developed are
available in FSL.

Standalone Fabber distribution
------------------------------

Standalone versions of Fabber including a selection of model libraries are available
for a number of platforms. These may be useful if you don't want the rest of FSL
or if you need a more up to date version of Fabber than the one included with FSL.

The current standalone release can be found at https://github.com/physimals/fabber_core/releases.

The standalone release can be used with tools requiring a Fabber installation such as 
the `Python API <https://pyfab.readthedocs.io/>`_, or Fabber-based plugins 
for `Quantiphyse <https://quantiphyse.readthedocs.io/>`_. You should set the environment
variable ``FABBERDIR`` to the unpacked distribution directory to ensure these tools
can find Fabber. Note that some Quantiphyse plugins require a full FSL installation,
notably the ASL plugin.

Building from source code
-------------------------

You can build Fabber from the source code available in the `Github repository <https://github.com/physimals/fabber_core.git>`_.
You will need an FSL installation for this. For instructions see `Building Fabber`_.

.. _Building Fabber: building.html

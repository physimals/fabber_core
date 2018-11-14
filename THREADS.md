Ideas for multithreading
========================

Difficulties
------------

Model is not threadsafe - pass data in before calling Evaluate
vectors in ThreadContext are not threadsafe when written to

Implementation
--------------

ThreadContext class to encapsulate the run data for a thread.
Contains:
  - All mutable data (e.g. output posterior, priors)
  - Non-mutable per-thread data:
    - Timeseries data, suppdata etc (Could be shared but easier to split up on per-voxel basis)
  - Model reference (model is not threadsafe so need 1 per thread)

ThreadContext to have a `Run()` method which takes no parameters.

Specialized ThreadContext implementations for each inference method:

  - NllsThreadContext
  - VbThreadContext
  - SpatialVbThreadContext

InferenceMethod may create ThreadContexts and run them as they see fit - e.g. SpatialVB needs to be more fine grained about threading to ensure spatial updates occur across the
threads.

`fabber_threads.h` contains very basic cross platform threading wrapper. Includes 
a 'no threads' option in which creating a thread simply runs the method synchronously
and joining a thread does nothing.

Spatial mode
------------
In spatial mode, updates to priors are not multi-threaded
Must be distributed to threads after each iteration. Suggests creating 
separate set of threads for each spatial iteration and updating the priors
in the main thread between each iteration. 

SaveResults
-----------
A bit tricky - need to combine output of multiple ThreadContext objects. Hacked 
for now but should be systematic way of doing this.

Alternatives
------------
Coarse parallelism based around multiple RunData instances each with own copy of
InferenceTechnique, FwdModel etc.
Disadvantage here is it will not work for spatial VB
Essentially the same approach as the Python based parallelism in Quantiphyse

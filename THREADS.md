Ideas for multithreading
========================

Difficulties
------------

Model is not threadsafe - pass data in before calling Evaluate
vectors in ThreadContext are not threadsafe when written to

Possible implementation
-----------------------

Put all mutable data into ThreadContext
Also non-mutable per-thread data:
  - Timeseries data, suppdata etc (Could be shared but easier to split up on per-voxel basis)
  - Model reference (model is not threadsafe so need 1 per thread)

Create separate ThreadContext for each thread containing the relevant data subset
Voxel processing only writes to ThreadContext
Copy of the model for each ThreadContext
ThreadContext(s) created by InferenceMethod - necessary since parallelism is method dependent

Spatial mode
------------
In spatial mode, updates to priors are not multi-threaded
Must be distributed to threads after each iteration. Suggests creating 
separate set of threads for each spatial iteration

SaveResults
-----------
A bit tricky - need to combine output of multiple ThreadContext objects

Alternatives
------------
Coarse parallelism based around multiple RunData instances each with own copy of
InferenceTechnique, FwdModel etc.
Disadvantage here is it will not work for spatial VB
Essentially the same approach as the Python based parallelism in Quantiphyse

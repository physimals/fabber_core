Ideas for multithreading
========================

Difficulties
------------

Model is not threadsafe - pass data in before calling Evaluate
vectors in RunContext are not threadsafe when written to

Possible implementation
-----------------------

Put all data into RunContext (including timeseries data, suppdata etc)
Create separate RunContext for each thread containing the relevant subset
Ensure that all voxel processing only writes to RunContext
Either separate copies of the model or locking around passmodeldata/evaluate

Spatial mode
------------
In spatial mode, updates to priors are not multi-threaded
Must be distributed to threads after each iteration. Suggests creating 
separate set of threads for each spatial iteration

Alternatives
------------
Coarse parallelism based around multiple RunData instances each with own copy of
InferenceTechnique, FwdModel etc.
Disadvantage here is it will not work for spatial VB
Essentially the same approach as the Python based parallelism in Quantiphyse

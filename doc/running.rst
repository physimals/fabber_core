Running Fabber
==============

Specifying options
------------------

The simplest way to run Fabber is as a command line program. It uses
the following syntax for options::

    -–option
    -–option=value

A simple example command line would be::

    fabber --data=fmri_volume.nii.gz --mask=roi.nii.gz \
           --model=poly --degree=2 \
           --method=vb --noise=white \
           --output=out --save-model-fit \

Common options
--------------

--output=OUTPUTDIR  
        Directory for output files (including logfile)
--method=METHOD     
        Use this inference method
--model=MODEL       
        Use this forward model
--data=DATAFILE     
        Specify a single input data file
--mask=MASKFILE     
        Mask file. Inference will only be performed where mask value > 0
--optfile
        File containing additional options, one per line, in the same form as specified on the command line
--overwrite
        If set will overwrite existing output. If not set, new output directories will be created by appending '+' to the directory name 
--suppdata
        'Supplemental' timeseries data, required for some models

Selecting what data to output
-----------------------------

--save-model-fit
        Output the model prediction as a 4d volume

--save-residuals
        Output the residuals (difference between the data and the model prediction)

--save-model-extras
        Output any additional model-specific timeseries data

--save-mvn
        Output the final MVN distributions.

--save-mean
        Output the parameter means.

--save-std
        Output the parameter standard deviations.

--save-zstat
        Output the parameter Zstats.

--save-noise-mean
        Output the noise means.

--save-noise-std
        Output the noise standard deviations. 

--save-free-energy
        Output the free energy, if calculated. 

Help and usage information
--------------------------

--help
        Print this usage method. If given with --method or --model, display relevant method/model usage information

--version
        Print version identifier. If given with --model will print the model's version identifier
        
--listmethods
        List all known inference methods

--listmodels
        List all known forward models

--listparams
        List model parameters (requires model configuration options to be given)

--listoutputs
        List additional model outputs (requires model configuration options to be given)


Advanced options
----------------

--simple-output
        Instead of usual standard output, simply output series of lines each giving progress as percentage

--data1=DATAFILE, --data2=DATAFILE
        Specify multiple data files for n=1, 2, 3...

--data-order
        If multiple data files are specified, how they will be handled: concatenate = one after the other,  interleave = first record from each file, then  second, etc.

--mt1=INDEX, --mt2=INDEX
        List of masked time points, indexed from 1. These will be ignored in the parameter updates

--debug
        Output large amounts of debug information. ONLY USE WITH VERY SMALL NUMBERS OF VOXELS

--link-to-latest
        Try to create a link to the most recent output directory with the prefix _latest

--loadmodels
        Load models dynamically from the specified filename, which should be a DLL/shared library

Variational Bayes options (used when method=vb)
-----------------------------------------------

--noise=NOISE
        Noise model to use (white or ar1)

--convergence=CONVERGENCE
        Name of method for detecting convergence - default maxits, other values are fchange, trialmode

--max-iterations=NITS
        number of iterations of VB to use with the maxits convergence detector

--min-fchange=FCHANGE
        When using the fchange convergence detector, the change in F to stop at

--max-trials=NTRIALS
        When using the trial mode convergence detector, the maximum number of trials after an initial reduction in F

--print-free-energy
        Output the free energy in the log file

--continue-from-mvn=MVNFILE
        Continue previous run from output MVN files

--output-only
        Skip model fitting, just output requested data based on supplied MVN. Can only be used with continue-from-mvn

--noise-initial-prior=MVNFILE
        MVN of initial noise prior

--noise-initial-posterior=MVNFILE
        MVN of initial noise posterior

--noise-pattern=PATTERN
        repeating pattern of noise variances for each point (e.g. 12 gives odd and even data points different variances)

--PSP_byname1=PARAMNAME, --PSP_byname2=PARAMNAME
        Name of model parameter to use for prior specification 1, 2, 3... 

--PSP_byname1_type=PRIORTYPE
        Type of prior to use for parameter 1 - I=image prior

--PSP_byname1_image=FILENAME
        File containing image for image prior for parameter 1

--PSP_byname1_prec
        Precision to apply to image prior for parameter 1

--PSP_byname1_transform
        Transform to apply to parameter 1

--allow-bad-voxels
        Continue if numerical error found in a voxel, rather than stopping

--ar1-cross-terms=TERMS
        For AR1 noise, type of cross-linking (dual, same or none)

--spatial-dims=NDIMS
        Number of spatial dimensions (1, 2 or 3). Default is 3.

--spatial-speed=SPEED
        Restrict speed of spatial smoothing

--param-spatial-priors=PRIORSTR
        Type of spatial priors for each parameter, as a sequence of characters. N=nonspatial, M=Markov random field, P=Penny, A=ARD

--locked-linear-from-mvn=MVNFILE
        MVN file containing fixed centres for linearization

Model-specific options
----------------------

These are usually quite extensive and control the fine details of the model that is being
implemented. For example the generic ASL model will need to be told the TIs/PLDs of the
sequence, the number of repeats, the structure of the data, bolus duration and what components
to include in the model (arterial as well as tissue, dispersion and exchange options, ...).

The best way to look at model options is to use ``--help``, e.g.::

    fabber_asl --help --model=aslrest

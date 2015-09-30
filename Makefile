include ${FSLCONFDIR}/default.mk

PROJNAME = fabbercore

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST}
#-I${HOME}/include 
#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__FABBER_MOTION
#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__FABBER_LIBRARYONLY
#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__FABBER_LIBRARYONLY -D__FABBER_LIBRARYONLY_TESTWITHNEWIMAGE
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz 

#
# Executables
#

XFILES = mvntool #remove fabber from here, since we do not ususally build a version of fabber without models
SCRIPTS = fabber_var


# Sets of objects separated into logical divisions
# Basic objects - things that have nothing directly to do with inference
BASICOBJS = tools.o dataset.o dist_mvn.o easylog.o easyoptions.o 
# Core objects - things that implement the framework for inference
COREOBJS =  noisemodel.o fwdmodel.o inference.o utils.o fwdmodel_linear.o

# Infernce methods
INFERENCEOBJS = inference_spatialvb.o inference_vb.o inference_nlls.o

# Noise models
NOISEOBJS = noisemodel_white.o noisemodel_ar.o

# Configuration
CONFIGOBJS = setup.o

# Client objects - contains the main interface
CLIENTOBJS = fabber_core.o

# Everything together
OBJS = ${BASICOBJS} ${COREOBJS} ${INFERENCEOBJS} ${NOISEOBJS} ${CONFIGOBJS}

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

#
# Build
#

all:	${XFILES} libfabbercore.a

mvntool: ${BASICOBJS} mvntool.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BASICOBJS} mvntool.o ${LIBS}

#
# Build a fabber exectuable, this will have nothing but the linear model so it not practically useful for data analysis
#
fabber: ${OBJS} fabber.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} fabber.o ${LIBS}

#
# Library build
#

# fabber core is the basic implementation with no models (except linear), but includes default inference and noise models
libfabbercore.a : ${BASICOBJS} ${COREOBJS} ${INFERENCEOBJS} ${NOISEOBJS} ${CONFIGOBJS} ${CLIENTOBJS}
	${AR} -r $@ ${BASICOBJS} ${COREOBJS} ${INFERENCEOBJS} ${NOISEOBJS} ${CONFIGOBJS} ${CLIENTOBJS}

# DO NOT DELETE

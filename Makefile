include ${FSLCONFDIR}/default.mk

PROJNAME = fabbercore

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L/lib64

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz -ldl

#
# Executables
#

XFILES = fabber mvntool
SCRIPTS = fabber_var

# Sets of objects separated into logical divisions

# Basic objects - things that have nothing directly to do with inference
BASICOBJS = tools.o rundata.o dist_mvn.o easylog.o fabber_capi.o version.o dist_gamma.o rundata_array.o

# Core objects - things that implement the framework for inference
COREOBJS =  noisemodel.o fwdmodel.o inference.o fwdmodel_linear.o fwdmodel_poly.o convergence.o motioncorr.o 

# Infernce methods
INFERENCEOBJS = inference_vb.o inference_nlls.o inference_spatialvb.o

# Noise models
NOISEOBJS = noisemodel_white.o noisemodel_ar.o

# Configuration
CONFIGOBJS = setup.o factories.o

# Library for executables
EXECOBJS = rundata_newimage.o fabber_core.o

# Executable main object
CLIENTOBJS =  fabber_main.o

# Everything together
OBJS = ${BASICOBJS} ${COREOBJS} ${INFERENCEOBJS} ${NOISEOBJS} ${CONFIGOBJS}

# For debugging:
OPTFLAGS = -ggdb -Wall

# Pass Git revision details
GIT_SHA1:=$(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

all:	${XFILES} libfabbercore.a libfabberexec.a

mvntool: ${OBJS} mvn_tool/mvntool.o rundata_newimage.o 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} mvn_tool/mvntool.o rundata_newimage.o ${LIBS}

#
# Build a fabber exectuable, this will have nothing but the linear model so it not practically useful for data analysis
#
fabber: ${OBJS} ${EXECOBJS} ${CLIENTOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} ${EXECOBJS} ${CLIENTOBJS} ${LIBS} 

#
# Library build
#
libfabbercore.a : ${OBJS}
	${AR} -r $@ ${OBJS}

#
# Library build
#
libfabberexec.a : ${EXECOBJS} 
	${AR} -r $@ ${EXECOBJS} 

# DO NOT DELETE

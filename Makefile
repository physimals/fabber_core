include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_core

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -DFABBER_SRC_DIR="\"${PWD}\"" -DFABBER_BUILD_DIR="\"${PWD}\""
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L/lib64

FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
  LIB_NEWMAT = ${LIB_NEWMAT} -lnewmat
else 
  NIFTILIB = -lNewNifti
endif

LIBS = -lutils -lnewimage -lmiscmaths -lprob -L${LIB_NEWMAT} ${NIFTILIB} -lznz -lz -ldl
TESTLIBS = -lgtest -lpthread

#
# Executables
#

XFILES = fabber mvntool
SCRIPTS = fabber_var

# Sets of objects separated into logical divisions

# Basic objects - things that have nothing directly to do with inference
BASICOBJS = tools.o rundata.o dist_mvn.o easylog.o fabber_capi.o version.o dist_gamma.o rundata_array.o

# Core objects - things that implement the framework for inference
COREOBJS =  noisemodel.o fwdmodel.o inference.o fwdmodel_linear.o fwdmodel_poly.o convergence.o motioncorr.o priors.o transforms.o

# Infernce methods
INFERENCEOBJS = inference_vb.o inference_nlls.o covariance_cache.o

# Noise models
NOISEOBJS = noisemodel_white.o noisemodel_ar.o

# Configuration
CONFIGOBJS = setup.o factories.o

# Library for executables
EXECOBJS = rundata_newimage.o fabber_core.o

# Executable main object
CLIENTOBJS =  fabber_main.o

# Unit tests
TESTOBJS = test/fabbertest.o test/test_inference.o test/test_priors.o test/test_vb.o test/test_convergence.o test/test_commandline.o test/test_rundata.o

# Everything together
OBJS = ${BASICOBJS} ${COREOBJS} ${INFERENCEOBJS} ${NOISEOBJS} ${CONFIGOBJS}

# For debugging:
#OPTFLAGS = -ggdb -Wall

# Pass Git revision details
GIT_SHA1:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

# Targets

all:	${XFILES} libfabbercore.a libfabberexec.a

clean:
	${RM} -f /tmp/fslgrot *.o mvn_tool/*.o *.a *.exe core depend.mk fabber_test

mvntool: ${OBJS} mvn_tool/mvntool.o rundata_newimage.o 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} mvn_tool/mvntool.o rundata_newimage.o ${LIBS}

# Build a fabber exectuable, this will have nothing but the generic models so it not practically useful for data analysis
fabber: ${OBJS} ${EXECOBJS} ${CLIENTOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} ${EXECOBJS} ${CLIENTOBJS} ${LIBS} 

# Library build
libfabbercore.a : ${OBJS}
	${AR} -r $@ ${OBJS}

libfabberexec.a : ${EXECOBJS} 
	${AR} -r $@ ${EXECOBJS} 

# Unit tests
test: ${OBJS} ${EXECOBJS} ${CLIENTOBJS} ${TESTOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${TESTINC} -o fabber_test ${OBJS} ${EXECOBJS} ${TESTOBJS} ${LIBS} ${TESTLIBS} 

# DO NOT DELETE

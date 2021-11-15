include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_core

USRINCFLAGS = -DFABBER_SRC_DIR="\"${PWD}\"" -DFABBER_BUILD_DIR="\"${PWD}\""

# The FSL build system changed
# substantially in FSL 6.0.6
# FSL >= 6.0.6
ifeq (${FSL_GE_606}, true)
  LIBS = -lfsl-newimage -lfsl-miscmaths -lfsl-utils \
         -lfsl-cprob -lfsl-NewNifti -lfsl-znz -ldl
# FSL <= 6.0.5
else
  ifeq ($(shell uname -s), Linux)
	MATLIB := -lopenblas
  endif

  USRINCFLAGS += -I${INC_NEWMAT} -I${INC_CPROB} -I${INC_BOOST} \
                 -I${FSLDIR}/extras/include/armawrap
  USRLDFLAGS   = -L${LIB_NEWMAT} -L${LIB_PROB}         \
                 -lnewimage -lmiscmaths -lutils -lprob \
                 -lNewNifti ${MATLIB} -lznz -lz -ldl
endif

TESTLIBS = -lgtest -lpthread

#
# Executables and libraries provided by this project
#

XFILES  = fabber mvntool
SOFILES = libfsl-fabbercore.so libfsl-fabberexec.so
AFILES  = libfabbercore.a libfabberexec.a
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
GIT_SHA1 := $(shell git describe --dirty)
GIT_DATE := $(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

clean:
	${RM} -f /tmp/fslgrot *.o mvn_tool/*.o *.a *.so *.exe core depend.mk fabber_test

# FSL >=606 uses dynamic linking
ifeq (${FSL_GE_606}, true)

all: ${XFILES} ${SOFILES}

# Dynamically linked libraries are compiled for FSL >= 606
libfsl-fabbercore.so : ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

libfsl-fabberexec.so : ${EXECOBJS} | libfsl-fabbercore.so
	${CXX} ${CXXFLAGS} -shared -o $@ $^ -lfsl-fabbercore ${LDFLAGS}

mvntool: mvn_tool/mvntool.o rundata_newimage.o | libfsl-fabbercore.so
	${CXX} ${CXXFLAGS} -o $@ $^ -lfsl-fabbercore ${LDFLAGS}

# Build a fabber exectuable, this will have nothing but the generic models so it not practically useful for data analysis
fabber: ${CLIENTOBJS} | libfsl-fabberexec.so libfsl-fabbercore.so
	${CXX} ${CXXFLAGS} -o $@ $^ -lfsl-fabberexec -lfsl-fabbercore ${LDFLAGS}

# Unit tests
test: ${CLIENTOBJS} ${TESTOBJS} | libfsl-fabberexec.so libfsl-fabbercore.so
	${CXX} ${CXXFLAGS} -o fabber_test $^ -lfsl-fabberexec -lfsl-fabbercore ${LDFLAGS} ${TESTLIBS}

# FSL <=605 uses static linking
else

all: ${XFILES} ${AFILES}

libfabbercore.a : ${OBJS}
	${AR} -r $@ $^

libfabberexec.a : ${EXECOBJS}
	${AR} -r $@ $^

mvntool: ${OBJS} mvn_tool/mvntool.o rundata_newimage.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

fabber: ${OBJS} ${EXECOBJS} ${CLIENTOBJS}
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

test: ${OBJS} ${EXECOBJS} ${CLIENTOBJS} ${TESTOBJS}
	${CXX} ${CXXFLAGS} -o fabber_test $^ ${LDFLAGS} ${TESTLIBS}
endif

# DO NOT DELETE

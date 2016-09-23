AUTOMAKE_OPTIONS=foreign

bin_PROGRAMS=fabber
check_PROGRAMS=testfabber
lib_LTLIBRARIES=libfabbercore.la

TESTS=testfabber

# Basic objects - things that have nothing directly to do with inference
BASIC_SRC = tools.cc dataset.cc dist_mvn.cc easylog.cc setup.cc fabber_core.cc

# Core objects - things that implement the framework for inference
CORE_SRC =  noisemodel.cc fwdmodel.cc inference.cc utils.cc fwdmodel_linear.cc fwdmodel_trivial.cc fwdmodel_poly.cc convergence.cc

# Inference methods
INFERENCE_SRC= inference_vb.cc inference_nlls.cc inference_spatialvb.cc

# Noise models
NOISE_SRC = noisemodel_white.cc noisemodel_ar.cc

CPPFLAGS=-I$(FSLDIR)/include -I$(FSLDIR)/extras/include/newmat -I$(FSLDIR)/extras/include/boost -DUSE_NEWIMAGE

fabber_SOURCES=fabber_main.cc 
fabber_CPPFLAGS= -I. $(CPPFLAGS)
fabber_LDADD= -lfabbercore \
	-L$(FSLDIR)/lib -L$(FSLDIR)/extras/lib \
	-lnewimage -lmiscmaths -lnewmat -lfslio -lniftiio -lutils \
	-lpthread -lznz -lz
#fabber_DEPENDENCISE: libfabbercore.a
#fabber_LDFLAGS=--all-static

testfabber_SOURCES=test/fabbertest.cc \
                   test/test_inference.cc \
                   test/test_spatialvb.cc \
                   test/test_vb.cc \
                   test/test_convergence.cc \
		   test/test_rundata.cc
testfabber_CPPFLAGS= -I. $(CPPFLAGS)
testfabber_LDADD= $(fabber_LDADD) -lgtest 
#testfabber_LDFLAGS=--all-static
	
libfabbercore_la_SOURCES=$(BASIC_SRC) $(CORE_SRC) $(INFERENCE_SRC) $(NOISE_SRC)
libfabbercore_la_CPPFLAGSas=$(CPPFLAGS)
libfabbercore_la_LIBADD=
	
all-local :
	$(DOXYGEN)

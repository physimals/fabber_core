include ${FSLCONFDIR}/default.mk

PROJNAME = fabber

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__DEVEL #*** NB these (D) tags should become redundant with the new arrangements
#-I${HOME}/include 
#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__FABBER_MOTION
#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__FABBER_LIBRARYONLY
#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -D__OXASL -D__FABBER_LIBRARYONLY -D__FABBER_LIBRARYONLY_TESTWITHNEWIMAGE
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

#LIBS = -lutils -lprob -lnewmat # Will report the MISCMATHS dependencies
#LIBS = -lutils -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz 
LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz 

#
# Executables
#

XFILES = fabber mvntool
SCRIPTS = fabber_var

#
# C++ build - original all-in-one approach
#

# Sets of objects separated into logical divisions
# Basic objects - things that have nothing directly to do with inference
BASICOBJS = tools.o dataset.o dist_mvn.o easylog.o easyoptions.o 
# Core objects - things that implement the framework for inference
COREOBJS =  noisemodel.o fwdmodel.o inference.o utils.o fwdmodel_linear.o
# Forward model groups
FWDOBJS_ASL =  fwdmodel_asl_multiphase.o fwdmodel_asl_grase.o asl_models.o fwdmodel_asl_rest.o
FWDOBJS_DUALECHO = fwdmodel_quipss2.o fwdmodel_q2tips.o fwdmodel_pcASL.o
# ***TEMP removed models: fwdmodel_custom.o  fwdmodel_simple.o fwdmodel_asl_buxton.o fwdmodel_asl_pvc.o fwdmodel_asl_satrecov.o fwdmodel_cest.o fwdmodel_dsc.o fwdmodel_dce.o fwdmodel_biexp.o
# all the forward models
FWDOBJS = ${FWDOBJS_ASL} ${FWDOBJS_DUALECHO}
# Infernce methods
INFERENCEOBJS = inference_spatialvb.o inference_vb.o inference_nlls.o
# Noise models
NOISEOBJS = noisemodel_white.o noisemodel_ar.o
# Configuration
CONFIGOBJS = setup.o

# Everything together
OBJS = ${BASICOBJS} ${COREOBJS} ${FWDOBJS} ${INFERENCEOBJS} ${NOISEOBJS} ${CONFIGOBJS}

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

all:	${XFILES}

fabber: ${OBJS} fabber.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} fabber.o ${LIBS}

mvntool: ${MVNOBJS} mvntool.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BASICOBJS} mvntool.o ${LIBS}



# DO NOT DELETE

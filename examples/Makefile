include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_exp

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -I..
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
else 
  NIFTILIB = -lNewNifti
endif

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat ${NIFTILIB} -lznz -lz -ldl

XFILES = fabber_exp

# Forward models
OBJS =  fwdmodel_exp.o

# For debugging:
#OPTFLAGS = -ggdb

# Build

all:	${XFILES} libfabber_models_exp.a

# Static library
libfabber_models_exp.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specified in this project
fabber_exp: fabber_main.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# DO NOT DELETE

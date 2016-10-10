cmake_minimum_required(VERSION 3.5)

# Basic objects - things that have nothing directly to do with inference
set(BASIC_SRC tools.cc dataset.cc dist_mvn.cc easylog.cc setup.cc fabber_core.cc)

# Core objects - things that implement the framework for inference
set(CORE_SRC noisemodel.cc fwdmodel.cc inference.cc utils.cc fwdmodel_linear.cc
	     fwdmodel_trivial.cc fwdmodel_poly.cc convergence.cc)

# Inference methods
set(INFERENCE_SRC inference_vb.cc inference_nlls.cc inference_spatialvb.cc)

# Noise models
set(NOISE_SRC noisemodel_white.cc noisemodel_ar.cc)

# Tests
set(TEST_SRC test/fabbertest.cc test/test_inference.cc test/test_spatialvb.cc test/test_vb.cc
	     test/test_convergence.cc test/test_commandline.cc test/test_rundata.cc)

# Find our dependencies
find_path(NEWMAT_INCLUDE_DIR newmat.h
          HINTS /home/martinc/dev/fsl/extras/src/newmat REQUIRED)
find_path(UTILS_INCLUDE_DIR tracer_plus.h
          HINTS /home/martinc/dev/fsl/src/utils REQUIRED)

if (DEFINED ENV{ZNZDIR})
  find_library(ZNZ_LIBRARY NAMES znz libznz
               HINTS $ENV{ZNZDIR} REQUIRED )
else ()
  find_library(ZNZ_LIBRARY NAMES znz libznz
               HINTS $ENV{FSLDIR}/lib REQUIRED )
endif()

find_library(UTILS_LIBRARY NAMES utils libutils
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(MISCMATHS_LIBRARY NAMES miscmaths libmiscmaths
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(FSLIO_LIBRARY NAMES fslio libfslio
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(NIFTIIO_LIBRARY NAMES niftiio libniftiio
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(NEWMAT_LIBRARY NAMES newmat libnewmat
             HINTS $ENV{FSLDIR}/extras/lib REQUIRED )
find_library(NEWIMAGE_LIBRARY NAMES newimage libnewimage
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(GTEST_LIBRARY NAMES gtest libgtest REQUIRED )
find_library(PTHREAD_LIBRARY NAMES pthread libpthread REQUIRED )
find_library(Z_LIBRARY NAMES z libz REQUIRED )

include_directories(SYSTEM ${CMAKE_SOURCE_DIR} ${NEWMAT_INCLUDE_DIR} ${UTILS_INCLUDE_DIR}/..)

set(LIBS ${NEWIMAGE_LIBRARY} 
         ${FSLIO_LIBRARY} 
         ${MISCMATHS_LIBRARY} 
         ${NIFTIIO_LIBRARY} 
         ${UTILS_LIBRARY}   
         ${NEWMAT_LIBRARY} ${ZNZ_LIBRARY} ${Z_LIBRARY})

add_definitions(-DUSE_NEWIMAGE -DUSE_NLLS -DVERSION="4.0pre" 
                -DFABBER_SRC_DIR="${CMAKE_SOURCE_DIR}"
                -DFABBER_BUILD_DIR="${CMAKE_BINARY_DIR}")

add_library(fabbercore ${BASIC_SRC} ${CORE_SRC} ${INFERENCE_SRC} ${NOISE_SRC})

add_executable(fabber fabber_main.cc fabber_core.cc)
target_link_libraries(fabber fabbercore ${LIBS})

add_executable(mvntool mvntool.cc )
target_link_libraries(mvntool fabbercore ${LIBS})

add_executable(testfabber ${TEST_SRC})
target_link_libraries(testfabber fabbercore ${LIBS} ${GTEST_LIBRARY} ${PTHREAD_LIBRARY})


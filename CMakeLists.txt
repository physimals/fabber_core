cmake_minimum_required(VERSION 3.5)

# Basic objects - things that have nothing directly to do with inference
set(BASIC_SRC tools.cc dataset.cc dist_mvn.cc easylog.cc setup.cc)

# Core objects - things that implement the framework for inference
set(CORE_SRC noisemodel.cc fwdmodel.cc inference.cc utils.cc fwdmodel_linear.cc
	     fwdmodel_trivial.cc fwdmodel_poly.cc convergence.cc)

# Inference methods
set(INFERENCE_SRC inference_vb.cc inference_nlls.cc inference_spatialvb.cc)

# Noise models
set(NOISE_SRC noisemodel_white.cc noisemodel_ar.cc)

# Main dependencies
find_path(NEWMAT_INCLUDE_DIR newmat.h
          HINTS /home/martinc/dev/fsl/extras/src/newmat REQUIRED)
find_path(UTILS_INCLUDE_DIR tracer_plus.h
          HINTS /home/martinc/dev/fsl/src/utils REQUIRED)

find_library(ZNZ_LIBRARY NAMES znz libznz
             HINTS $ENV{FSLDIR}/lib REQUIRED )
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
find_library(Z_LIBRARY NAMES z libz REQUIRED )

Message("-- Using libznz: ${ZNZ_LIBRARY}")
Message("-- Using libutils: ${UTILS_LIBRARY} ${UTILS_INCLUDE_DIR}")
Message("-- Using miscmaths: ${MISCMATHS_LIBRARY}")
Message("-- Using fslio: ${FSLIO_LIBRARY}")
Message("-- Using niftiio: ${NIFTIIO_LIBRARY}")
Message("-- Using newmat: ${NEWMAT_LIBRARY} ${NEWMAT_INCLUDE_DIR}")
Message("-- Using newimage: ${NEWIMAGE_LIBRARY}")
Message("-- Using libz: ${Z_LIBRARY}")

# Basic build configuration
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

# Main Targets
add_library(fabbercore ${BASIC_SRC} ${CORE_SRC} ${INFERENCE_SRC} ${NOISE_SRC})

add_executable(fabber fabber_main.cc)
target_link_libraries(fabber fabbercore ${LIBS})

add_executable(mvntool mvntool/mvntool.cc )
target_link_libraries(mvntool fabbercore ${LIBS})

INSTALL(TARGETS fabber mvntool fabbercore
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# Conditional Targets

find_library(GTEST_LIBRARY NAMES gtest libgtest OPTIONAL )
find_library(PTHREAD_LIBRARY NAMES pthread libpthread OPTIONAL )

if (GTEST_LIBRARY)
  Message("-- Using gtest: ${GTEST_LIBRARY}")
  set(TEST_SRC test/fabbertest.cc test/test_inference.cc test/test_spatialvb.cc test/test_vb.cc
               test/test_convergence.cc test/test_commandline.cc test/test_rundata.cc)
  add_executable(testfabber ${TEST_SRC})
  target_link_libraries(testfabber fabbercore ${LIBS} ${GTEST_LIBRARY} ${PTHREAD_LIBRARY})
  enable_testing()
  add_test(COMMAND testfabber)
else(GTEST_LIBRARY)
  Message("-- Gtest not found - will not build unit tests")
endif(GTEST_LIBRARY)

find_package(Doxygen)

if(DOXYGEN_FOUND)
 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
 add_custom_target(doc
 ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
 COMMENT "Generating API documentation with Doxygen" VERBATIM
)
else(DOXYGEN_FOUND)
  Message("-- Doxygen not found - will not build API documentation")
endif(DOXYGEN_FOUND)

find_package(PythonLibs)
find_program(PYSIDE_UIC pyside-uic)

if (PYSIDE_UIC) 
  Message("-- Pyside found: ${PYSIDE_UIC}")
  set(PY_SRC pyfab/pyfab/ui.py pyfab/pyfab/views.py pyfab/pyfab/model.py pyfab/pyfab/plotviews.py pyfab/pyfab/__init__.py)
  configure_file(pyfab/pyfabber.py.in pyfab/pyfabber.py)
  add_custom_target(pyfab_ui ALL
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_mainwin.ui >ui_qtd.py
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_optdialog.ui >>ui_qtd.py
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_logview.ui >>ui_qtd.py
                   )

  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pyfab/pyfabber.py DESTINATION bin)
  INSTALL(FILES ${PY_SRC} DESTINATION lib/python/pyfab)
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/ui_qtd.py DESTINATION lib/python/pyfab)

else(PYSIDE_UIC)
  Message("-- Pyside not found - will not build GUI")
endif(PYSIDE_UIC)


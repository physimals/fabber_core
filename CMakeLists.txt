cmake_minimum_required(VERSION 2.7)

# Basic objects - things that have nothing directly to do with inference
set(BASIC_SRC fabber_core.cc tools.cc dataset.cc dist_mvn.cc easylog.cc setup.cc)

# Core objects - things that implement the framework for inference
set(CORE_SRC noisemodel.cc fwdmodel.cc inference.cc utils.cc fwdmodel_linear.cc
	     fwdmodel_trivial.cc fwdmodel_poly.cc convergence.cc)

# Inference methods
set(INFERENCE_SRC inference_vb.cc inference_nlls.cc inference_spatialvb.cc)

# Noise models
set(NOISE_SRC noisemodel_white.cc noisemodel_ar.cc)

# Main dependencies

if (FSL_BUILD)
  set(LIBS newimage fslio miscmaths niftiio utils newmat znz ${ZLIB_LIBRARY})
  include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR})
else (FSL_BUILD)

  find_path(FSL_INCLUDE_DIR utils/tracer_plus.h
            HINTS $ENV{FSLDIR}/include REQUIRED)
  find_path(FSL_EXTRAS_INCLUDE_DIR cprob/cprob.h
            HINTS $ENV{FSLDIR}/extras/include REQUIRED)
  find_path(NEWMAT_INCLUDE_DIR newmat.h
            HINTS $ENV{FSLDIR}/extras/include/newmat REQUIRED)
  
  find_path(MCFLIRT_INCLUDE_DIR /mcflirt/rigidreglib.h
            HINTS $ENV{FSLDIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/.. OPTIONAL)
  if(NOT MCFLIRT_INCLUDE_DIR)
     Message("-- MCFLIRT not found - will not enable motion corrrection")
  else(NOT MCFLIRT_INCLUDE_DIR)
     Message("-- Found MCFLIRT ${MCFLIRT_INCLUDE_DIR}")
     find_library(WARPFNS_LIBRARY NAMES warpfns libwarpfns
               HINTS $ENV{FSLDIR}/lib OPTIONAL )
     if(WARPFNS_LIBRARY)
       Message("-- Found warpfns ${WARPFNS_LIBRARY}")
       include_directories(${MCFLIRT_INCLUDE_DIR})
       add_definitions(-D__FABBER_MOTION -DEXPOSE_TREACHEROUS)
       set(CORE_SRC ${CORE_SRC} Update_deformation.cc)
     endif(WARPFNS_LIBRARY)
  endif(NOT MCFLIRT_INCLUDE_DIR)

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
  find_library(Z_LIBRARY NAMES z libz zlib zlibd 
               HINTS $ENV{FSLDIR}/lib REQUIRED )

  find_package(Boost)
  if(NOT Boost_FOUND)
      find_path(Boost_INCLUDE_DIRS boost/shared_ptr.hpp
                HINTS $ENV{FSLDIR}/extras/include/boost REQUIRED)
  endif(NOT Boost_FOUND)
  
  find_package(ZLIB)
  if (NOT ZLIB_FOUND)
    find_library(ZLIB_LIBRARY NAMES z libz zlib zlibd 
                 HINTS $ENV{FSLDIR}/lib REQUIRED )
  endif(NOT ZLIB_FOUND)
   
Message("-- FSL headers in ${FSL_INCLUDE_DIR}")
Message("-- Using libznz: ${ZNZ_LIBRARY}")
Message("-- Using libutils: ${UTILS_LIBRARY} ${UTILS_INCLUDE_DIR}")
Message("-- Using miscmaths: ${MISCMATHS_LIBRARY}")
Message("-- Using fslio: ${FSLIO_LIBRARY}")
Message("-- Using niftiio: ${NIFTIIO_LIBRARY}")
Message("-- Using newmat: ${NEWMAT_LIBRARY} ${NEWMAT_INCLUDE_DIR}")
Message("-- Using newimage: ${NEWIMAGE_LIBRARY}")
Message("-- Using zlib: ${ZLIB_LIBRARY}")

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR} ${NEWMAT_INCLUDE_DIR} ${FSL_INCLUDE_DIR} ${FSL_EXTRAS_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

set(LIBS ${NEWIMAGE_LIBRARY} ${WARPFNS_LIBRARY}
         ${FSLIO_LIBRARY} 
         ${MISCMATHS_LIBRARY} 
         ${NIFTIIO_LIBRARY} 
         ${UTILS_LIBRARY}   
         ${NEWMAT_LIBRARY} ${ZNZ_LIBRARY} ${ZLIB_LIBRARY})
endif(FSL_BUILD)

# Basic build configuration

# the commit's SHA1, and whether the building workspace was dirty or not
execute_process(COMMAND
  git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the date of the commit
execute_process(COMMAND
  git log -1 --format=%ad --date=local
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_DATE
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)


if (NOT GIT_SHA1)
  set(GIT_SHA1 unknown)
endif(NOT GIT_SHA1)

if (NOT GIT_DATE)
  set(GIT_DATE unknown)
endif(NOT GIT_DATE)

add_definitions(-DGIT_SHA1="${GIT_SHA1}" -DGIT_DATE="${GIT_DATE}"
                -DFABBER_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
                -DFABBER_BUILD_DIR="${CMAKE_CURRENT_BINARY_DIR}")

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

#INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" )
INSTALL_FILES(/include/fabber_core .*\\.h$)

# Conditional Targets

find_library(GTEST_LIBRARY NAMES gtest libgtest OPTIONAL )
find_library(PTHREAD_LIBRARY NAMES pthread libpthread OPTIONAL )
find_path(GTEST_INCLUDE_DIR gtest/gtest.h OPTIONAL )

if (GTEST_LIBRARY)
  Message("-- Using gtest: ${GTEST_LIBRARY} ${GTEST_INCLUDE_DIR}")  
  if (PTHREAD_LIBRARY)
    set(PTH_LIB ${PTHREAD_LIBRARY})
  endif(PTHREAD_LIBRARY)
  include_directories(${GTEST_INCLUDE_DIR})
  
  set(TEST_SRC test/fabbertest.cc test/test_inference.cc test/test_spatialvb.cc test/test_vb.cc
               test/test_convergence.cc test/test_commandline.cc test/test_rundata.cc test/test_priors.cc)
  add_executable(testfabber ${TEST_SRC})
  target_link_libraries(testfabber fabbercore ${LIBS} ${GTEST_LIBRARY} ${PTH_LIB})
  enable_testing()
  add_test(COMMAND testfabber)
else(GTEST_LIBRARY)
  Message("-- Gtest NOT found - will not build unit tests")
endif(GTEST_LIBRARY)

find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/html DESTINATION share/doc/fabber)

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
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_chooserun.ui >>ui_qtd.py
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_matrix.ui >>ui_qtd.py
                   )

  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pyfab/pyfabber.py DESTINATION bin)
  INSTALL(DIRECTORY pyfab/pyfab DESTINATION lib/python/)
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/ui_qtd.py DESTINATION lib/python/pyfab)

else(PYSIDE_UIC)
  Message("-- Pyside not found - will not build GUI")
endif(PYSIDE_UIC)


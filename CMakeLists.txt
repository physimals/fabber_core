cmake_minimum_required(VERSION 2.7)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(APPLE)
  # Uncomment to force use of libstdc++ on Mac as
  # FSL uses this by default. Not required if building FSL from
  # Cmake sources
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libstdc++")
  # Get rid of CMake warning
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long")
endif(UNIX)

# Basic objects - things that have nothing directly to do with inference
set(BASIC_SRC tools.cc rundata.cc dist_mvn.cc easylog.cc setup.cc fabber_capi.cc rundata_array.cc dist_gamma.cc version.cc)

# Core objects - things that implement the framework for inference
set(CORE_SRC noisemodel.cc fwdmodel.cc inference.cc factories.cc fwdmodel_linear.cc
	           fwdmodel_poly.cc convergence.cc motioncorr.cc)

# Inference methods
set(INFERENCE_SRC inference_vb.cc inference_nlls.cc inference_spatialvb.cc)

# Noise models
set(NOISE_SRC noisemodel_white.cc noisemodel_ar.cc)

# Find FSL dependencies

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

  #find_path(MCFLIRT_INCLUDE_DIR /mcflirt/rigidreglib.h
  #          HINTS $ENV{FSLDIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/.. OPTIONAL)
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
     endif(WARPFNS_LIBRARY)
  endif(NOT MCFLIRT_INCLUDE_DIR)

  find_library(ZNZ_LIBRARY NAMES znz libznz
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

  find_package(Boost)
  if(NOT Boost_FOUND)
      find_path(Boost_INCLUDE_DIRS boost/shared_ptr.hpp
                HINTS $ENV{FSLDIR}/extras/include/boost REQUIRED)
  endif(NOT Boost_FOUND)

  if(NOT WIN32)
    find_package(ZLIB)
  endif(NOT WIN32)
  if (NOT ZLIB_FOUND)
      find_library(ZLIB_LIBRARY NAMES z libz zlib zlibd
                   HINTS $ENV{FSLDIR}/lib REQUIRED NO_DEFAULT_PATH)
  endif(NOT ZLIB_FOUND)

  Message("-- FSL headers in ${FSL_INCLUDE_DIR} ${FSL_EXTRAS_INCLUDE_DIR} ${NEWMAT_INCLUDE_DIR}")
  Message("-- Using newimage: ${NEWIMAGE_LIBRARY}")
  Message("-- Using fslio: ${FSLIO_LIBRARY}")
  Message("-- Using miscmaths: ${MISCMATHS_LIBRARY}")
  Message("-- Using niftiio: ${NIFTIIO_LIBRARY}")
  Message("-- Using libznz: ${ZNZ_LIBRARY}")
  Message("-- Using zlib: ${ZLIB_LIBRARY}")
  Message("-- Using newmat: ${NEWMAT_LIBRARY}")
  Message("-- Using Boost: ${Boost_INCLUDE_DIRS}")

  include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR} ${NEWMAT_INCLUDE_DIR} ${FSL_INCLUDE_DIR} ${FSL_EXTRAS_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

  set(LIBS ${NEWIMAGE_LIBRARY} ${WARPFNS_LIBRARY}
           ${FSLIO_LIBRARY}
           ${MISCMATHS_LIBRARY}
           ${NIFTIIO_LIBRARY}
           ${NEWMAT_LIBRARY} ${ZNZ_LIBRARY} ${ZLIB_LIBRARY})
endif(FSL_BUILD)

if (UNIX)
  set(LIBS ${LIBS} dl)
endif(UNIX)

# Versioning information

# the commit's SHA1, and whether the building workspace was dirty or not
execute_process(COMMAND
  git describe --dirty
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

add_library(fabbercore STATIC ${BASIC_SRC} ${CORE_SRC} ${INFERENCE_SRC} ${NOISE_SRC})
add_library(fabbercore_shared SHARED ${BASIC_SRC} ${CORE_SRC} ${INFERENCE_SRC} ${NOISE_SRC})
target_link_libraries(fabbercore_shared ${LIBS})
add_library(fabberexec rundata_newimage.cc fabber_core.cc )

add_executable(fabber fabber_main.cc)
target_link_libraries(fabber fabberexec fabbercore ${LIBS})

add_executable(mvntool mvn_tool/mvntool.cc )
target_link_libraries(mvntool fabberexec fabbercore ${LIBS})

add_executable(niftidiff test/niftidiff.cc )
target_link_libraries(niftidiff ${LIBS})

INSTALL(TARGETS fabber mvntool fabbercore fabbercore_shared fabberexec
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

#INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" )
INSTALL_FILES(/include/fabber_core .*\\.h$)
INSTALL(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/fabber_var DESTINATION bin)

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
  target_link_libraries(testfabber fabbercore fabberexec ${LIBS} ${GTEST_LIBRARY} ${PTH_LIB})
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
  configure_file(pyfab/pyfabber.py.in pyfab/pyfabber.py)
  add_custom_target(pyfab_ui ALL
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_mainwin.ui >ui_qtd.py
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_optdialog.ui >>ui_qtd.py
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_logview.ui >>ui_qtd.py
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_chooserun.ui >>ui_qtd.py
                   COMMAND ${PYSIDE_UIC} ${CMAKE_CURRENT_SOURCE_DIR}/pyfab/pyfab_matrix.ui >>ui_qtd.py
                   )

  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pyfab/pyfabber.py DESTINATION bin)
  INSTALL(DIRECTORY pyfab/fabber DESTINATION lib/python/)
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/ui_qtd.py DESTINATION lib/python/fabber)

else(PYSIDE_UIC)
  Message("-- Pyside not found - will not build GUI")
endif(PYSIDE_UIC)

# Packaging
include(CPackConfig.cmake)

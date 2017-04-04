cmake_minimum_required(VERSION 2.7)

# For single model you should only need to change the line below, if your files keep to the naming convention!
set(MODELS_NAME sine)

# For multiple models you may need to add your model source code files below
set(MODELS_SRC fwdmodel_${MODELS_NAME}.cc)

add_library(fabber_models_${MODELS_NAME} SHARED ${MODELS_SRC} ${MODELS_NAME}_models.cc)
add_executable(fabber_${MODELS_NAME} ${MODELS_SRC} fabber_main.cc)

if(APPLE)
  # Uncomment to force use of libstdc++ on Mac (if FSL build uses this and you cannot link otherwise)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libstdc++")
  # Get rid of CMake warning
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

find_path(FABBER_INCLUDE_DIR fabber_core/fwdmodel.h
            HINTS $ENV{FSLDIR}/include .. REQUIRED)
find_path(NEWMAT_INCLUDE_DIR newmat.h
            HINTS $ENV{FSLDIR}/extras/include/newmat REQUIRED)
find_path(FSL_INCLUDE_DIR newimage/newimage.h
            HINTS $ENV{FSLDIR}/include REQUIRED)
find_path(FSL_EXT_INCLUDE_DIR libprob/libprob.h
            HINTS $ENV{FSLDIR}/extras/include REQUIRED)
find_path(BOOST_INCLUDE_DIR boost/shared_ptr.hpp
          HINTS $ENV{FSLDIR}/extras/include/boost $ENV{FSLDIR}/extras/include/ REQUIRED)

Message("-- FSL headers in ${FSL_INCLUDE_DIR} ${NEWMAT_INCLUDE_DIR} ${FSL_EXT_INCLUDE_DIR} ${BOOST_INCLUDE_DIR}")
Message("-- Fabber headers in ${FABBER_INCLUDE_DIR}/fabber_core")

find_library(FABBER_LIBRARY NAMES fabbercore libfabbercore
             HINTS $ENV{FSLDIR}/lib $ENV{FABBERDIR}/lib ../fabber_core/build_release ../fabber_core/build_debug ../fabber_core/build_relwithdebinfo REQUIRED )

find_library(FABBER_EXEC_LIBRARY NAMES fabberexec libfabberexec
             HINTS $ENV{FSLDIR}/lib $ENV{FABBERDIR}/lib ../fabber_core/build_release ../fabber_core/build_debug ../fabber_core/build_relwithdebinfo REQUIRED )

find_library(NEWMAT_LIBRARY NAMES newmat libnewmat
             HINTS $ENV{FSLDIR}/extras/lib REQUIRED )
find_library(MISCMATHS_LIBRARY NAMES miscmaths libmiscmaths
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(NEWIMAGE_LIBRARY NAMES newimage libnewimage
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(FSLIO_LIBRARY NAMES fslio libfslio
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(NIFTIIO_LIBRARY NAMES niftiio libniftiio
             HINTS $ENV{FSLDIR}/lib REQUIRED )
find_library(ZNZ_LIBRARY NAMES znz libznz
             HINTS $ENV{FSLDIR}/lib REQUIRED )

if (WINDOWS)
  find_library(Z_LIBRARY NAMES z libz zlib zlibd 
                 HINTS $ENV{FSLDIR}/lib REQUIRED NO_DEFAULT_PATH)
else(WINDOWS)
  find_library(Z_LIBRARY NAMES z libz zlib zlibd 
                 HINTS $ENV{FSLDIR}/lib REQUIRED)
endif(WINDOWS)

set(LIBS ${FABBER_LIBRARY} ${FABBER_EXEC_LIBRARY}  ${NEWIMAGE_LIBRARY} ${MISCMATHS_LIBRARY} ${FSLIO_LIBRARY} ${NIFTIIO_LIBRARY}  ${ZNZ_LIBRARY} ${Z_LIBRARY} ${NEWMAT_LIBRARY})

if (UNIX)
  set(LIBS ${LIBS} dl)
endif(UNIX)

Message("-- Using Fabber libraries: ${FABBER_LIBRARY} ${FABBER_EXEC_LIBRARY}")
Message("-- Using libznz: ${ZNZ_LIBRARY}")
Message("-- Using miscmaths: ${MISCMATHS_LIBRARY}")
Message("-- Using fslio: ${FSLIO_LIBRARY}")
Message("-- Using newimage: ${NEWIMAGE_LIBRARY}")
Message("-- Using niftiio: ${NIFTIIO_LIBRARY}")
Message("-- Using newmat: ${NEWMAT_LIBRARY} ${NEWMAT_INCLUDE_DIR}")
Message("-- Using newimage: ${NEWIMAGE_LIBRARY}")
Message("-- Using zlib: ${Z_LIBRARY}")

target_link_libraries(fabber_models_${MODELS_NAME} ${LIBS})
target_link_libraries(fabber_${MODELS_NAME} ${LIBS})

include_directories(SYSTEM ${FABBER_INCLUDE_DIR} ${NEWMAT_INCLUDE_DIR} ${FSL_INCLUDE_DIR} ${FSL_EXT_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})

# Extract GIT revision for versioning

execute_process(COMMAND
  git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

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

add_definitions(-DGIT_SHA1="${GIT_SHA1}" -DGIT_DATE="${GIT_DATE}")

INSTALL(TARGETS fabber_${MODELS_NAME} fabber_models_${MODELS_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
       
# Packaging
#include(CPackConfig.cmake)


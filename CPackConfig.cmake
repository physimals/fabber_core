# Generic
SET(CPACK_PACKAGE_NAME "fabber")
SET(CPACK_PACKAGE_VENDOR "University of Oxford, IBME")
SET(CPACK_PACKAGE_VERSION "3.9.0")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_PACKAGE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

# Debian package
set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/fabber")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "martin.craig@eng.ox.ac.uk")

# Windows installer
if(BUILD_WIN32_NSIS)
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_PACKAGING_INSTALL_PREFIX "c:\Program Files\fabber")
endif(BUILD_WIN32_NSIS)

# OSX app bundle
if(BUILD_OSX_BUNDLE)
  set(CPACK_BUNDLE_NAME "fabber")
endif(BUILD_OSX_BUNDLE)

# Have to do this last for the above to have any effect
include(CPack)


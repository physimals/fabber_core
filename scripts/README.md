# Build scripts

Fabber now uses CMake with a view to building on Win32

A few simple shell scripts are provided for Unix-style environments:

  build.sh <type>

Rebuild from scratch. `<type>` is the build type - `Debug` or `Release`. The
default is `Debug`. The built files will be placed in the directory
`build_<type>`

  install.sh <type> <prefix>

Build and install. Build type is as above. The installation
prefix defaults to `$FSLDIR` if set, `$HOME` if not. To install globally 
you might do:

  install.sh Release /usr/local/

  package.sh

This is to build a binary package, and is generally only used when preparing
a new release.

The scripts are very simple and you may prefer to run cmake directly instead.





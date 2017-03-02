# Build scripts

Fabber now uses CMake with a view to building on Win32. A few simple shell
scripts are provided to help. The scripts are very simple and you may prefer to
run cmake directly instead.

## Unix entironments

    build.sh <type>

Rebuild from scratch. `<type>` is the build type - `debug` or `release`. The
default is `debug`. The built files will be placed in the directory
`build_<type>`

    install.sh <type> <prefix>

Build and install. Build type is as above. The installation
prefix defaults to `$FSLDIR` if set, `$HOME` if not. To install globally
you might do:

    install.sh release /usr/local/

This will install under ``/usr/local/bin`, etc.

    package.sh

This is to build binary packages, and is generally only used when preparing
a new release.

## Windows

Windows build scripts should be run from a Visual Studio command prompt.

    build.bat <arch> <type>

Build the executables from scratch. `arch` is the build architecture, `x86` or
`x64`. `type` is `debug` or `release`

Note that the Windows ABI is sensitive to debug and release builds so if you
want to link to other software you must make sure you are building in the same
mode, and are using the same compiler version.

    install.bat <arch> <type>

Build and install the code. By default installs into `%HOME%/fsl_<arch>_<type>`,
if this is unsuitable, edit the script directly.

    package.bat <arch> <type>

Build packages. Used when preparing a new release.

    build_all.bat
    install_all.bat
    package_all.bat

For building all architecture variants. Used for preparing releases.

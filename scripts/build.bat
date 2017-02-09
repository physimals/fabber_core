echo "You must run this from a Visual Studio developer commmand line"

set ARCH=%1
set TYPE=%2
set ORIGFSLDIR=%FSLDIR%
set ORIGDIR=%cd%
set ORIGPATH=%PATH%

set FSLDIR=%HOMEDRIVE%%HOMEPATH%\fsl_%ARCH%_%TYPE%
set BUILDDIR= %0\..\..\build_%ARCH%_%TYPE%

call "%VCINSTALLDIR%\vcvarsall" %ARCH%

rd /s /q %BUILDDIR%
mkdir %BUILDDIR%
cd %BUILDDIR%
cmake .. -DCMAKE_BUILD_TYPE=%TYPE% -G "NMake Makefiles"
nmake
cd ..

set PATH=%ORIGPATH%
set FSLDIR=%ORIGFSLDIR%
cd %ORIGDIR%

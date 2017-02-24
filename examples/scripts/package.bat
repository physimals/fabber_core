
set ARCH=%1
set TYPE=%2

set FSLDIR=%HOMEDRIVE%%HOMEPATH%/fsl_%ARCH%_%TYPE%
set BUILDDIR=%0\..\..\build_%ARCH%_%TYPE%
set ORIGPATH=%PATH%
set ORIGDIR=%cd%

call "%VCINSTALLDIR%\vcvarsall" %ARCH%

rd /s /q %BUILDDIR%
mkdir %BUILDDIR%
cd %BUILDDIR%
cmake .. -DCMAKE_BUILD_TYPE=%TYPE% -G "NMake Makefiles"
nmake package

set PATH=%ORIGPATH%
cd %ORIGDIR%

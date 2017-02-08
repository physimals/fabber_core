
set ARCH=%1
set TYPE=%2

set FSLDIR=c:\Users\ctsu0221\fsl_vs2013_%ARCH%_%TYPE%
set BUILDDIR=build_%ARCH%_%TYPE%
set ORIGPATH=%PATH%
call "%VCINSTALLDIR%\vcvarsall" %ARCH%

rd /s /q %BUILDDIR%
mkdir %BUILDDIR%
cd %BUILDDIR%
cmake .. -DCMAKE_BUILD_TYPE=%TYPE% -G "NMake Makefiles"
nmake
nmake package
cd ..

set PATH=%ORIGPATH%

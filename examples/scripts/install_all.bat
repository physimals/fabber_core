echo "You must run this from a Visual Studio developer commmand line"

call %0\..\install.bat x64 debug

call %0\..\install.bat x64 release

call %0\..\install.bat x86 debug

call %0\..\install.bat x86 release

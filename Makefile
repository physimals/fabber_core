.PHONY: all clean install doc check mvntool
default:    all
all:    
	scons
clean:
	scons -c
install:
	scons install
doc:
	scons doc
check:
	scons test
mvntool:
	scons mvntool

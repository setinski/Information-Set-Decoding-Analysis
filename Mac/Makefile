CC=clang++
IPATHS=-I<TOOLS_DIR>/platform/<PLATFORM>/h -I<BOOST_DIR>
LPATHS=-L<TOOLS_DIR>/platform/<PLATFORM>/bin
OBJECTS = main.o isd.o space.o misc.o

InformationSetDecoding: $(OBJECTS)
	$(CC) -Wl,-headerpad_max_install_names -std=c++11 -stdlib=libc++ -g $(IPATHS) $(LPATHS) -o InformationSetDecoding $(OBJECTS) -lfusion64 -lmosek64
	install_name_tool -change libfusion64.9.2.dylib <TOOLS_DIR>/platform/<PLATFORM>/bin/libfusion64.9.2.dylib InformationSetDecoding || rm -f InformationSetDecoding
	install_name_tool -change libmosek64.9.2.dylib <TOOLS_DIR>/platform/<PLATFORM>/bin/libmosek64.9.2.dylib InformationSetDecoding || rm -f InformationSetDecoding
	rm -f $(OBJECTS)
	mkdir -p ../bin
	mv InformationSetDecoding ../bin/InformationSetDecoding

../src/main.o: ../src/isd.h ../src/space.h ../src/misc.h
	$(CC) -std=c++11 -stdlib=libc++ -g $(IPATHS) -c ../src/main.cpp

../src/isd.o: ../src/isd.h ../src/space.h ../src/misc.h
	$(CC) -std=c++11 -stdlib=libc++ -g $(IPATHS) -c ../src/isd.cpp

../src/space.o: ../src/space.h ../src/misc.h
	$(CC) -std=c++11 -stdlib=libc++ -g $(IPATHS) -c ../src/space.cpp

../src/misc.o: ../src/misc.h
	$(CC) -std=c++11 -stdlib=libc++ -c ../src/misc.cpp

.PHONY: clean

clean:
	rm -rf ../bin

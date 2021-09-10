CC=clang++
IPATHS=-I<TOOLS_PATH>/platform/<PLATFORM>/h -I<BOOST_PATH>
LPATHS=-L<TOOLS_PATH>/platform/<PLATFORM>/bin
OBJECTS = main.o isd.o space.o misc.o

InformationSetDecoding: $(OBJECTS)
	$(CC) -Wl,-headerpad_max_install_names -std=c++11 -stdlib=libc++ -g $(IPATHS) $(LPATHS) -o InformationSetDecoding $(OBJECTS) -lfusion64 -lmosek64
	install_name_tool -change libfusion64.9.2.dylib <TOOLS_PATH>/platform/<PLATFORM>/bin/libfusion64.9.2.dylib InformationSetDecoding || rm -f InformationSetDecoding
	install_name_tool -change libmosek64.9.2.dylib <TOOLS_PATH>/platform/<PLATFORM>/bin/libmosek64.9.2.dylib InformationSetDecoding || rm -f InformationSetDecoding
	rm -f $(OBJECTS)

main.o: isd.h space.h misc.h
	$(CC) -std=c++11 -stdlib=libc++ -g $(IPATHS) -c main.cpp

isd.o: isd.h space.h misc.h
	$(CC) -std=c++11 -stdlib=libc++ -g $(IPATHS) -c isd.cpp

space.o: space.h misc.h
	$(CC) -std=c++11 -stdlib=libc++ -g $(IPATHS) -c space.cpp

misc.o: misc.h
	$(CC) -std=c++11 -stdlib=libc++ -c misc.cpp

.PHONY: clean

clean:
	rm -f InformationSetDecoding *.txt


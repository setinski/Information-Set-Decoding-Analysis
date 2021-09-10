CC=clang++
IPATHS=-I<TOOLS_PATH>/platform/<PLATFORM>/h -I<BOOST_PATH>
LPATHS=-L$<TOOLS_PATH>/platform/<PLATFORM>/bin -Wl,-rpath-link,<TOOLS_PATH>/platform/<PLATFORM>/bin '-Wl,-rpath=<TOOLS_PATH>/platform/<PLATFORM>/bin'
OBJECTS = main.o isd.o space.o misc.o

InformationSetDecoding: $(OBJECTS)
	$(CC) -pthread -std=c++11 -g $(IPATHS) $(LPATHS) -o InformationSetDecoding $(OBJECTS) -lfusion64 -lmosek64
	rm -f $(OBJECTS)

main.o: isd.h space.h misc.h
	$(CC) -pthread -std=c++11 -g $(IPATHS) -c main.cpp

isd.o: isd.h space.h misc.h
	$(CC) -pthread -std=c++11 -g $(IPATHS) -c isd.cpp

space.o: space.h misc.h
	$(CC) -pthread -std=c++11 -g $(IPATHS) -c space.cpp

misc.o: misc.h
	$(CC) -pthread -std=c++11 -c misc.cpp


.PHONY: clean

clean:
	rm -f InformationSetDecoding *.txt

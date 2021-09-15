CC=clang++
IPATHS=-I<TOOLS_DIR>/platform/<PLATFORM>/h -I<BOOST_DIR>
LPATHS=-L$<TOOLS_DIR>/platform/<PLATFORM>/bin -Wl,-rpath-link,<TOOLS_DIR>/platform/<PLATFORM>/bin '-Wl,-rpath=<TOOLS_DIR>/platform/<PLATFORM>/bin'
OBJECTS = main.o isd.o space.o misc.o

InformationSetDecoding: $(OBJECTS)
	$(CC) -pthread -std=c++11 -g $(IPATHS) $(LPATHS) -o InformationSetDecoding $(OBJECTS) -lfusion64 -lmosek64
	rm -f $(OBJECTS)
	mkdir -p ../bin
	mv InformationSetDecoding ../bin/InformationSetDecoding


main.o: ../src/isd.h ../src/space.h ../src/misc.h
	$(CC) -pthread -std=c++11 -g $(IPATHS) -c ../src/main.cpp

isd.o: ../src/isd.h ../src/space.h ../src/misc.h
	$(CC) -pthread -std=c++11 -g $(IPATHS) -c ../src/isd.cpp

space.o: ../src/space.h ../src/misc.h
	$(CC) -pthread -std=c++11 -g $(IPATHS) -c ../src/space.cpp

misc.o: ../src/misc.h
	$(CC) -pthread -std=c++11 -c ../src/misc.cpp


.PHONY: clean

clean:
	rm -rf ../bin


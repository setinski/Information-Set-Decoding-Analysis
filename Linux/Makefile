CC=clang++
IPATHS=-I<TOOLS_DIR>/platform/<PLATFORM>/h -I<BOOST_DIR>
LPATHS=-L$<TOOLS_DIR>/platform/<PLATFORM>/bin -Wl,-rpath-link,<TOOLS_DIR>/platform/<PLATFORM>/bin '-Wl,-rpath=<TOOLS_DIR>/platform/<PLATFORM>/bin'
OBJECTS = main.o isd.o space.o misc.o

InformationSetDecoding: $(OBJECTS)
	$(CC) -pthread -std=c++11 -g $(IPATHS) $(LPATHS) -o InformationSetDecoding $(OBJECTS) -lfusion64 -lmosek64
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


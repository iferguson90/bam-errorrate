BINARY=bam-errorrate

$(BINARY): bam-errorrate.cpp
	g++ -Wall -std=c++0x -Ibamtools/include -lz -I/opt/boost-tgi/1.51/include bam-errorrate.cpp -o $(BINARY) bamtools/lib/libbamtools.a bamtools/lib/libbamtools-utils.a bamtools/lib/libjsoncpp.a

clean:
	rm -f $(BINARY)

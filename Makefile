bam-errorrate: bam-errorrate.cpp
	g++ -Wall -std=c++0x -Ibamtools/include -lz -I/opt/boost-tgi/1.51/include bam-errorrate.cpp -o bam-errorrate bamtools/lib/libbamtools.a bamtools/lib/libbamtools-utils.a bamtools/lib/libjsoncpp.a

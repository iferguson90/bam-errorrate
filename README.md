```
sudo apt-get install libboost1.51-tgi-dev
g++ -Wall -std=c++0x -Ibamtools/include -lz -I/opt/boost-tgi/1.51/include bam-errorrate.cpp -o bam-errorrate bamtools/lib/libbamtools.a bamtools/lib/libbamtools-utils.a bamtools/lib/libjsoncpp.a
```

```
clang++ -Wall -std=c++11 -isystem bamtools/include -Lbamtools/lib -lbamtools -lboost_regex bam-errorrate.cpp -o bam-errorrate
```

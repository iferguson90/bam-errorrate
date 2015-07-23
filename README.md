Build on Linux (on MGI workstations):
```
sudo apt-get install libboost1.51-tgi-dev cmake build-essential git-core
git submodule update --init --recursive
mkdir -p bamtools/build
pushd bamtools/build
cmake ..
make -j8 
popd
g++ -Wall -std=c++0x -Ibamtools/include -lz -I/opt/boost-tgi/1.51/include bam-errorrate.cpp -o bam-errorrate bamtools/lib/libbamtools.a bamtools/lib/libbamtools-utils.a bamtools/lib/libjsoncpp.a
```

Adjust your dependencies and includes for other distros accordingly.

This ought to work on OS X:
```
clang++ -Wall -std=c++11 -isystem bamtools/include -Lbamtools/lib -lbamtools -lboost_regex bam-errorrate.cpp -o bam-errorrate
```

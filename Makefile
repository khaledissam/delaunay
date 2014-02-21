test: test.o delaunay.o
	g++ -o test test.o delaunay.o

test.o: test.cpp
	g++ -c test.cpp

delaunay.o: delaunay.cpp
	g++ -c delaunay.cpp
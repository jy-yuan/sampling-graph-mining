all: main

main: main.o Graph.o IEStop.o
	mpicxx -pthread main.o Graph.o IEStop.o -o main

main.o:
	mpicxx -pthread -O2 -c main.cpp

Graph.o: Graph.cpp Graph.hpp
	g++ -c -O3 Graph.cpp

IEStop.o: IEStop.cpp IEStop.hpp
	g++ -c -O3 IEStop.cpp

clean:
	rm -f *.o main

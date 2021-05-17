all: main

main: main.o graph.o iestop.o
	mpicxx -pthread main.o graph.o iestop.o -o main

main.o:
	mpicxx -pthread -c main.cpp

graph.o: Graph.cpp Graph.hpp
	g++ -c Graph.cpp

iestop.o: IEStop.cpp IEStop.hpp
	g++ -c IEStop.cpp

clean:
	rm -f *.o main
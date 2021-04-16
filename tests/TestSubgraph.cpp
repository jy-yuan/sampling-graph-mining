#include <iostream>

#include "ZGraphInstance.hpp"

using namespace std;

int main(int argc, char ** argv) {
	Edge<Empty> edge;

	SubGraph<Empty> sub_graph_0;
	edge.src = 0; edge.dst = 1; sub_graph_0.add_edge(edge);
	edge.src = 0; edge.dst = 2; sub_graph_0.add_edge(edge);
	edge.src = 1; edge.dst = 2; sub_graph_0.add_edge(edge);
	cout << "subgraph_0: " << sub_graph_0 << endl;

	SubGraph<Empty> sub_graph_1;
	edge.src = 0; edge.dst = 1; sub_graph_1.add_edge(edge);
	edge.src = 0; edge.dst = 2; sub_graph_1.add_edge(edge);
	cout << "subgraph_1: " << sub_graph_1 << endl; 

	cout << "subgraph_0 - subgraph_1 = " << sub_graph_0 - sub_graph_1 << endl;
	return 0;
}

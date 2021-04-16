#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <string>

#include "Type.hpp"

bool file_exists(std::string filename);
long file_size(std::string filename);
void load_graph_meta_data(std::string meta_path, VertexId &num_vertices);
unsigned int get_globally_consistent_seed();

#endif

#ifndef TYPE_HPP
#define TYPE_HPP

#include <stdint.h>

typedef uint32_t VertexId;
typedef uint64_t EdgeId;

struct Empty {
};

template<typename EdgeData> 
struct EdgeUnit {
	VertexId src;
	VertexId dst;
	EdgeData edge_data;
} __attribute__((packed));

template<>
struct EdgeUnit <Empty> {
	VertexId src;
	union {
		VertexId dst;
		Empty edge_data;
	};
} __attribute__((packed));

template<typename EdgeData>
struct AdjUnit {
	EdgeId edge_id;
	VertexId src;
	VertexId dst;
	EdgeData edge_data;
} __attribute__((packed));

template<>
struct AdjUnit<Empty> {
	EdgeId edge_id;
	VertexId src;
	union {
		VertexId dst;
		Empty edge_data;
	};
} __attribute__((packed));

template<typename EdgeData> 
struct VertexAdjList {
	AdjUnit<EdgeData> * begin;
	AdjUnit<EdgeData> * end;
	VertexAdjList(): begin(nullptr), end(nullptr) {
	}
	VertexAdjList(AdjUnit<EdgeData> * begin_, AdjUnit<EdgeData> * end_): begin(begin_), end(end_) {
	}
};

#endif

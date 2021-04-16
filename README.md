# ZGraph 

## Introduction 
ZGraph is an open-source approximate graph mining implementation of the ASAP system mentioned in this paper[1].

## Todolist
[ ] Improve caching performance and implement 3-motif, 4-motif applications
[ ] Implement the sampling generator
[ ] Add complex predicate matching support
[ ] conduct experiments:
	[ ] comparing Arabesque: 6 x 96 GB machines: wiki-vote, youtube, patents, live-journal: triangle, 4-clique, 3-motif
	[ ] scalability:
		[ ] single node scalability: shared memory (threads)
		[ ] distributed scalability
	[ ] experiment on very large graph (clueweb12)
	[ ] experiment comparing ZGraph && ASAP (although the machines are differnet)
	[ ] experiment about implementation details
		[ ] improvement on time-profiling && time-profiling performance (expected VS actual)
		[ ] improvemnt on error-profile && error-profile performance (expected VS actual)
		[ ] improvement on conditonal_sample_edge (binary search)
		[ ] improvement of numa-aware sub-partition
	[ ] auto-generator (two-level sampling) // that's bad... change to canonicality check solution ... 
		[ ] to show that auto-generated sampling programs is comparable with the hand-write ones 
		[ ] to compare our generator with the naive implementation  // without neighbour samling or quick-pattern aggregation
		[ ] (possible) degree assisted optimization
	[ ] fix some problem about (at-least one) predicate matching 
	[ ] complex predicate matching (results) 
		[ ] comparing the performance using complex predicate matching && those with (all, at-least) predicate matching (show that is comparable)

## Reference
[1] Anand, I., et al. "ASAP: Fast, Approximate Graph Pattern Mining at Scale." Proceedings of the USENIX conference. 2018.
# sampling-graph-mining

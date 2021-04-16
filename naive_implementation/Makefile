all: TriangleCounting ThreeChain ThreeChainApprox WeightedTriangleCounting AtLeastOneWeightedTriangleCounting FiveChain FourChain FiveStar FiveHouse

TriangleCounting: TriangleCounting.cpp
	g++ TriangleCounting.cpp -o TriangleCounting

ThreeChain: ThreeChain.cpp
	g++ ThreeChain.cpp -o ThreeChain

ThreeChainApprox: ThreeChainApprox.cpp
	g++ ThreeChainApprox.cpp -o ThreeChainApprox

WeightedTriangleCounting: WeightedTriangleCounting.cpp
	g++ WeightedTriangleCounting.cpp -o WeightedTriangleCounting

AtLeastOneWeightedTriangleCounting: AtLeastOneWeightedTriangleCounting.cpp
	g++ AtLeastOneWeightedTriangleCounting.cpp -o AtLeastOneWeightedTriangleCounting

FiveChain: FiveChain.cpp
	g++ FiveChain.cpp -o FiveChain

FourChain: FourChain.cpp
	g++ FourChain.cpp -o FourChain

#FourChainApprox: FourChainApprox.cpp
#	g++ FourChainApprox.cpp -o FourChainApprox

FiveStar: FiveStar.cpp
	g++ FiveStar.cpp -o FiveStar

FiveHouse: FiveHouse.cpp
	g++ FiveHouse.cpp -o FiveHouse -fopenmp -lnuma -std=c++11 -O3

clean:
	rm TriangleCounting ThreeChain ThreeChainApprox WeightedTriangleCounting AtLeastOneWeightedTriangleCounting FiveChain FourChain FiveStar FiveHouse

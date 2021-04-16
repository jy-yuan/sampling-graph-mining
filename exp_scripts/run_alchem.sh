for application in {"Triangle","ThreeMotifWithoutCache","FourChain","FiveStar"}
do
	mkdir result_alchem/${application}
	for graph in {"wiki-vote","youtube","patents","live-journal","twitter-2010"}
	do
		mkdir result_alchem/${application}/${graph}
		for i in {0..2}
		do
			echo "mpirun -n 8 -N 1 -hostfile ~/hosts ./build/${application} ~/zgraph_datasets/${graph} > result_alchem/${application}/${graph}/${i}.txt"
			mpirun -n 8 -N 1 -hostfile ~/hosts ./build/${application} ~/zgraph_datasets/${graph} > result_alchem/${application}/${graph}/${i}.txt
		done
	done
done

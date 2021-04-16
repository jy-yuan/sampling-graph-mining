#for application in {"FourChain","FiveStar"}
#for application in {"Triangle","ThreeMotifWithoutCache"}
application="Triangle2"
base_dir=result_triangle2
#do
	mkdir ${base_dir}/${application}
	#graph="twitter-2010"
	#for graph in {"wiki-vote","youtube","patents","live-journal","twitter-2010","friendster"}
	for graph in {"wiki-vote","youtube","patents","live-journal"}
	do
		mkdir ${base_dir}/${application}/${graph}
		for i in {0..2}
		do
			echo "mpirun --allow-run-as-root -n 6 -N 1 -hostfile ~/hosts ./build/${application} ~/ZGraph/zgraph_datasets/${graph} > ${base_dir}/${application}/${graph}/${i}.txt"
			mpirun --allow-run-as-root -n 6 -N 1 -hostfile ~/hosts ./build/${application} ~/ZGraph/zgraph_datasets/${graph} > ${base_dir}/${application}/${graph}/${i}.txt
		done
	done
#done

#application="FiveStar"
for nodes in {1..6}
do
	base_dir=result_node_${nodes}
	mkdir ${base_dir}
	for application in {"FourChain","FiveStar"}
	do
		mkdir ${base_dir}/${application}
		for graph in {"wiki-vote","youtube","patents","live-journal"}
		do
			mkdir ${base_dir}/${application}/${graph}
			for i in {0..2}
			do
				echo "mpirun --allow-run-as-root -n ${nodes} -N 1 -hostfile ~/hosts ./build/${application} ~/ZGraph/zgraph_datasets/${graph} > ${base_dir}/${application}/${graph}/${i}.txt"
				mpirun --allow-run-as-root -n ${nodes} -N 1 -hostfile ~/hosts ./build/${application} ~/ZGraph/zgraph_datasets/${graph} > ${base_dir}/${application}/${graph}/${i}.txt
			done
		done
	done
done

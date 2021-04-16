#application="FiveStar"
base_dir=result_thread_12
for application in {"FourChain","FiveStar"}
do
	mkdir ${base_dir}/${application}
	for graph in {"wiki-vote","youtube","patents","live-journal"}
	do
		mkdir ${base_dir}/${application}/${graph}
		for i in {0..2}
		do
			echo "mpirun --allow-run-as-root -n 6 -N 1 -hostfile ~/hosts ./build/${application} ~/ZGraph/zgraph_datasets/${graph} > ${base_dir}/${application}/${graph}/${i}.txt"
			mpirun --allow-run-as-root -n 6 -N 1 -hostfile ~/hosts ./build/${application} ~/ZGraph/zgraph_datasets/${graph} > ${base_dir}/${application}/${graph}/${i}.txt
		done
	done
done

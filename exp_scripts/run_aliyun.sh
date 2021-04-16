for application in {"Triangle","ThreeMotifWithoutCache"}
do
	mkdir result_aliyun/${application}
	graph=youtube
		mkdir result_aliyun/${application}/${graph}
		for i in {0..2}
		do
			echo "mpirun --allow-run-as-root -n 6 -N 1 -hostfile ~/hosts ./build/${application} /storage/datasets/${graph} > result_aliyun/${application}/${graph}/${i}.txt"
			mpirun --allow-run-as-root -n 6 -N 1 -hostfile ~/hosts ./build/${application} /storage/datasets/${graph} > result_aliyun/${application}/${graph}/${i}.txt
		done
done

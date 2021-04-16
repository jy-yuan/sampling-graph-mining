cd build 
for i in {3..6}:
do
	mpirun -n 4 -N 1 -host nerv5,nerv6,nerv7,nerv8 ./ChainMining ../datasets/wiki-vote/dataset/wiki-vote ${i}
	mpirun -n 4 -N 1 -host nerv5,nerv6,nerv7,nerv8 ./ChainMining ../datasets/youtube/dataset/youtube ${i}
	mpirun -n 4 -N 1 -host nerv5,nerv6,nerv7,nerv8 ./ChainMining ../datasets/soc-livejournal/dataset/live-journal ${i}
	mpirun -n 4 -N 1 -host nerv5,nerv6,nerv7,nerv8 ./ChainMining ../datasets/twitter/dataset/twitter-2010 ${i}
done
scancel 3497

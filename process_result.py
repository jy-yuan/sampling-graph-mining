import math

applications = [
       # "Triangle2",
       # "ThreeMotifWithoutCache",
       # "FourChain",
       # "FiveStar"
       "ChainMining"
        ]

graphs = [
        "wiki-vote",
        "youtube",
        "patents",
        "live-journal",
        #"twitter-2010",
	#"friendster"
        ]

for application in applications:
    for graph in graphs:
        results = []
        time = 0
        for i in range(3):
            d = "./result_elp_opt/%s/%s/%s.txt" % (application, graph, i)
            result = None
            with open (d, "r") as f:
                while True:
                    line = f.readline()
                    if line == None or len(line) == 0:
                        break
                    key = "estimated number of "
                    if application == "ThreeMotifWithoutCache":
                        key = "estimated number of 3-chains:"
                    if key in line:
                        tmp = line.strip().split(" ")[-1]
                        tmp = float(tmp)
                        if result == None:
                            result = tmp
                        else:
                            assert abs(result - tmp) < 1e-6
                    if "estimation," in line:
                        tmp = line.strip().split("\t")[-1]
                        tmp = float(tmp)
                        time += tmp
            results.append(result)
        time /= 3
        print "%s-%s: "%(application, graph), " ", time, " ", results
        #print "%s-%s: "%(application, graph), " %.3f" %(time / 1000.) 


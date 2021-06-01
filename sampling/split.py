# file like
# u1 v1
# u2 v2
# ....
# transfrom into several files with same from
# every node should contain every adjacent edge
# if not directed, double the edge

from tqdm import tqdm
import argparse

def genLines(u, vs, mapping):
    printlis = []
    vs = sorted(list(vs))
    for v in vs:
        edge = str(mapping[u]) + ' ' + str(mapping[v]) + '\n'
        printlis.append(edge)
    return printlis

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--num", help="num of slices", required=True)
    parser.add_argument("-t", "--target", help="target graph", required=True)
    args = parser.parse_args()
    n = int(args.num)
    target = args.target
    directed = False
    # default wiki-vote
    path = "../datasets/wiki-vote/dataset/wiki-Vote.txt"
    givennodes = 7115
    if target == 'fr':
        path = "../datasets/friendster/dataset/com-friendster.ungraph.txt"
        givennodes = 65608366
    elif target == 'tw':
        path = "../datasets/twitter/dataset/twitter-2010.txt"
        givennodes = 41652230
    elif target == 'lj':
        path = "../datasets/soc-livejournal/dataset/soc-LiveJournal1.txt"
        givennodes = 4847571
    elif target == 'pt':
        path = "../datasets/patent/dataset/cit-Patents.txt"
        givennodes = 3774768
    elif target == 'yt':
        path = "../datasets/youtube/dataset/com-youtube.ungraph.txt"
        givennodes = 1134890
    graph = {}
    with open(path, "r") as f:
        print("reading raw graph")
        while 1:
            lines = f.readlines(1 << 30)
            if not lines:
                break
            else:
                lines = tqdm(lines)
                for line in lines:
                    uv = line.split()
                    if not uv[0].isnumeric() or len(uv) != 2:
                        continue
                    u = int(uv[0])
                    v = int(uv[1])
                    if u == v:
                        continue
                    if u not in graph.keys():
                        graph[u] = set()
                    graph[u].add(v)
                    if not directed:
                        if v not in graph.keys():
                            graph[v] = set()
                        graph[v].add(u)

    keys = sorted(list(graph.keys()))
    maxnode = keys[-1]
    mapping = [0] * (maxnode+1)
    for i in range(len(keys)):
        mapping[keys[i]] = i
    size = len(keys) // n
    for i in range(n):
        edges_count = 0
        print("saving graph %d"%i)
        with open(str(i), "w") as f:
            for j in tqdm(range(i*size, (i+1)*size)):
                u = keys[j]
                f.writelines(genLines(u, graph[u], mapping))
                edges_count += len(graph[u])
            if i == n-1:
                for j in range((i+1)*size, len(keys)):
                    u = keys[j]
                    f.writelines(genLines(u, graph[u], mapping))
                    edges_count += len(graph[u])
            
        with open(str(i), "r+") as f:
            content = f.read()
            f.seek(0,0)
            f.write(str(len(keys)) + ' ' + str(edges_count) +'\n' + content)

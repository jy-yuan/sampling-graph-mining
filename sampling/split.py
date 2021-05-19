# file like
# u1 v1
# u2 v2
# ....
# transfrom into several files with same from
# every node should contain every adjacent edge
# if not directed, double the edge

from tqdm import tqdm

def genLines(u, vs):
    printlis = []
    for v in vs:
        edge = str(u) + ' ' + str(v) + '\n'
        printlis.append(edge)
    return printlis

if __name__ == "__main__":
    n = 2
    directed = True
    path = "../datasets/wiki-vote/dataset/wiki-Vote.txt"
    graph = {}
    with open(path, "r") as f:
        lines = tqdm(f.readlines())
        for line in lines:
            uv = line.split()
            if not uv[0].isnumeric() or len(uv) != 2:
                continue
            u = int(uv[0])
            v = int(uv[1])
            if u not in graph.keys():
                graph[u] = set()
            graph[u].add(v)
            if not directed:
                if v not in graph.keys():
                    graph[v] = set()
                graph[v].add(u)

    keys = sorted(list(graph.keys()))
    maxnode = keys[-1]
    size = len(keys) // n
    for i in range(n):
        edges_count = 0
        with open(str(i) + ".graph", "w") as f:
            for j in range(i*size, (i+1)*size):
                u = keys[j]
                f.writelines(genLines(u, graph[u]))
                edges_count += len(graph[u])
            if i == n-1:
                for j in range((i+1)*size, len(keys)):
                    u = keys[j]
                    f.writelines(genLines(u, graph[u]))
                    edges_count += len(graph[u])
            
        with open(str(i) + ".graph", "r+") as f:
            content = f.read()
            f.seek(0,0)
            f.write(str(maxnode+1) + ' ' + str(edges_count) +'\n' + content)
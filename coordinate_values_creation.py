import networkx as nx
filename = "pmed39.txt"

def readGraph(filename):
    input_file = open(filename, 'r')
    lines = input_file.readlines()
    input_file.close()

    header = lines.pop(0).strip().split(' ')
    p = int(header[2])
    G = nx.Graph()
    for line in lines:
        edge = line.strip().split(' ')
        i = int(edge[0]) - 1
        j = int(edge[1]) - 1
        cost = int(edge[2])
        edge = (i, j)
        G.add_edge(*edge, weight=cost)

    distances = nx.floyd_warshall(G)
    #print(distances)
    #print(distances[899])
    output_file = open('B900_coord.txt', 'w')
    coordinate_values = str(distances[899])
    output_file.write(coordinate_values)
    output_file.close()


    input_file.close()
    return G, p, distances

readGraph("pmed39.txt")
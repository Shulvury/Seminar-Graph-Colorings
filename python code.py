import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import interp1d


def create_randomized_graph(a: float, n: int):
    graph = np.zeros((n, n), dtype=bool)
    count_edges = 0
    for i in range(n):
        for j in range(i + 1, n):
            if a > np.random.uniform(0, 1):
                count_edges += 1
                graph[i][j] = True
                graph[j][i] = True
    print("number of edges: ", count_edges)
    print("overall prob: ", count_edges / (n * (n - 1) / 2))
    # print(graph)
    return graph, n, count_edges


def graph_to_dimacs(graph: np.array, n, m):
    text = "p edge " + str(n) + ' ' + str(m) + '\n'
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if graph[i-1][j-1]:
                text += "e " + str(i) + ' ' + str(j) + '\n'
                #print("e ", i + 1, j + 1)
    return text


def create_graph_at_most_k(n: int, k: int):
    graph = np.zeros((n, n), dtype=bool)
    deg_list = np.zeros(n, dtype=int)
    count_edges = 0
    for l in range(n * k):
        if l > n * k / 2:
            break
        # i, j = np.random.uniform(0,n, 2)
        # i = int(i)
        # j = int(j)
        i = np.random.randint(0, n)
        j = np.random.randint(0, n)
        if i == j:
            continue
        elif deg_list[i] >= k:
            continue
        elif deg_list[j] >= k:
            continue
        elif not graph[i][j]:
            count_edges += 1
            deg_list[i] += 1
            deg_list[j] += 1
            graph[i][j] = True
            graph[j][i] = True

    # print("number of edges: ", count_edges)
    # print("overall prob: ", count_edges / (n * (n - 1) / 2))
    # print(graph)
    return graph, n, count_edges


def create_graph_at_least_k(n: int, k: int):
    graph = np.zeros((n, n), dtype=bool)
    deg_list = np.zeros(n, dtype=int)
    count_edges = 0
    #for i in range(0, n):
    #    for j in range(0, n):
    #        if np.random.uniform(0, 1) < (k / n) * 1.5:
    #            graph[i][j] = True
    #            graph[j][i] = True
    #            deg_list[i] += 1
    #            deg_list[j] += 1

    unvisited_vertices = list(range(0, n))
    for count in range(0, n):
        i = unvisited_vertices[np.random.randint(0, (n - count))]
        unvisited_vertices.remove(i)

        if deg_list[i] < k:
            missing_edges = list()
            for j in range(0, n):
                if i == j:
                    continue
                if not graph[i][j]:
                    missing_edges.append(j)
            missing_edges_len = len(missing_edges)
            while deg_list[i] < k:
                j = missing_edges[np.random.randint(0, missing_edges_len)]
                graph[i][j] = True
                graph[j][i] = True
                deg_list[i] += 1
                deg_list[j] += 1
                count_edges += 1
    count_edges = 0
    for i in range(n):
        for j in range(i+1, n):
            if graph[i][j]:
                count_edges += 1
    # print("number of edges: ", count_edges)
    # print("overall prob: ", count_edges / (n * (n - 1) / 2))
    # print(graph)
    # print((deg_list[:] - k))
    return graph, n, count_edges


def create_flat_graph_equidistant(n: int, q: int, a: float):
    graph = np.zeros((n, n), dtype=bool)
    count_edges = 0

    ind_set_list = []
    set_size = math.ceil(n / q)

    for i in range(0, n, set_size):
        vert_set = []
        if i + set_size < n:
            for k in range(i, i + set_size):
                vert_set.append(k)
            ind_set_list.append(vert_set)
        else:
            for k in range(i, n):
                vert_set.append(k)
            ind_set_list.append(vert_set)

    #print(ind_set_list)


    for i in range(len(ind_set_list)):
        min_vert = min(ind_set_list[i])
        max_vert = max(ind_set_list[i])

        for j in range(min_vert, max_vert+1):
            for k in range(max_vert+1, n):
                if k >= n:
                    continue
                if a > np.random.uniform(0, 1):
                    count_edges += 1
                    graph[k][j] = True
                    graph[j][k] = True

    print("number of edges: ", count_edges)
    print("overall prob: ", count_edges / (n * (n - 1) / 2))
    # print(graph)
    return graph, n, count_edges


def create_flat_inhomogeneous(n: int, size_limit: float, a: float):
    graph = np.zeros((n, n), dtype=int)
    count_edges = 0

    ind_set_list = []
    num_of_sets = 0
    i = 0

    while i < n:
        upper_bound = (n - i + 1) * size_limit
        if upper_bound < 10:
            upper_bound = 15
        q = np.random.randint(5, upper_bound)
        vert_set = []
        if i + q - 1 <= n:
            for k in range(i, i+q):
                vert_set.append(k)
            ind_set_list.append(vert_set)
        else:
            for k in range(i, n):
                vert_set.append(k)
            ind_set_list.append(vert_set)
        i = i + q
        num_of_sets += 1

    #print(ind_set_list)
    i = 0
    #len_list = []
    #for k in ind_set_list:
    #    len_list.append(len(k))
    #len_list.sort()
    #print(len(len_list), len_list)

    for i in range(len(ind_set_list)):
        min_vert = min(ind_set_list[i])
        max_vert = max(ind_set_list[i])

        for j in range(min_vert, max_vert + 1):
            for k in range(max_vert + 1, n):
                if k >= n:
                    continue
                if a > np.random.uniform(0, 1):
                    count_edges += 1
                    graph[k][j] = True
                    graph[j][k] = True

    #print("number of edges: ", count_edges)
    #print("overall prob: ", count_edges / (n * (n - 1) / 2))
    #print("number of sets: ", num_of_sets)
    #for v in ind_set_list:
    #    print(len(v))

    # print(graph)
    return graph, n, count_edges


def write_graphs(graph, n, m, _name):
    myfile = open(_name, "w")
    myfile.write(graph_to_dimacs(graph, n, m))
    myfile.close()


# G1 = create_randomized_graph(A, N)
# graph_to_dimacs(G1[0], G1[1], G1[2])

#G2 = create_graph_at_most_k(N, K)
# print(G2[1], G2[2])
# print(G2[3])
#graph_to_dimacs(G2[0], G2[1], G2[2])

#G3 = create_graph_at_least_k(N, K)
#print(G3)

#K = 100
#N = 800
#i = 0
#for A in np.linspace(0.05, 1, num=19, endpoint=False):
#    i = i+1
#    for n_files in range(10):
#        file_name = "random_graph_" + str(i) + "_A_" + str(n_files) + ".txt"
#        G = create_randomized_graph(A, N)
#        write_graphs(G[0], G[1], G[2], file_name)

mode = "graph"
if mode == "gen":
    N = 500
    Q = 19
    count = 0
    Quot_size = 0.3
    for p in np.linspace(0.05, 1, num=19, endpoint=False):
        count = count+1
        num_of_instances = 20
        for n_files in range(num_of_instances):
            file_name = "flat_graph_" + str(count) + "_A_" + str(n_files) + ".txt"
            #file_name = "random_graph_" + str(count) + "_A_" + str(n_files) + ".txt"
            #file_name = "inhom_graph_" + str(count) + "_A_" + str(n_files) + ".txt"
            #G = create_randomized_graph(p, N)
            G = create_flat_graph_equidistant(N, Q, p)
            #G = create_flat_inhomogeneous(N, Quot_size, p)
            write_graphs(G[0], G[1], G[2], file_name)

elif mode == "graph":
    #file = open('Rand N100.log', 'r')
    file = open('Inhom.log', 'r')
    for alg_num in range(6):

        tries = 20

        lines = list()
        j = 0
        #temp = file.readline()
        for j in range(tries*19):
            temp = file.readline()
            lines.append(temp.split())
            #print(j, lines[j])

            lines[j][0] = int(lines[j][0])
            lines[j][1] = int(lines[j][1])
            lines[j][2] = int(lines[j][2])
            #print(lines[i])

        temp = file.readline()
        average_of_lines_1 = list()
        #average_of_lines_2 = list()
        average_of_lines_3 = list()

        deviation_list_col = 0
        deviation_list_ops = 0

        for i in range(0, tries*19, tries):
            avrg_1 = 0
            avrg_3 = 0

            for j in range(tries):
                avrg_1 += lines[i+j][0]
                #avrg_2 += lines[i+j][1]
                avrg_3 += lines[i+j][2]

            avrg_1 = int(avrg_1 / tries)
            avrg_3 = int(avrg_3 / tries)

            average_of_lines_1.append(avrg_1)
            average_of_lines_3.append(avrg_3)

        slice_stage = 8
        for j in range(tries):
            deviation_list_col += math.fabs(average_of_lines_1[slice_stage] - int(lines[slice_stage*tries + j][0]))
            deviation_list_ops += math.fabs(average_of_lines_3[slice_stage] - int(lines[slice_stage*tries + j][2]))
        deviation_list_ops = deviation_list_ops / tries
        deviation_list_col = deviation_list_col / tries

        print("Col deviation: ", deviation_list_col)
        print("Ops deviation: ", deviation_list_ops, '\n\n')

        palette = []

        plt.figure(1)
        plt.title("n = 500\n Inhomogeneous")
        #plt.plot(np.linspace(0.05, 0.6, num=11, endpoint=False), average_of_lines_1[0:11])
        #plt.plot(np.linspace(0.05, 1, num=len(average_of_lines_1), endpoint=False), average_of_lines_1[:])
        plt.plot(np.linspace(0.05, 1, num=len(average_of_lines_1), endpoint=False), average_of_lines_1[:])
        #plt.axhline(y=25, color='black', linewidth='0.5')
        plt.xlabel("Edge probability")
        plt.ylabel("Colors")
        plt.legend(['Greedy', 'LDO', 'DSatur', 'IDO', 'WP', 'RLF'])
        plt.ylim(0, 100)
        plt.grid(True)


        #plt.figure(2)
        #plt.title("n = 500\n Flat Q25")
        #smooth_data = interp1d(np.linspace(0.05, 0.95, num=len(average_of_lines_1), endpoint=True), average_of_lines_1[:],
        #                       kind="cubic")
        #new_x = np.linspace(0.05, 0.95, 500, endpoint=True)
        #plt.plot(new_x, smooth_data(new_x))
        #plt.xlabel("Edge probability")
        #plt.ylabel("Colors")
        #plt.yscale('log')
        #plt.legend(['Greedy', 'LDO', 'DSatur', 'IDO', 'WP', 'RLF'])
        #plt.grid(True)

        plt.figure(3)
        plt.title("n = 500\n Q = 19")
        #plt.title("n = 800\n Randomized")
        #plt.plot(np.linspace(0.05, 0.6, num=11, endpoint=False), average_of_lines_3[0:11])
        #plt.plot(np.linspace(0.05, 1, num=len(average_of_lines_3), endpoint=False), average_of_lines_3[:])
        plt.plot(np.linspace(0.05, 1, num=len(average_of_lines_3), endpoint=False), average_of_lines_3[:])
        plt.xlabel("Edge probability")
        plt.ylabel("Number of operations")
        plt.legend(['Greedy', 'LDO', 'DSatur', 'IDO', 'WP', 'RLF'])
        plt.yscale('log')
        plt.grid(True)

    file.close()
    plt.show()

elif mode == "try":

    N = 500
    Q = 3
    p = 0.5
    Cap_size = 0.3

    sizes_list = []
    mean_size_val = 0
    mean_set_val = 0

    runs = 100

    for i in range(runs):
        G = create_flat_inhomogeneous(N, Cap_size, p)
        mean_set_val += len(G[3])
        for j in G[3]:
            sizes_list.append(j)
            mean_size_val += j

    mean_set_val = mean_set_val / runs
    mean_size_val = mean_size_val / len(sizes_list)

    dev = 0
    Disp = 0
    for i in range(len(sizes_list)):
        dev += math.fabs(sizes_list[i] - mean_size_val)
        Disp += (sizes_list[i] - mean_size_val)**2
    dev = dev / len(sizes_list)
    Disp = Disp / len(sizes_list)
    print("Mean color sets value: ", mean_set_val)
    print("Mean size of set value: ", mean_size_val)
    print("Mean deviation of set size: ", dev)
    print("Dispersion: ", Disp)
    print("Coefficient varation: ", math.sqrt(Disp)/mean_size_val)
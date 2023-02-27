#include <cstring>
#include <string>
#include <climits>
#include <fstream>
#include <iostream>
#include <vector>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include <tuple>
#include <set>
#include <unordered_set>
#include <cmath>
using namespace std;

//Global Variables-----------------------------------------------------------------------
unsigned long long numConfChecks;
int verbose;

//Struct used for holding the graph------------------------------------------------------
struct Graph {
	int n = 0;
	int m = 0;
	int maxDeg = 0;
	vector<vector<int>> AList;
	vector<int> deg;
};

//Struct used in conjunction with sorting by degree--------------------------------------
struct degItem {
	int deg;
	int vertex;
};
struct maxDeg {
	bool operator() (const degItem& lhs, const degItem& rhs) const {
		//Compares two degItems by degree, then vertex label
		numConfChecks += 2;
		if (lhs.deg > rhs.deg) return true;
		numConfChecks += 2;
		if (lhs.deg < rhs.deg) return false;
		//if we are here we know that lhs.deg == rhs.deg. Our choice can be arbitrary
		if (lhs.vertex > rhs.vertex) return true;
		else return false;
	}
};

//Struct used in conjunction with the DSatur priority queue-------------------------------
struct satItem {
	int sat;
	int deg;
	int vertex;
};
struct maxSat {
	bool operator() (const satItem& lhs, const satItem& rhs) const {
		//Compares two satItems sat deg, then degree, then vertex label
		numConfChecks += 2;
		if (lhs.sat > rhs.sat) return true;
		numConfChecks += 2;
		if (lhs.sat < rhs.sat) return false;
		//if we are we know that lhs.sat == rhs.sat
		numConfChecks += 2;
		if (lhs.deg > rhs.deg) return true;
		numConfChecks += 2;
		if (lhs.deg < rhs.deg) return false;
		//if we are here we know that lhs.sat == rhs.sat and lhs.deg == rhs.deg. Our choice can be arbitrary
		if (lhs.vertex > rhs.vertex) return true;
		else return false;
	}
};

//Procedures used by all algorithms------------------------------------------------------
void logAndExit(string s) {
	//Writes message s to screen and log file and then exits the program
	ofstream resultsLog("resultsLog.log", ios::app);
	resultsLog << s;
	cout << s;
	resultsLog.close();
	//exit(1);
}
Graph readInputFile(string fname) {
	//Reads a DIMACS format file and return the corresponding Graph struct
	ifstream inStream;
	inStream.open(fname.c_str());
	if (inStream.fail()) logAndExit("Error. Unrecognized argument/filename <" + fname + "> in command.\n");
	char c, str[1000];
	int line = 0, u, v, edgeCnt = 0;
	vector<vector<bool>> A;
	Graph G;
	try {
		while (!inStream.eof()) {
			line++;
			inStream.get(c);
			if (inStream.eof()) break;
			switch (c) {
			case 'c':
				//ignore a comment line
				inStream.putback('c');
				inStream.get(str, 999, '\n');
				break;
			case 'p':
				//read the parameter line of the file and set up the adjacency matrix
				inStream.get(c);
				inStream.getline(str, 999, ' ');
				if (G.n != 0 || G.m != 0) logAndExit("Invalid input file. Line starting with 'p' at line " + to_string(line) + " defined more than once.\n");
				if (strcmp(str, "edge")) logAndExit("Invalid input file. Problem at line " + to_string(line) + ". No 'edge' keyword found.\n");
				inStream >> G.n >> G.m;
				A.clear();
				A.resize(G.n, vector<bool>(G.n, false));
				for (u = 0; u < G.n; u++) A[u][u] = true;
				break;
			case 'e':
				//Read an edge
				inStream >> u >> v;
				if (u < 1 || u > G.n || v < 1 || v > G.n || u == v) logAndExit("Invalid input file. Problem at line " + to_string(line) + ". Invalid edge.\n");
				if (!A[u - 1][v - 1]) edgeCnt++;
				else logAndExit("Invalid input file. Problem at line " + to_string(line) + ". Edge defined previously.\n");
				A[u - 1][v - 1] = true;
				A[v - 1][u - 1] = true;
				break;
			default:
				logAndExit("Invalid input file. Problem at line " + to_string(line) + "\n");
			}
			inStream.get();
		}
		//Finished reading the file
		inStream.close();
		if (edgeCnt != G.m) logAndExit("Invalid input file. Number of read edges does not equal number specified at the top of the file.\n");
	}
	catch (...) {
		logAndExit("Invalid input file. Unidentified error near line " + to_string(line) + ".\n");
	}
	inStream.close();
	//Check to see if there are no edges. If so, exit straight away
	if (G.m <= 0) {
		logAndExit("Graph has no edges. Optimal solution is obviously using one colour.\n");
	}
	//Now use the adjacency matrix to construct the graph G
	G.deg.clear();
	G.AList.clear();
	G.deg.resize(G.n, 0);
	G.AList.resize(G.n, vector<int>());
	for (u = 0; u < G.n; u++) {
		for (v = 0; v < G.n; v++) {
			if (A[u][v] && u != v) {
				G.AList[u].push_back(v);
				G.deg[u]++;
			}
		}
	}
	G.maxDeg = *max_element(G.deg.begin(), G.deg.end());
	return(G);
}
void prettyPrintSolution(vector<vector<int>>& S) {
	int i, count = 0, col;
	cout << "\n\n";
	for (col = 0; col < S.size(); col++) {
		cout << "C-" << col << "\t= {";
		if (S[col].size() == 0) cout << "empty}\n";
		else {
			for (i = 0; i < S[col].size() - 1; i++) cout << S[col][i] << ", ";
			cout << S[col][S[col].size() - 1] << "}\n";
			count = count + S[col].size();
		}
	}
	cout << "Total Number of Nodes = " << count << endl;
}

//Functions for Greedy and DSatur Algorithm----------------------------------------------
int getFirstFeasCol(Graph& G, int v, vector<int>& c, vector<bool>& used) {
	int i;
	for (int u : G.AList[v]) {
		if (c[u] != -1) used[c[u]] = true;
	}
	for (i = 0; i < used.size(); i++) {
		if (used[i] == false) break;
	}
	for (int u : G.AList[v]) {
		if (c[u] != -1) used[c[u]] = false;
	}
	numConfChecks += G.deg[v] + G.deg[v];
	return i;
}
int greedycol(Graph& G, bool sortByDegree) {
	int i;
	vector<bool> used(G.maxDeg + 1, false);
	vector<int> c(G.n, -1), perm;
	if (sortByDegree) {
		//Sort vertices by degree in O(n lg n) time
		set<degItem, maxDeg> L;
		numConfChecks += G.n;
		for (i = 0; i < G.n; i++) L.insert({ G.deg[i], i });
		for (degItem el : L) perm.push_back(el.vertex);
	}
	else {
		//Shuffle the vertices in O(n) time
		for (i = 0; i < G.n; i++) perm.push_back(i);
		random_shuffle(perm.begin(), perm.end());
	}
	//Do the greedy algorithm using perm
	for (int v : perm) {
		i = getFirstFeasCol(G, v, c, used);
		c[v] = i;
	}
	return *max_element(c.begin(), c.end());
}
int DSatur(Graph& G) {
	int u, i;
	vector<bool> used(G.maxDeg + 1, false);
	vector<int> c(G.n), d(G.n);
	vector<set<int>> adjCols(G.n);
	set<satItem, maxSat> Q;
	set<satItem, maxSat>::iterator maxPtr;
	//Initialise the the data structures. These are a (binary heap) priority queue, a set of colours adjacent to each uncoloured vertex (initially empty)
	//and the degree d(v) of each uncoloured vertex in the graph induced by uncoloured vertices
	numConfChecks += G.n;
	for (u = 0; u < G.n; u++) {
		c[u] = -1;
		d[u] = G.deg[u];
		adjCols[u] = set<int>();
		Q.emplace(satItem{ 0, d[u], u });
	}
	//DSatur algorithm
	while (!Q.empty()) {
		//Get the vertex u with highest saturation degree, breaking ties with d. Remove it from the priority queue and colour it
		numConfChecks++;
		maxPtr = Q.begin();
		u = (*maxPtr).vertex;
		Q.erase(maxPtr);
		i = getFirstFeasCol(G, u, c, used);
		c[u] = i;
		//Update the saturation degrees and d-value of all uncoloured neighbours; hence modify their corresponding elements in the priority queue
		numConfChecks += G.deg[u];
		for (int v : G.AList[u]) {
			if (c[v] == -1) {
				Q.erase({ int(adjCols[v].size()), d[v], v });
				adjCols[v].insert(i);
				d[v]--;
				Q.emplace(satItem{ int(adjCols[v].size()), d[v], v });
			}
		}
	}
	return *max_element(c.begin(), c.end());
}

//Functions for Greedy IS algorithm------------------------------------------------------
Graph relabelVertices(Graph& G, vector<int>& P) {
	//Create a copy of G with vertices relabelled according to P
	Graph H;
	H.n = G.n, H.m = G.m, H.maxDeg = G.maxDeg;
	H.AList.resize(G.n, vector<int>());
	H.deg.resize(G.n, 0);
	for (int u = 0; u < G.n; u++) {
		numConfChecks += G.deg[u];
		for (int v : G.AList[u]) {
			H.AList[P[u]].push_back(P[v]);
			H.deg[P[u]]++;
		}
	}
	return H;
}
void updateSets(Graph& G, unordered_set<int>& X, unordered_set<int>& Y, vector<int>& c, int u) {
	//Remove u from X (it is now coloured) and move all uncoloured neighbours of u from X to Y
	X.erase(u);
	numConfChecks += log2(X.size() + 1);
	for (int v : G.AList[u]) {
		if (c[v] == -1) {
			X.erase(v);
			Y.insert(v);
			numConfChecks += log2(Y.size());
		}
	}
}
int greedyIS(Graph& G) {
	//Create a graph H. This is a copy of G with vertices randomly relabelled.
	vector<int> P;
	//for (int v = 0; v < G.n; v++) P[v] = v;
	//random_shuffle(P.begin(), P.end());

	set<degItem, maxDeg> L;
    numConfChecks += G.n;
    for (int j = 0; j < G.n; j++) L.insert({ G.deg[j], j });
    for (degItem el : L) P.push_back(el.vertex);

	Graph H = relabelVertices(G, P);
	//Now colour H
	vector<int> c(H.n, -1);
	unordered_set<int> X, Y;
	for (int v = 0; v < H.n; v++) X.insert(v);
	int i = 0, u;
	while (!X.empty()) {
		//Constructing colour class i.
		while (!X.empty()) {
			u = *X.begin();
			c[u] = i;
			updateSets(H, X, Y, c, u);
		}
		X.swap(Y);
		i++;
	}
	//Convert the vertex labels used in H to correspond to a solution for G.
	//vector<int> final(G.n);
	//for (int v = 0; v < G.n; v++) final[v] = c[P[v]];
	return *max_element(c.begin(), c.end());
}

//Functions for RLF Algorithm------------------------------------------------------------
void populateNeighbourArrays(Graph& G, unordered_set<int>& X, vector<int>& NInX, vector<int>& NInY) {
	for (int u : X) {
		NInX[u] = 0;
		NInY[u] = 0;
	}
	for (int u : X) {
		numConfChecks += G.deg[u];
		for (int v : G.AList[u]) {
			if (X.count(v) == 1) NInX[u]++;
		}
	}
}
int chooseFirstVertex(unordered_set<int>& X, vector<int>& NInX) {
	//Select the vertex in (non-empty) X that has the maximum number of neighbours in X
	int v = -1, max = -1;
	for (int u : X) {
		if (NInX[u] > max) {
			max = NInX[u];
			v = u;
		}
	}
	return v;
}
int chooseNextVertex(unordered_set<int>& X, vector<int>& NInY, vector<int>& NInX) {
	//Select vertex in (non-empty) X with max neighbours in Y; break ties according to minimum neighbours within X
	int v = -1, max = -1, min = INT_MAX;
	for (int u : X) {
		if ((NInY[u] > max) || (NInY[u] == max && NInX[u] < min)) {
			max = NInY[u];
			min = NInX[u];
			v = u;
		}
	}
	return(v);
}
void updateSetsAndDegrees(Graph& G, unordered_set<int>& X, unordered_set<int>& Y, vector<int>& NInX, vector<int>& NInY, unordered_set<int>& D2, vector<int>& c, int u) {
	//Remove u from X (it is now coloured)
	X.erase(u);
	//Move all uncoloured neighbours of u from X to Y
	numConfChecks += G.deg[u];
	for (int v : G.AList[u]) {
		if (c[v] == -1) {
			X.erase(v);
			Y.insert(v);
		}
	}
	//The remaining parts of this procedure now recalculate the contets of NinX and NinY. First calculate a set D2 of all uncoloured vertites within distance two of u.
	D2.clear();
	numConfChecks += G.deg[u];
	for (int v : G.AList[u]) {
		if (c[v] == -1) {
			D2.insert(v);
			numConfChecks += G.deg[v];
			for (int w : G.AList[v]) {
				if (c[w] == -1) {
					D2.insert(w);
				}
			}
		}
	}
	//For each vertex v in D2, now recalculate the number of (uncoloured) neighbours in X and Y
	for (int v : D2) {
		NInX[v] = 0;
		NInY[v] = 0;
		numConfChecks += G.deg[v];
		for (int w : G.AList[v]) {
			if (c[w] == -1) {
				if (X.count(w) == 1) NInX[v]++;
				else if (Y.count(w) == 1) NInY[v]++;
			}
		}
	}
}
int RLF(Graph& G) {
	vector<int> NInX(G.n), NInY(G.n), c(G.n, -1);
	unordered_set<int> X, Y, D2;
	for (int v = 0; v < G.n; v++) X.insert(v);
	int u, i = 0;
	while (!X.empty()) {
		//Constructing colour class i. First calculate the contents of the neighbours arrays, then colour the vertex u in X that has the most neighbours in X
		populateNeighbourArrays(G, X, NInX, NInY);
		u = chooseFirstVertex(X, NInX);
		c[u] = i;
		updateSetsAndDegrees(G, X, Y, NInX, NInY, D2, c, u);
		while (!X.empty()) {
			//Colour the vertex u in X that has the largest number of neighbours in Y. Break ties according to the minimum neighbours within X
			u = chooseNextVertex(X, NInY, NInX);
			c[u] = i;
			updateSetsAndDegrees(G, X, Y, NInX, NInY, D2, c, u);
		}
		//Have finished constructing colour i
		X.swap(Y);
		i++;
	}
	return *max_element(c.begin(), c.end());
}

//Incidence Coloring---------------------------------------------------------------------
int incidensecol(Graph& G){

    vector<int> c(G.n, -1); //color
    vector<int> incidence_degree(G.n, 0);
    vector<bool> used(G.n, false);
    vector<bool> used_colors(G.maxDeg + 1, false);

    // find vertex with largest degree
    int current_point = 0;
    for(int i = 0; i < G.n; i++){
        if(G.deg[i] = G.maxDeg){
            current_point = i;
            break;
        }
    }
    numConfChecks += G.n;

    //initialize first value
    c[current_point] = 0;
    used[current_point] = true;
    int inc_deg = 0;                        //maximum inc deg in current run
    for(int u : G.AList[current_point]){    //update inc deg list
        incidence_degree[u] += 1;
    }
    numConfChecks += G.AList[current_point].size();

    //start algorithm
    vector<int> inc_deg_list;
    int n = G.n - 1;
    for(int i = 0; i < n; i++){
        //find biggest inc deg
        for(int j = 0; j < G.n; j++){
            if(used[j]) continue;
            if(inc_deg < incidence_degree[j]){
               inc_deg = incidence_degree[j];
            }
            numConfChecks += 1;
        }

        //find all vertices with biggest inc deg
        for(int j = 0; j < G.n; j++){
            if(used[j]) continue;
            if(inc_deg == incidence_degree[j]){
               inc_deg_list.push_back(j);
            }
        }

        int deg_comp = -1;
        if(inc_deg_list.size() > 1){
            //break ties with degree
            for(int j: inc_deg_list){
                if(deg_comp < G.deg[j]){
                    deg_comp = G.deg[j];
                    current_point = j;
                }
            }
        }
        else{
            current_point = inc_deg_list[0];
        }

        numConfChecks += inc_deg_list.size();
        inc_deg_list.clear();
        deg_comp = -1;
        inc_deg = -1;


        //color current point with smallest feasible color
        int clr = getFirstFeasCol(G, current_point, c, used_colors);
        c[current_point] = clr;

        used[current_point] = true;

        //update inc deg list
        set<int> pallete;
        for(int j = 0; j < G.n; j++){
            if( used[j] ) continue;
            for(int v: G.AList[j]){
                if(c[v] == -1) continue;
                pallete.insert(c[v]);
                numConfChecks += log2(pallete.size());
            }
            incidence_degree[j] = pallete.size();
            pallete.clear();
        }

    }

    return *max_element(c.begin(), c.end());
}

//Welsh-Powell---------------------------------------------------------------------------


//Main algorithm-------------------------------------------------------------------------
int main(int argc, char** argv) {

    const int algorithms_number = 6;
    //const int number_of_runs = 20;
    const int first_num = 19;
    const int second_num = 19;
	for(int algChoice = 1; algChoice <= algorithms_number; algChoice++){

        //logAndExit("Algorithm " + std::to_string(algChoice) + "\n\n");

        for(int i = 1; i <= first_num; i++){
            for(int j = 0; j <= second_num; j++){
                //string _file_name = "random_graph_" + std::to_string(i) + "_A_" + std::to_string(j) + ".txt";
                //string _file_name = "inhom_graph_" + std::to_string(i) + "_A_" + std::to_string(j) + ".txt";
                string _file_name = "flat_graph_" + std::to_string(i) + "_A_" + std::to_string(j) + ".txt";
                //string _file_name = "graph.txt";
                numConfChecks = 0;
                verbose = 0;
                Graph G = readInputFile(_file_name);

                int c;

                clock_t runStart = clock();

                if      (algChoice == 1)	c = greedycol(G, false);    // Greedy
                else if (algChoice == 2)	c = greedycol(G, true);     // LDO
                else if (algChoice == 3)	c = DSatur(G);              // DSatur
                else if (algChoice == 4)	c = incidensecol(G);        // IDO
                else if (algChoice == 5)	c = greedyIS(G);            // Welsh Powell
                else if (algChoice == 6)	c = RLF(G);                 // RLF

                int duration = (int)(((clock() - runStart) / double(CLOCKS_PER_SEC)) * 1000);

                if      (algChoice == 1)	logAndExit(to_string(c+1) + "\t" + to_string(duration) + "\t" + to_string(numConfChecks) + "\n");
                else if (algChoice == 2)	logAndExit(to_string(c+1) + "\t" + to_string(duration) + "\t" + to_string(numConfChecks) + "\n");
                else if (algChoice == 3)	logAndExit(to_string(c+1) + "\t" + to_string(duration) + "\t" + to_string(numConfChecks) + "\n");
                else if (algChoice == 4)	logAndExit(to_string(c+1) + "\t" + to_string(duration) + "\t" + to_string(numConfChecks) + "\n");
                else if (algChoice == 5)	logAndExit(to_string(c+1) + "\t" + to_string(duration) + "\t" + to_string(numConfChecks) + "\n");
                else if (algChoice == 6)	logAndExit(to_string(c+1) + "\t" + to_string(duration) + "\t" + to_string(numConfChecks) + "\n");

                numConfChecks = 0;
            }
        }
        logAndExit("\n");

	}

	return(0);
}

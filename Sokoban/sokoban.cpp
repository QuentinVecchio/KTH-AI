#include <iostream>
#include <climits>
#include <vector>
#include <deque>

using namespace std;

int main() {
    // Get the map from input
    vector<string> map;
    string m;
    while (getline(cin, m)) {
        map.push_back(m);
    }

    // Create the maze and label all vortex
    vector<vector<int> > matrice;
    int start = -1;
    vector<int> ends;
    int currentVortex = 0;
    for(int i=0;i<map.size();i++) {
        vector<int> matr;
        for(int j=0;j<map[i].size();j++) {
            switch(map[i][j]) {
                case ' ':
                    matr.push_back(currentVortex);
                    currentVortex++;
                break;
                case '@': 
                    start = currentVortex;
                    matr.push_back(currentVortex);
                    currentVortex++;
                break;
                case '+':
                    cout << endl;
                    return 0;
                break;
                case '.': 
                    ends.push_back(currentVortex);
                    matr.push_back(currentVortex);
                    currentVortex++;
                break;
                default:
                matr.push_back(-1);
            }
        }
        matrice.push_back(matr);
    }

    if(start == -1 || ends.size() == 0) {
        cout << "no path" << endl;
        return 0;
    }

    /*for(int i=0;i<matrice.size();i++) {
        for(int j=0;j<matrice[i].size();j++) {
            cout << matrice[i][j] << " ";
        }
        cout << endl;
    }*/

    // Create edges between vortex
    vector<vector<int> > adjList;
    for(int i=0;i<matrice.size();i++) {
        for(int j=0;j<matrice[i].size();j++) {
            if(matrice[i][j] >= 0) {
                vector<int> adj;
                if((i-1) >= 0 && matrice[i-1][j] >= 0 && j < matrice[i+1].size()) {
                    adj.push_back(matrice[i-1][j]);
                }
                if((i+1) < matrice.size() && matrice[i+1][j] >= 0 && j < matrice[i+1].size()) {
                    adj.push_back(matrice[i+1][j]);
                }
                if((j-1) >= 0 && matrice[i][j-1] >= 0) {
                    adj.push_back(matrice[i][j-1]);
                }
                if((j+1) < matrice[i].size() && matrice[i][j+1] >= 0) {
                    adj.push_back(matrice[i][j+1]);
                }
                adjList.push_back(adj);
            }
        }
    }

    /*for(int i=0;i<adjList.size();i++) {
        cout << i << " : ";
        for(int j=0;j<adjList[i].size();j++) {
            cout << adjList[i][j] << " ";
        }
        cout << endl;
    }*/

    // BFS on the graph
    vector<int> d;
    vector<int> pi;
    vector<int> color; // 0: WHITE, 1: GRAY; 2: BLACK
    deque<int> F;
    for(int i=0;i<matrice.size();i++) {
        for(int j=0;j<matrice[i].size();j++) {
            if(matrice[i][j] >= 0 && matrice[i][j] != start) {
                color.push_back(0);
                d.push_back(INT_MAX);
                pi.push_back(-1);
            } else if(matrice[i][j] == start) {
                color.push_back(1);
                d.push_back(0);
                pi.push_back(-1);
            }
        }
    }
    F.push_back(start);

    while(F.size() > 0) {
        int u = F[0];
        F.pop_front();
        /*cout << "F : ";
        for(int i=0;i<F.size();i++) {
            cout << F[i] << " ";
        }
        cout << endl;
        cout << "u : " << u << endl;*/
        for(int i=0;i<adjList[u].size();i++) {
            int v = adjList[u][i];
            //cout << "Adj : " << v << " , c : " << color[v] << endl;
            if(color[v] == 0) {
                color[v] = 1;
                d[v] = d[u] + 1;
                pi[v] = u;
                F.push_back(v);
            }
        }
        /*cout << "After F : ";
        for(int i=0;i<F.size();i++) {
            cout << F[i] << " ";
        }
        cout << endl;*/
        color[u] = 2;
    }

    /*// Get the best way
    for(int i=0;i<d.size();i++) {
        cout << "S to " << i << " : " << d[i] << endl;
    }

    for(int i=0;i<pi.size();i++) {
        cout << "Predecesor of " << i << " is " << pi[i] << endl;
    }*/

    // Print the answer
    string bestPath = "";
    int best = -1;
    for(int i=0;i<ends.size();i++) {
        int current = ends[i];
        int distance = 0;
        string path = "";
        while(current != start && current != -1) {
            int xCurrent;
            int yCurrent;
            int xPiCurrent;
            int yPiCurrent;
            for(int i=0;i<matrice.size();i++) {
                for(int j=0;j<matrice[i].size();j++) {
                    if(matrice[i][j] == current) {
                        xCurrent = i;
                        yCurrent = j;
                    }
                    if(matrice[i][j] == pi[current]) {
                        xPiCurrent = i;
                        yPiCurrent = j;
                    }
                }
            }
            if(yCurrent-yPiCurrent == 1) {
                path += "R ";
            } else if(yCurrent-yPiCurrent == -1) {
                path += "L ";
            } else if(xCurrent-xPiCurrent < 0) {
                path += "U ";
            }
            else if(xCurrent-xPiCurrent > 0) {
                path += "D ";
            }
            current = pi[current];
            distance++;
        }
        if(current == start && (distance < best || best == -1)) {
            best = distance;
            bestPath = path;
        }
    }
    if(best > 0) {
        bestPath.pop_back();
        for (string::reverse_iterator rit=bestPath.rbegin(); rit!=bestPath.rend(); ++rit)
            cout << *rit;
        cout << endl;
        return 0;
    } else {
        cout << "no path" << endl;
    }

    return 0;
}
#include <iostream>
#include <vector>
#include <cmath>
#include <deque>

using namespace std;

int main() {
    int Ar, Ac, Br, Bc, PIr, PIc, nS;
    vector< vector<double> > A;
    vector< vector<double> > B;
    vector<double> PI;
    vector<int> Sequence;
    // -------- GET A ----------
    cin >> Ar >> Ac;
    for(int i=0;i<Ar;i++) {
        vector<double> row;
        for(int j=0;j<Ac;j++) {
            double n;
            cin >> n;
            row.push_back(n);
        }
        A.push_back(row);
    }
    // -------- GET B ----------
    cin >> Br >> Bc;
    for(int i=0;i<Br;i++) {
        vector<double> row;
        for(int j=0;j<Bc;j++) {
            double n;
            cin >> n;
            row.push_back(n);
        }
        B.push_back(row);
    }
    // -------- GET PI ----------
    cin >> PIr >> PIc;
    for(int j=0;j<PIc;j++) {
        double n;
        cin >> n;
        PI.push_back(n);
    }
    // -------- GET Sequence ----------
    cin >> nS;
    for(int j=0;j<nS;j++) {
        int n;
        cin >> n;
        Sequence.push_back(n);
    }

    // -------- COMPUTE SEQUENCE ----------
    vector<vector<double> > delta;
    vector<vector<int> > deltaCarac;
    double max = 0;
    double ind = -1;
    // Init Delta and DeltaState
    for(int i=0;i<Sequence.size();i++) {
        vector<double> a;
        vector<int> b;
        for(int j=0;j<PI.size();j++) {
            a.push_back(-1);
            b.push_back(-1);
        }
        delta.push_back(a);
        deltaCarac.push_back(b);
    }
    // Compute Delta 0
    for(int i=0;i<PI.size();i++) {
        delta[0][i] = PI[i] * B[i][Sequence[0]];
        if(delta[0][i] > max) {
            max = delta[0][i];
            ind = i;
        }
        //cout << delta[0][i]  << " ";
    }
    //cout << endl << " ---------" << endl;
    // Compute Delta 1..T
    vector<double> mat;
    for(int t=1;t<Sequence.size();t++) {
        mat.clear();
        for(int i=0;i<A.size();i++) {
            max = 0;
            ind = -1;
            for(int j=0;j<A[i].size();j++) {
                //cout << A[j][i] << " * " << delta[j][t-1] << " * " << B[i][Sequence[t]] << " + ";
                double c = A[j][i] * delta[t-1][j] * B[i][Sequence[t]];
                //cout << c << endl;
                if(c > max) {
                    max = c;
                    ind = j;
                }
            }
            //cout << endl;
            delta[t][i] = max;
            deltaCarac[t][i] = ind;
            //cout << delta[t][i] << " ";
            //cout << " (" << deltaState[i][t] << ") " << endl;
        }
        //cout << endl << " ---------" << endl;
    }
    // Backtrack
    deque<int> best;
    // Compute the last term
    int l = delta.size()-1;
    max = 0;
    ind = -1;
    for(int j=0;j<delta[l].size();j++) {
        if(delta[l][j] > max) {
            max = delta[l][j];
            ind = j;
        }
    }
    best.push_back(ind);
    // Compute the others terms by backtring
    int b = ind;
    for(int i=delta.size()-1;i>=1;i--) {
        max = 0;
        ind = -1;
        for(int j=0;j<delta[i].size();j++) {
            //cout << delta[i][j] << " ";
            //cout << delta[i-1][j] << " * " << A[j][b] << " = " << delta[i-1][j] * A[j][b] << endl;
            //cout << (delta[i-1][j] * A[i][j] > max) << endl;
            if(delta[i-1][j] * A[j][b] > max) {
                max = delta[i-1][j] * A[j][b];
                ind = j;
            }
        }
        b = ind;
        //cout << "BEST : " << ind << endl;
        best.push_front(ind);
    }
    // Print answer
    for(int i=0;i<best.size();i++) {
        cout << best[i] << " ";
    }
    cout << endl;
    return 0;
}
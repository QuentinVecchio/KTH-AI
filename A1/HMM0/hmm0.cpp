#include <iostream>
#include <vector>

using namespace std;

int main() {
    int Ar, Ac, Br, Bc, PIr, PIc;
    vector< vector<double> > A;
    vector< vector<double> > B;
    vector<double> PI;
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
    // -------- COMPUTE NEW PI ----------
    vector<double> PI2;
    for(int i=0;i<PI.size();i++) {
        PI2.push_back(0);
    }
    for(int i=0;i<PI2.size();i++) {
        for(int j=0;j<A[i].size();j++) {
             PI2[j] += PI[i] * A[i][j];
        }
    }

    // -------- COMPUTE NEW OBSERVATION ----------
    vector<double> observation;
    for(int i=0;i<Bc;i++) {
        observation.push_back(0);
    }
    for(int i=0;i<Bc;i++) {
        for(int j=0;j<PI2.size();j++) {
            observation[i] +=  PI2[j] * B[j][i];
        }
    }
    /*cout << "----- A -----" << endl;
    for(int i=0;i<A.size();i++) {
        for(int j=0;j<A[i].size();j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    cout << "----- B -----" << endl;
    for(int i=0;i<B.size();i++) {
        for(int j=0;j<B[i].size();j++) {
            cout << B[i][j] << " ";
        }
        cout << endl;
    }

    cout << "----- PI -----" << endl;
    for(int i=0;i<PI.size();i++) {
        cout << PI[i] << " ";
    }
    cout << endl;
    
    cout << "----- PI2 -----" << endl;
    for(int i=0;i<PI2.size();i++) {
        cout << PI2[i] << " ";
    }
    
    cout << "----- Obs -----" << endl;
    for(int i=0;i<observation.size();i++) {
        cout << observation[i] << " ";
    }
    */
    cout << 1 << " " << Bc;
    for(int i=0;i<observation.size();i++) {
        cout << " " << observation[i];
    }
    cout << endl;
    return 0;
}
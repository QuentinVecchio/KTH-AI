#include <iostream>
#include <vector>

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
    // Creation treilli
    vector<vector<double> > alpha;
    for(int i=0;i<PI.size();i++) {
        vector<double> a;
        for(int j=0;j<Sequence.size();j++) {
            a.push_back(-1);
        }
        alpha.push_back(a);
    }
    // Initialisation
    for(int i=0;i<PI.size();i++) {
        alpha[i][0] = PI[i] * B[i][Sequence[0]];
    }
    /*cout << "------ " << endl;
    for(int i=0;i<alpha.size();i++) {
        cout << alpha[i][0] << endl;
    }*/
    // Récursivité
    double p = 0;
    vector<double> mat;
    for(int t=1;t<Sequence.size();t++) {
        p = 0;
        mat.clear();
        for(int i=0;i<A.size();i++) {
            double res = 0;
            for(int j=0;j<A[i].size();j++) {
                res += A[j][i] * alpha[j][t-1];
            }
            mat.push_back(res);
        }
        for(int i=0;i<PI.size();i++) {
            int choice = Sequence[t];
            double m = mat[i];
            alpha[i][t] = m * B[i][choice];
            p += alpha[i][t];
        }
    }
    cout << p << endl;

    // -------- COMPUTE NEW OBSERVATION ----------
    /*vector<double> observation;
    for(int i=0;i<Bc;i++) {
        observation.push_back(0);
    }
    for(int i=0;i<Bc;i++) {
        for(int j=0;j<PI2.size();j++) {
            observation[i] +=  PI2[j] * B[j][i];
        }
    }*/

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

    cout << "----- Sequence -----" << endl;
    for(int i=0;i<Sequence.size();i++) {
        cout << Sequence[i] << " ";
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
    
    cout << 1 << " " << Bc;
    for(int i=0;i<observation.size();i++) {
        cout << " " << observation[i];
    }
    cout << endl;*/
    return 0;
}
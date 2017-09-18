#include <iostream>
#include <vector>
#include <climits>
#include <deque>

using namespace std;

vector<vector<long double> > getMatrice(int n, int m);
vector<long double> getdoubleVector(int n);
vector<int> getIntVector(int n);
void printMatrice(vector<vector<long double> > matrice);
void printdoubleVector(vector<long double> vec);
void printIntVector(vector<int> vec);
vector<vector<long double> > initMatrice(int n, int m, long double value = -1);
vector<int> initIntVector(int n, int value = -1);
vector<long double> initdoubleVector(int n, long double value = -1);
vector<long double> computeTransition(vector<vector<long double> > A, vector<long double> pi);
vector<long double> computeObservation(vector<vector<long double> > B, vector<long double> pi);
vector<vector<long double> > computeAlpha(vector<vector<long double> > A, vector<vector<long double> > B, vector<long double> pi, vector<int> sequence);
vector<vector<long double> > computeDelta(vector<vector<long double> > A, vector<vector<long double> > B, vector<long double> pi, vector<int> sequence);
deque<int> computeBestLikely(vector<vector<long double> > A, vector<vector<long double> > delta);
vector<vector<long double> > computeBeta(vector<vector<long double> > A, vector<vector<long double> > B, vector<long double> pi, vector<int> sequence);
vector<vector<vector<long double> > > computeDiGamma(vector<vector<long double> > A, vector<vector<long double> > B, vector<int> sequence, vector<vector<long double> > alpha, vector<vector<long double> > beta);
long double computeProbability(vector<vector<long double> > alpha);
vector<vector<long double> > transitionEstimates(vector<vector<vector<long double> > > diGamma, vector<vector<long double> > gamma, vector<vector<long double> > alpha, vector<vector<long double> > A);
vector<vector<long double> > emissionEstimates(vector<vector<vector<long double> > > diGamma, vector<vector<long double> > gamma, vector<int> sequence, vector<vector<long double> > B);
vector<vector<long double> > computeGamma(vector<vector<vector<long double> > > diGamma, vector<vector<long double> > alpha);

int main() {
    int Ar, Ac, Br, Bc, PIr, PIc, nS;
    vector< vector<long double> > A, B, alpha, beta, gamma;
    vector<vector<vector<long double> > > diGamma;
    vector<long double> PI;
    vector<int> Sequence;
    // -------- GET A ----------
    cin >> Ar >> Ac;
    A = getMatrice(Ar, Ac);
    // -------- GET B ----------
    cin >> Br >> Bc;
    B = getMatrice(Br, Bc);
    // -------- GET PI ----------
    cin >> PIr >> PIc;
    PI = getdoubleVector(PIc);
    // -------- GET Sequence ----------
    cin >> nS;
    Sequence = getIntVector(nS);
    int maxIteration = 100;
    int i = 0;
    long double oldProbability = 0;
    alpha = computeAlpha(A, B, PI, Sequence);
    long double previousProbability = computeProbability(alpha);
    long double currentProbability = previousProbability;
    long double logProbability;
    //cout << "############ > " << currentProbability << endl;
    do {
        previousProbability = currentProbability;
        //cout << "Compute beta ..." << endl;
        beta = computeBeta(A, B, PI, Sequence);
        //cout << "Compute di-gamma ..." << endl;
        diGamma = computeDiGamma(A, B, Sequence, alpha, beta);
        //cout << "Compute Gamma ..." << endl;
        gamma = computeGamma(diGamma, alpha);
        //cout << "Compute new PI ..." << endl;
        PI = gamma[0];
        //cout << "Compute new A ..." << endl;
        A = transitionEstimates(diGamma, gamma, alpha, A);
        //cout << "Compute new B ..." << endl;
        B = emissionEstimates(diGamma, gamma, Sequence, B);
        //cout << "Compute alpha with new lambda ..." << endl;
        alpha = computeAlpha(A, B, PI, Sequence);
        currentProbability = computeProbability(alpha);
        //cout << i << " ############ > " << currentProbability << endl;
        i++;
    } while(currentProbability > previousProbability && i < maxIteration);
    cout << A.size() << " " << A[0].size() << " ";
    for(int i=0;i<A.size();i++) {
        for(int j=0;j<A[i].size();j++) {
            cout << A[i][j] << " ";
        }
    }
    cout << endl;
    cout << B.size() << " " << B[0].size() << " ";
    for(int i=0;i<B.size();i++) {
        for(int j=0;j<B[i].size();j++) {
            cout << B[i][j] << " ";
        }
    }
    cout << endl;
    return 0;
}

vector<long double> computeTransition(vector<vector<long double> > A, vector<long double> pi) {
    vector<long double> newPi = initdoubleVector(pi.size(), 0);
    for(int i=0;i<newPi.size();i++) {
        for(int j=0;j<A[i].size();j++) {
             newPi[j] += pi[i] * A[i][j];
        }
    }
    return newPi;
}

vector<long double> computeObservation(vector<vector<long double> > B, vector<long double> pi) {
    vector<long double> obs = initdoubleVector(B[0].size(), 0);
    for(int i=0;i<B[0].size();i++) {
        for(int j=0;j<pi.size();j++) {
            obs[i] +=  pi[j] * B[j][i];
        }
    }
    return obs;
}

long double computeProbability(vector<vector<long double> > alpha) {
    long double p = 0;
    for(int i=0;i<alpha[alpha.size()-1].size();i++) {
        p += alpha[alpha.size()-1][i];
    }
    return p;
} 

vector<vector<long double> > transitionEstimates(vector<vector<vector<long double> > > diGamma, vector<vector<long double> > gamma, vector<vector<long double> > alpha, vector<vector<long double> > A) {
    vector<vector<long double> > newA = initMatrice(A.size(), A.size());
    for(int i=0;i<A.size();i++) {
        for(int j=0;j<A[i].size();j++) {
            long double top = 0;
            long double bottom = 0;
            for(int t=0;t<diGamma.size()-1;t++) {
                top += diGamma[t][i][j];
                bottom += gamma[t][i];
            }
            if(bottom > 0)
                newA[i][j] = top/bottom;
        }
    }
    return newA;
}

vector<vector<long double> > emissionEstimates(vector<vector<vector<long double> > > diGamma, vector<vector<long double> > gamma, vector<int> sequence, vector<vector<long double> > B) {
    vector<vector<long double> > newB = initMatrice(B.size(), B[0].size());
    for(int i=0;i<B.size();i++) {
        for(int k=0;k<B[i].size();k++) {
            long double top = 0;
            long double bottom = 0;
            for(int t=0;t<diGamma.size();t++) {
                if(sequence[t] == k) {
                    top += gamma[t][i];
                }
                bottom += gamma[t][i];
            }
            if(bottom > 0)
                newB[i][k] = top/bottom;
        }
    }
    return newB;
}

vector<vector<long double> > computeGamma(vector<vector<vector<long double> > > diGamma, vector<vector<long double> > alpha) {
    vector<vector<long double> > gamma = initMatrice(diGamma.size(), diGamma[0].size(), 0);
    for(int t=0;t<diGamma.size()-1;t++) {
        for(int i=0;i<diGamma[t].size();i++) {
            for(int j=0;j<diGamma[t][i].size();j++) {
                gamma[t][i] += diGamma[t][i][j];
            }
        }
    }
    long double bottom = 0;
    int t = gamma.size()-1;
    for(int i=0;i<diGamma[t].size();i++) {
        bottom += alpha[t][i];
    }
    for(int i=0;i<diGamma[t].size();i++) {
        gamma[t][i] += alpha[t][i]/bottom;
    }
    return gamma;
}

vector<vector<vector<long double> > > computeDiGamma(vector<vector<long double> > A, vector<vector<long double> > B, vector<int> sequence, vector<vector<long double> > alpha, vector<vector<long double> > beta) {
    vector<vector<vector<long double> > > diGamma;
    for(int i=0;i<sequence.size();i++) {
        diGamma.push_back(initMatrice(A.size(), A.size()));
    }
    for(int t=0;t<diGamma.size()-1;t++) {
        long double bottom = 0;
        for(int i=0;i<diGamma[t].size();i++) {
            for(int j=0;j<diGamma[t][i].size();j++) {
                bottom += alpha[t][i] * A[i][j] * B[j][sequence[t+1]] * beta[t+1][j];
            }
        }
        for(int i=0;i<diGamma[t].size();i++) {
            for(int j=0;j<diGamma[t][i].size();j++) {
                long double top = alpha[t][i] * A[i][j] * B[j][sequence[t+1]] * beta[t+1][j];
                if(bottom > 0)
                    diGamma[t][i][j] = top/bottom;
            }
        }
    }
    return diGamma;
}

vector<vector<long double> > computeBeta(vector<vector<long double> > A, vector<vector<long double> > B, vector<long double> pi, vector<int> sequence) {
    vector<vector<long double> > beta = initMatrice(sequence.size(), pi.size());
    for(int i=0;i<beta[beta.size()-1].size();i++) {
        beta[beta.size()-1][i] = 1;
    }
    for(int t = beta.size()-2;t>=0;t--) {
        for(int i=0;i<beta[t].size();i++) {
            beta[t][i] = 0;
            for(int j=0;j<beta[t+1].size();j++) {
                beta[t][i] +=  beta[t+1][j] * A[i][j] * B[j][sequence[t+1]];
            }
        }
    }
    return beta;
}

vector<vector<long double> > computeAlpha(vector<vector<long double> > A, vector<vector<long double> > B, vector<long double> pi, vector<int> sequence) {
    vector<vector<long double> > alpha = initMatrice(sequence.size(), pi.size());
    for(int i=0;i<pi.size();i++) {
        alpha[0][i] = pi[i] * B[i][sequence[0]];
    }
    /*cout << "------ " << endl;
    for(int i=0;i<alpha.size();i++) {
        cout << alpha[0][i] << endl;
    }*/
    // Récursivité
    vector<long double> mat;
    for(int t=1;t<sequence.size();t++) {
        mat.clear();
        for(int i=0;i<A.size();i++) {
            long double res = 0;
            for(int j=0;j<A[i].size();j++) {
                res += A[j][i] * alpha[t-1][j];
            }
            mat.push_back(res);
        }
        for(int i=0;i<pi.size();i++) {
            alpha[t][i] = mat[i] * B[i][sequence[t]];
        }
    }
    return alpha;
}

vector<vector<long double> > computeDelta(vector<vector<long double> > A, vector<vector<long double> > B, vector<long double> pi, vector<int> sequence) {
    vector<vector<long double> > delta = initMatrice(sequence.size(), pi.size());
    long double max = 0;
    int ind = -1;
    // Compute Delta 0
    for(int i=0;i<pi.size();i++) {
        delta[0][i] = pi[i] * B[i][sequence[0]];
        if(delta[0][i] > max) {
            max = delta[0][i];
            ind = i;
        }
        //cout << delta[0][i]  << " ";
    }
    //cout << endl << " ---------" << endl;
    // Compute Delta 1..T
    vector<long double> mat;
    for(int t=1;t<sequence.size();t++) {
        mat.clear();
        for(int i=0;i<A.size();i++) {
            max = 0;
            ind = -1;
            for(int j=0;j<A[i].size();j++) {
                //cout << A[j][i] << " * " << delta[j][t-1] << " * " << B[i][Sequence[t]] << " + ";
                long double c = A[j][i] * delta[t-1][j] * B[i][sequence[t]];
                //cout << c << endl;
                if(c > max) {
                    max = c;
                    ind = j;
                }
            }
            //cout << endl;
            delta[t][i] = max;
            //cout << delta[t][i] << " ";
            //cout << " (" << deltaState[i][t] << ") " << endl;
        }
        //cout << endl << " ---------" << endl;
    }
    return delta;
}

deque<int> computeBestLikely(vector<vector<long double> > A, vector<vector<long double> > delta) {
    deque<int> best;
    // Compute the last term
    int l = delta.size()-1;
    long double max = 0;
    int ind = -1;
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
    return best;
}

vector<vector<long double> > getMatrice(int n, int m) {
    vector<vector<long double> > matrice;
    for(int i=0;i<n;i++) {
        vector<long double> row;
        for(int j=0;j<m;j++) {
            long double nb;
            cin >> nb;
            row.push_back(nb);
        }
        matrice.push_back(row);
    }
    return matrice;
}

vector<vector<long double> > initMatrice(int n, int m, long double value) {
    vector<vector<long double> > matrice;
    for(int i=0;i<n;i++) {
        vector<long double> row;
        for(int j=0;j<m;j++) {
            row.push_back(value);
        }
        matrice.push_back(row);
    }
    return matrice;
}

vector<int> initIntVector(int n, int value) {
    vector<int> vec;
    for(int i=0;i<n;i++) {
        vec.push_back(value);
    }
    return vec;
}

vector<long double> initdoubleVector(int n, long double value) {
    vector<long double> vec;
    for(int i=0;i<n;i++) {
        vec.push_back(value);
    }
    return vec;
}

vector<long double> getdoubleVector(int n) {
    vector<long double> vect;
    for(int j=0;j<n;j++) {
        long double nb;
        cin >> nb;
        vect.push_back(nb);
    }
    return vect;
}

vector<int> getIntVector(int n) {
    vector<int> vect;
    for(int j=0;j<n;j++) {
        int nb;
        cin >> nb;
        vect.push_back(nb);
    }
    return vect;
}

void printMatrice(vector<vector<long double> > matrice) {
    for(int i=0;i<matrice.size();i++) {
        for(int j=0;j<matrice[i].size();j++) {
            cout << matrice[i][j] << " ";
        }
        cout << endl;
    }
}

void printdoubleVector(vector<long double> vec) {
    for(int i=0;i<vec.size();i++) {
        cout << vec[i] << " ";
    }
    cout << endl;
}

void printIntVector(vector<int> vec) {
    for(int i=0;i<vec.size();i++) {
        cout << vec[i] << " ";
    }
    cout << endl;
}
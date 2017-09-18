#include <iostream>
#include <cmath>

using namespace std;

double** getMatrice(const unsigned n, const unsigned m);
double* getDoubleVector(const unsigned n);
int* getIntVector(const unsigned n);
void printMatrice(double** matrice, const unsigned n, const unsigned m);
void printDoubleVector(double* vect, const unsigned n);
void printIntVector(int* vect, const unsigned n);
double** initMatrice(const unsigned n, const unsigned m, double value = -1);
int* initIntVector(const unsigned n, int value = -1);
double* initDoubleVector(const unsigned n, double value = -1);
double computeAlpha(double** alpha, double* C, double** A, double** B, double* pi, int* sequence, const unsigned nbSuiteObservation, const unsigned nbState);
double computeProbability(double** alpha, int n, int m);
void computeDelta(double** delta, double** A, double** B, double* pi, int* sequence, const unsigned An, const unsigned Am, const unsigned DeltaN, const unsigned DeltaM);
int* computeBestLikely(double** A, double** delta, int Sn, const unsigned An, const unsigned Am, const unsigned Dn, const unsigned Dm);
void computeBeta(double** beta, double* C, double** A, double** B, int* sequence, const unsigned BetaN, const unsigned BetaM);
double* computeObservation(double** B, double* pi, const unsigned nbState, const unsigned nbObservation);
double* computeTransition(double** A, double* pi, const unsigned nbState);
double computeLog(double * C, const unsigned Sn);
void transitionAndEmissionEstimates(double** A, double** B, double*** diGamma, double** gamma, double** alpha, const unsigned nbState, const unsigned diGammaT, int* sequence, const unsigned nbObservation);
void computeDiGammaAndGamma(double*** diGamma, double** gamma, double** A, double** B, int* sequence, double** alpha, double** beta, const unsigned nbSuiteObservation, const unsigned nbState);

int main() {
    std::ios::sync_with_stdio(false);
    unsigned nbState, nbObservation, nbSuiteObservation, fuck;
    double **A, **B, **alpha, **beta, **gamma, **delta;
    double *** diGamma;
    double *PI, *C;
    int *Sequence;
    // -------- GET A ----------
    cin >> nbState >> fuck;
    A = getMatrice(nbState, nbState);
    // -------- GET B ----------
    cin >> nbState >> nbObservation;
    B = getMatrice(nbState, nbObservation);
    // -------- GET PI ----------
    cin >> fuck >> nbState;
    PI = getDoubleVector(nbState);
    // -------- GET Sequence ----------
    cin >> nbSuiteObservation;
    Sequence = getIntVector(nbSuiteObservation);
    // -------- Compute sizes ----------
    const unsigned maxIteration = 30;
    double oldLogProb = -1000000, logProb = 0;
    unsigned i = 0, j;
    // -------- Init all matrices ----------
    alpha = initMatrice(nbSuiteObservation, nbState);
    beta = initMatrice(nbSuiteObservation, nbState);
    diGamma = new double**[nbSuiteObservation];
    for(i=0;i<nbSuiteObservation;++i) {
        diGamma[i] = initMatrice(nbState, nbState);
    }
    i = 0;
    gamma = initMatrice(nbSuiteObservation, nbState, 0);
    C = initDoubleVector(nbSuiteObservation, 0);
    // -------- Learning  ----------
    do {
        if(i > 0)
            oldLogProb = logProb;
        logProb = computeAlpha(alpha, C, A, B, PI, Sequence, nbSuiteObservation, nbState);
        if(!((logProb - oldLogProb) > 0.01))
            break;
        computeBeta(beta, C, A, B, Sequence, nbSuiteObservation, nbState);
        //printMatrice(beta, nbSuiteObservation, nbState);
        
        computeDiGammaAndGamma(diGamma, gamma, A, B, Sequence, alpha, beta, nbSuiteObservation, nbState);
        PI = gamma[0];
        transitionAndEmissionEstimates(A, B, diGamma, gamma, alpha, nbState, nbSuiteObservation, Sequence, nbObservation);
        ++i;
        //cout << i << " " << logProb << " " << (logProb - oldLogProb) << endl;
    } while(i < maxIteration);
    cout << nbState << " " << nbState << " ";
    for(i=0;i<nbState;++i) {
        for(j=0;j<nbState;++j) {
            cout << A[i][j] << " ";
        }
    }
    cout << endl;
    cout << nbState << " " << nbObservation << " ";
    for(i=0;i<nbState;++i) {
        for(j=0;j<nbObservation;++j) {
            cout << B[i][j] << " ";
        }
    }
    cout << endl;
    return 0;
}

/*double computeLog(double * C, const unsigned Sn) {
    double logProb = 0;
    for(unsigned a=0;a<Sn;a++) {
        logProb += log(C[a]);
    }
    return -logProb;
}*/

double* computeTransition(double** A, double* pi, const unsigned An, const unsigned Am) {
    double* newPi = initDoubleVector(An, 0);
    for(unsigned i=0;i<An;++i) {
        for(unsigned j=0;j<Am;++j) {
             newPi[j] += pi[i] * A[i][j];
        }
    }
    return newPi;
}

double* computeObservation(double** B, double* pi, const unsigned Bn, const unsigned Bm) {
    double* obs = initDoubleVector(Bm, 0);
    for(unsigned i=0;i<Bm;++i) {
        for(unsigned j=0;j<Bn;++j) {
            obs[i] +=  pi[j] * B[j][i];
        }
    }
    return obs;
}

double computeProbability(double** alpha, const unsigned n, const unsigned m) {
    double p = 0;
    for(unsigned i=0;i<m;++i) {
        p += alpha[n-1][i];
    }
    return p;
} 

void transitionAndEmissionEstimates(double** A, double** B, double*** diGamma, double** gamma, double** alpha, const unsigned nbState, const unsigned nbSuiteObservation, int* sequence, const unsigned nbObservation) {
    double top = 0;
    double top2 = 0;
    double bottom = 0;
    for(unsigned i=0;i<nbState;++i) {
        for(unsigned j=0;j<nbState||j<nbObservation;++j) {
            top = top2 = bottom = 0;
            for(unsigned t=0;t<nbSuiteObservation-1;++t) {
                if(j<nbState)
                    top += diGamma[t][i][j];
                bottom += gamma[t][i];
                if(sequence[t] == j) {
                    top2 += gamma[t][i];
                }
            }
            if(bottom > 0) {
                if(j<nbObservation)
                    B[i][j] = top2/bottom;
                if(j<nbState)
                    A[i][j] = top/bottom;
            } 
        }
    }
}

void computeDiGammaAndGamma(double*** diGamma, double** gamma, double** A, double** B, int* sequence, double** alpha, double** beta, const unsigned nbSuiteObservation, const unsigned nbState) {
    double bottom = 0;
    double top;
    unsigned i, j;
    for(unsigned t=0;t<nbSuiteObservation-1;++t) {
        bottom = 0;
        for(i=0;i<nbState;++i) {
            for(j=0;j<nbState;++j) {
                bottom += alpha[t][i] * A[i][j] * B[j][sequence[t+1]] * beta[t+1][j];
            }
        }
        for(i=0;i<nbState;++i) {
            gamma[t][i] = 0;
            for(j=0;j<nbState;++j) {
                top = alpha[t][i] * A[i][j] * B[j][sequence[t+1]] * beta[t+1][j];
                if(bottom > 0)
                    diGamma[t][i][j] = top/bottom;
                gamma[t][i] += diGamma[t][i][j];
            }
        }
    }
    bottom = 0;
    for(i=0;i<nbState;++i) {
        bottom += alpha[nbSuiteObservation-1][i];
    }
    for(i=0;i<nbState;++i) {
        gamma[nbSuiteObservation-1][i] = alpha[nbSuiteObservation-1][i]/bottom;
    }
}

void computeBeta(double** beta, double* C, double** A, double** B, int* sequence, const unsigned nbSuiteObservation, const unsigned nbState) {
    unsigned i;
    for(i=0;i<nbState;++i) {
        beta[nbSuiteObservation-1][i] = C[nbSuiteObservation-1];
    }
    for(int t = nbSuiteObservation-2;t>=0;t--) {
        for(i=0;i<nbState;++i) {
            beta[t][i] = 0;
            for(unsigned j=0;j<nbState;++j) {
                beta[t][i] +=  beta[t+1][j] * A[i][j] * B[j][sequence[t+1]];
            }
            beta[t][i] *= C[t];
        }
    }
}

double computeAlpha(double** alpha, double* C, double** A, double** B, double* pi, int* sequence, const unsigned nbSuiteObservation, const unsigned nbState) {    
    C[0] = 0;
    unsigned i;
    for(i=0;i<nbState;++i) {
        alpha[0][i] = pi[i] * B[i][sequence[0]];
        C[0] += alpha[0][i];
    }
    C[0] = 1 / C[0];
    double logProb = log(C[0]);
    for(i=0;i<nbState;++i) {
        alpha[0][i] *= C[0];
    }
    // Récursivité
    for(unsigned t=1;t<nbSuiteObservation;++t) {
        C[t] = 0;
        for(i=0;i<nbState;++i) {
            alpha[t][i] = 0;
            for(unsigned j=0;j<nbState;++j) {
                alpha[t][i] += A[j][i] * alpha[t-1][j];
            }
        }
        for(i=0;i<nbState;++i) {
            alpha[t][i] *= B[i][sequence[t]];
            C[t] += alpha[t][i];
        }
        C[t] = 1/C[t];
        logProb += log(C[t]);
        for(i=0;i<nbState;++i) {
            alpha[t][i] *= C[t];
        }
    }
    return -logProb;
}

void computeDelta(double** delta, double** A, double** B, double* pi, int* sequence, const unsigned An, const unsigned Am, const unsigned DeltaN, const unsigned DeltaM) {
    double max = 0;
    double c;
    int ind = -1;
    unsigned i;
    for(i=0;i<DeltaM;++i) {
        delta[0][i] = pi[i] * B[i][sequence[0]];
        if(delta[0][i] > max) {
            max = delta[0][i];
            ind = i;
        }
    }
    for(unsigned t=1;t<DeltaN;++t) {
        for(i=0;i<An;++i) {
            max = 0;
            ind = -1;
            for(unsigned j=0;j<Am;++j) {
                c = A[j][i] * delta[t-1][j] * B[i][sequence[t]];
                if(c > max) {
                    max = c;
                    ind = j;
                }
            }
            delta[t][i] = max;
        }
    }
}

int* computeBestLikely(double** A, double** delta, const unsigned Sn, const unsigned An, const unsigned Am, const unsigned Dn, const unsigned Dm) {
    int *best = new int[Sn];
    int l = Dn - 1;
    double max = 0;
    int ind = -1;
    for(unsigned j=0;j<Dm;++j) {
        if(delta[l][j] > max) {
            max = delta[l][j];
            ind = j;
        }
    }
    int index = Sn-1;
    best[index] = ind;
    unsigned b = ind;
    for(unsigned i=Dn-1;i>=1;i--) {
        index--;
        max = 0;
        ind = -1;
        for(unsigned j=0;j<Dm;++j) {
            if(delta[i-1][j] * A[j][b] > max) {
                max = delta[i-1][j] * A[j][b];
                ind = j;
            }
        }
        b = ind;
        best[index] = ind;
    }
    return best;
}

double** getMatrice(const unsigned n, const unsigned m) {
    double** matrice = new double*[n];
    for(unsigned i=0;i<n;++i) {
        matrice[i] = new double[m];
        for(unsigned j=0;j<m;++j) {
            cin >> matrice[i][j];
        }
    }
    return matrice;
}

double** initMatrice(const unsigned n, const unsigned m, double value) {
    double** matrice = new double*[n];
    for(unsigned i=0;i<n;++i) {
        matrice[i] = new double[m];
        for(unsigned j=0;j<m;++j) {
            matrice[i][j] = value;
        }
    }
    return matrice;
}

int* initIntVector(const unsigned n, int value) {
    int* vect = new int[n];
    for(unsigned i=0;i<n;++i) {
        vect[i] = value;
    }
    return vect;
}

double* initDoubleVector(const unsigned n, double value) {
    double *vect = new double[n];
    for(unsigned i=0;i<n;++i) {
        vect[i] = value;
    }
    return vect;
}

double* getDoubleVector(const unsigned n) {
    double * vect = new double[n];
    double nb;
    for(unsigned j=0;j<n;++j) {
        cin >> nb;
        vect[j] = nb;
    }
    return vect;
}

int* getIntVector(const unsigned n) {
    int * vect = new int[n];
    unsigned nb;
    for(unsigned j=0;j<n;++j) {
        cin >> nb;
        vect[j] = nb;
    }
    return vect;
}

void printMatrice(double ** matrice, const unsigned N, const unsigned M) {
    for(unsigned i=0;i<N;++i) {
        for(unsigned j=0;j<M;++j) {
            cout << matrice[i][j] << " ";
        }
        cout << endl;
    }
}

void printDoubleVector(double *  vect, const unsigned size) {
    for(unsigned i=0;i<size;++i) {
        cout << vect[i] << " ";
    }
    cout << endl;
}

void printIntVector(int* vect, const unsigned size) {
    for(unsigned i=0;i<size;++i) {
        cout << vect[i] << " ";
    }
    cout << endl;
}
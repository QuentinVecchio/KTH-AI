#include <iostream>
#include <cmath>

using namespace std;

double* getMatrice(const unsigned n, const unsigned m);
double* getDoubleVector(const unsigned n);
unsigned char* getIntVector(const unsigned n);
void printMatrice(double** matrice, const unsigned n, const unsigned m);
void printDoubleVector(double* vect, const unsigned n);
void printIntVector(int* vect, const unsigned n);
double* initMatrice(const unsigned n, const unsigned m, double value = -1);
unsigned char* initIntVector(const unsigned n, int value = -1);
double* initDoubleVector(const unsigned n, double value = -1);
double computeAlpha();
void computeBeta();
void transitionAndEmissionEstimates();
void computeDiGammaAndGamma();

unsigned nbState, nbObservation, nbSuiteObservation, fuck;
double *A, *B, *alpha, *beta, *Gamma, **diGamma;
double *PI, *C;
unsigned char  *Sequence;

int main() {
    std::ios::sync_with_stdio(false);
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
    diGamma = new double*[nbSuiteObservation];//(double*)malloc(nbSuiteObservation*nbState*nbState*sizeof(double));
    for(unsigned i=0;i<nbSuiteObservation;i++) {
        diGamma[i] = initMatrice(nbState, nbState, 0);
    }
    i = 0;
    Gamma = initMatrice(nbSuiteObservation, nbState, 0);
    C = initDoubleVector(nbSuiteObservation, 0);
    // -------- Learning  ----------
    do {
        if(i > 0)
            oldLogProb = logProb;
        logProb = computeAlpha();
        if(!((logProb - oldLogProb) > 0.010))
            break;
        computeBeta();
        computeDiGammaAndGamma();
        for(unsigned i=0;i<nbState;i++)
            PI[i] = Gamma[i];
        transitionAndEmissionEstimates();
        ++i;
        //cout << i << " " << logProb << " " << (logProb - oldLogProb) << endl;
    } while(i < maxIteration);
    printf("%u %u ", nbState, nbState);
    for(i=0;i<nbState;++i) {
        for(j=0;j<nbState;++j) {
            printf("%lf ", A[i*nbState+j]);
        }
    }
    printf("\n%u %u ", nbState, nbObservation);
    for(i=0;i<nbState;++i) {
        for(j=0;j<nbObservation;++j) {
            printf("%lf ", B[i*nbObservation+j]);
        }
    }
    printf("\n");
    return 0;
}

void transitionAndEmissionEstimates() {
    double top = 0;
    double top2 = 0;
    double bottom = 0;
    for(unsigned i=0;i<nbState;++i) {
        for(unsigned j=0;j<nbState||j<nbObservation;++j) {
            top = top2 = bottom = 0;
            for(unsigned t=0;t<nbSuiteObservation-1;++t) {
                if(j<nbState)
                    top += diGamma[t][i*nbState+j];
                bottom += Gamma[t*nbState+i];
                if(Sequence[t] == j) {
                    top2 += Gamma[t*nbState+i];
                }
            }
            if(bottom > 0) {
                if(j<nbObservation)
                    B[i*nbObservation+j] = top2/bottom;
                if(j<nbState)
                    A[i*nbState+j] = top/bottom;
            } 
        }
    }
}

void computeDiGammaAndGamma() {
    double bottom = 0;
    double top;
    unsigned i, j;
    for(unsigned t=0;t<nbSuiteObservation-1;++t) {
        bottom = 0;
        for(i=0;i<nbState;++i) {
            for(j=0;j<nbState;++j) {
                bottom += alpha[t*nbState+i] * A[i*nbState+j] * B[j*nbObservation+Sequence[t+1]] * beta[(t+1)*nbState+j];
            }
        }
        for(i=0;i<nbState;++i) {
            Gamma[t*nbState+i] = 0;
            for(j=0;j<nbState;++j) {
                top = alpha[t*nbState+i] * A[i*nbState+j] * B[j*nbObservation+Sequence[t+1]] * beta[(t+1)*nbState+j];
                if(bottom > 0)
                    diGamma[t][i*nbState+j] = top/bottom;
                Gamma[t*nbState+i] += diGamma[t][i*nbState+j];
            }
        }
    }
    bottom = 0;
    for(i=0;i<nbState;++i) {
        bottom += alpha[(nbSuiteObservation-1)*nbState+i];
    }
    for(i=0;i<nbState;++i) {
        Gamma[(nbSuiteObservation-1)*nbState+i] = alpha[(nbSuiteObservation-1)*nbState+i]/bottom;
    }
}

void computeBeta() {
    unsigned i;
    for(i=0;i<nbState;++i) {
        beta[(nbSuiteObservation-1)*nbState+i] = C[nbSuiteObservation-1];
    }
    for(int t = nbSuiteObservation-2;t>=0;t--) {
        for(i=0;i<nbState;++i) {
            beta[t*nbState+i] = 0;
            for(unsigned j=0;j<nbState;++j) {
                beta[t*nbState+i] +=  beta[(t+1)*nbState+j] * A[i*nbState+j] * B[j*nbObservation+Sequence[t+1]];
            }
            beta[t*nbState+i] *= C[t];
        }
    }
}

double computeAlpha() {    
    C[0] = 0;
    unsigned i;
    for(i=0;i<nbState;++i) {
        alpha[i] = PI[i] * B[i*nbObservation+Sequence[0]];
        C[0] += alpha[i];
    }
    C[0] = 1 / C[0];
    double logProb = log(C[0]);
    for(i=0;i<nbState;++i) {
        alpha[i] *= C[0];
    }
    // Récursivité
    for(unsigned t=1;t<nbSuiteObservation;++t) {
        C[t] = 0;
        for(i=0;i<nbState;++i) {
            alpha[t*nbState+i] = 0;
            for(unsigned j=0;j<nbState;++j) {
                alpha[t*nbState+i] += A[j*nbState+i] * alpha[(t-1)*nbState+j];
            }
        }
        for(i=0;i<nbState;++i) {
            alpha[t*nbState+i] *= B[i*nbObservation+Sequence[t]];
            C[t] += alpha[t*nbState+i];
        }
        C[t] = 1/C[t];
        logProb += log(C[t]);
        for(i=0;i<nbState;++i) {
            alpha[t*nbState+i] *= C[t];
        }
    }
    return -logProb;
}

double* getMatrice(const unsigned n, const unsigned m) {
    double* matrice = (double*)malloc(n*m*sizeof(double));
    for(unsigned i=0;i<n;++i) {
        for(unsigned j=0;j<m;++j) {
            cin >> matrice[i*m+j];
        }
    }
    return matrice;
}

double* initMatrice(const unsigned n, const unsigned m, double value) {
    double* matrice = (double*)malloc(n*m*sizeof(double));
    for(unsigned i=0;i<n;++i) {
        for(unsigned j=0;j<m;++j) {
            matrice[i*m+j] = value;
        }
    }
    return matrice;
}

unsigned char* initIntVector(const unsigned n, unsigned value) {
    unsigned char* vect = new unsigned char[n];
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

unsigned char* getIntVector(const unsigned n) {
    unsigned char * vect = new unsigned char[n];
    for(unsigned j=0;j<n;++j) {
        cin >> vect[j];
        vect[j] -= '0';
    }
    return vect;
}

void printMatrice(double ** matrice, const unsigned N, const unsigned M) {
    for(unsigned i=0;i<N;++i) {
        for(unsigned j=0;j<M;++j) {
            cout << matrice[i*M+j] << " ";
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

void printIntVector(char* vect, const unsigned size) {
    for(unsigned i=0;i<size;++i) {
        cout << vect[i] << " ";
    }
    cout << endl;
}
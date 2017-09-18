#include <iostream>
#include <climits>
#include <cmath>
#include <deque>
#include <thread>

using namespace std;

long double** getMatrice(int n, int m);
long double* getDoubleVector(int n);
int* getIntVector(int n);
void printMatrice(long double** matrice, int n, int m);
void printDoubleVector(long double* vect, int n);
void printIntVector(int* vect, int n);
long double** initMatrice(int n, int m, long double value = -1);
int* initIntVector(int n, int value = -1);
long double* initDoubleVector(int n, long double value = -1);
void computeAlpha(long double** alpha, long double* C, long double** A, long double** B, long double* pi, int* sequence, int n, int m);
long double computeProbability(long double** alpha, int n, int m);
void computeDelta(long double** delta, long double** A, long double** B, long double* pi, int* sequence, int An, int Am, int DeltaN, int DeltaM);
int* computeBestLikely(long double** A, long double** delta, int Sn, int An, int Am, int Dn, int Dm);
void computeBeta(long double** beta, long double* C, long double** A, long double** B, long double* pi, int* sequence, int BetaN, int BetaM);
void computeDiGamma(long double*** diGamma, long double** A, long double** B, int* sequence, long double** alpha, long double** beta, int diGammaT, int diGammaN, int diGammaM);
void computeGamma(long double** gamma, long double*** diGamma, long double** alpha, int gammaT, int gammaN, int gammaM);
void transitionEstimates(long double** A, long double*** diGamma, long double** gamma, long double** alpha, int An, int Am, int diGammaT, int diGammaN, int diGammaM);
void emissionEstimates(long double** B, long double*** diGamma, long double** gamma, int* sequence, int diGammaT, int Bn, int BM);
long double* computeObservation(long double** B, long double* pi, int Bn, int Bm);
long double* computeTransition(long double** A, long double* pi, int An, int Am);
float computeLog(long double * C, int Sn);

int main() {
    int An, Am, Bn, Bm, PIn, PIm, Sn;
    long double **A, **B, **alpha, **beta, **gamma, **delta;
    long double *** diGamma;
    long double *PI, *C;
    int *Sequence;
    // -------- GET A ----------
    cin >> An >> Am;
    A = getMatrice(An, Am);
    // -------- GET B ----------
    cin >> Bn >> Bm;
    B = getMatrice(Bn, Bm);
    // -------- GET PI ----------
    cin >> PIn >> PIm;
    PI = getDoubleVector(PIm);
    // -------- GET Sequence ----------
    cin >> Sn;
    Sequence = getIntVector(Sn);
    // -------- Compute sizes ----------
    const int AlphaN = Sn, AlphaM = PIm; // Alpha Size
    const int BetaN = Sn, BetaM = PIm; // Beta Size;
    const int diGammaT = Sn, diGammaN = An, diGammaM = Am; // DiGamma Size;
    const int gammaN = diGammaT, gammaM = diGammaM; // DiGamma Size;
    const int maxIteration = 50;
    float oldLogProb = -1000000, logProb = 0;
    int index = 0, j, i;
    // -------- Init all matrices ----------
    alpha = initMatrice(AlphaN, AlphaM);
    beta = initMatrice(BetaN, BetaM);
    diGamma = new long double**[diGammaT];
    for(i=0;i<diGammaT;i++) {
        diGamma[i] = initMatrice(diGammaN, diGammaM);
    }
    index = 0;
    gamma = initMatrice(diGammaT, diGammaN, 0);
    C = initDoubleVector(Sn, 0);
    // -------- Learning  ----------
    do {
        if(index > 0)
            oldLogProb = logProb;
        computeAlpha(alpha, C, A, B, PI, Sequence, AlphaN, AlphaM);
		cout << "ALPHA" << endl;
		/*for(j=0;j<AlphaM;j++) {
			cout << alpha[AlphaN-1][j] << " ";
		}
		cout << endl;*/
        computeBeta(beta, C, A, B, PI, Sequence, BetaN, BetaM);
		/*cout << "BETA" << endl;
		for(j=0;j<BetaM;j++) {
			cout << beta[0][j] << " ";
		}
		cout << endl;*/
        computeDiGamma(diGamma, A, B, Sequence, alpha, beta, diGammaT, diGammaN, diGammaM);
        
		/*cout << "DiGamma" << endl;
		cout << diGamma[0][diGammaN-1][diGammaN-1] << endl;
		cout << diGamma[1][1][1] << endl;
		cout << diGamma[diGammaT-2][diGammaN-1][diGammaN*-1] << endl;
		cout << diGamma[diGammaT-1][diGammaN-1][diGammaN*-1] << endl;
		cout << endl;*/
		
		computeGamma(gamma, diGamma, alpha, diGammaT, diGammaN, diGammaM);
        PI = gamma[0];
        //
        transitionEstimates(A, diGamma, gamma, alpha, An, Am, diGammaT, diGammaN, diGammaM);
        //
        emissionEstimates(B, diGamma, gamma, Sequence, diGammaT, Bn, Bm);
        //
        logProb = computeLog(C, Sn);
        index++;
        //cout << i << " " << logProb << " " << (logProb - oldLogProb) << endl;
    } while((logProb - oldLogProb) > 0.015 && index < maxIteration);
    cout << An << " " << Am << " ";
    for(i=0;i<An;i++) {
        for(j=0;j<Am;j++) {
            cout << A[i][j] << " ";
        }
    }
    cout << endl;
    cout << Bn << " " << Bm << " ";
    for(i=0;i<Bn;i++) {
        for(j=0;j<Bm;j++) {
            cout << B[i][j] << " ";
        }
    }
    cout << endl;
    return 0;
}

float computeLog(long double * C, int Sn) {
    float logProb = 0;
    for(int a=0;a<Sn;a++) {
        logProb += log(C[a]);
    }
    return -logProb;
}

long double* computeTransition(long double** A, long double* pi, int An, int Am) {
    long double* newPi = initDoubleVector(An, 0);
    for(int i=0;i<An;i++) {
        for(int j=0;j<Am;j++) {
             newPi[j] += pi[i] * A[i][j];
        }
    }
    return newPi;
}

long double* computeObservation(long double** B, long double* pi, int Bn, int Bm) {
    long double* obs = initDoubleVector(Bm, 0);
    for(int i=0;i<Bm;i++) {
        for(int j=0;j<Bn;j++) {
            obs[i] +=  pi[j] * B[j][i];
        }
    }
    return obs;
}

long double computeProbability(long double** alpha, int n, int m) {
    long double p = 0;
    for(int i=0;i<m;i++) {
        p += alpha[n-1][i];
    }
    return p;
} 

void transitionEstimates(long double** A, long double*** diGamma, long double** gamma, long double** alpha, int An, int Am, int diGammaT, int diGammaN, int diGammaM) {
    long double top = 0;
    long double bottom = 0;
    for(int i=0;i<An;i++) {
        for(int j=0;j<Am;j++) {
            top = bottom = 0;
            for(int t=0;t<diGammaT-1;t++) {
                top += diGamma[t][i][j];
                bottom += gamma[t][i];
            }
            if(bottom > 0)
                A[i][j] = top/bottom;
        }
    }
}

void emissionEstimates(long double** B, long double*** diGamma, long double** gamma, int* sequence, int diGammaT, int Bn, int Bm) {
    long double top = 0;
    long double bottom = 0;
    for(int i=0;i<Bn;i++) {
        for(int k=0;k<Bm;k++) {
            top = bottom = 0;
            for(int t=0;t<diGammaT;t++) {
                if(sequence[t] == k) {
                    top += gamma[t][i];
                }
                bottom += gamma[t][i];
            }
            if(bottom > 0)
                B[i][k] = top/bottom;
        }
    }
}

void computeGamma(long double** gamma, long double*** diGamma, long double** alpha, int diGammaT, int diGammaN, int diGammaM) {
    int i;
    for(int t=0;t<diGammaT-1;t++) {
        for(i=0;i<diGammaN;i++) {
            gamma[t][i] = 0;
            for(int j=0;j<diGammaM;j++) {
                gamma[t][i] += diGamma[t][i][j];
            }
        }
    }
    long double bottom = 0;
    for(i=0;i<diGammaN;i++) {
        bottom += alpha[diGammaT-1][i];
    }
    for(i=0;i<diGammaN;i++) {
        gamma[diGammaT-1][i] = alpha[diGammaT-1][i]/bottom;
    }
}

void computeDiGamma(long double*** diGamma, long double** A, long double** B, int* sequence, long double** alpha, long double** beta, int diGammaT, int diGammaN, int diGammaM) {
    long double bottom = 0;
    long double top;
    int i, j;
    for(int t=0;t<diGammaT-1;t++) {
        bottom = 0;
        for(i=0;i<diGammaN;i++) {
            for(j=0;j<diGammaM;j++) {
                bottom += alpha[t][i] * A[i][j] * B[j][sequence[t+1]] * beta[t+1][j];
            }
        }
        for(i=0;i<diGammaN;i++) {
            for(j=0;j<diGammaM;j++) {
                top = alpha[t][i] * A[i][j] * B[j][sequence[t+1]] * beta[t+1][j];
                if(bottom > 0)
                    diGamma[t][i][j] = top/bottom;
            }
        }
    }
}

void computeBeta(long double** beta, long double* C, long double** A, long double** B, long double* pi, int* sequence, int BetaN, int BetaM) {
    int i;
    for(i=0;i<BetaM;i++) {
        beta[BetaN-1][i] = C[BetaN-1];
    }
    for(int t = BetaN-2;t>=0;t--) {
        for(i=0;i<BetaM;i++) {
            beta[t][i] = 0;
            for(int j=0;j<BetaM;j++) {
                beta[t][i] +=  beta[t+1][j] * A[i][j] * B[j][sequence[t+1]];
            }
            beta[t][i] *= C[t];
        }
    }
}

void computeAlpha(long double** alpha, long double* C, long double** A, long double** B, long double* pi, int* sequence, int n, int m) {
    C[0] = 0;
    int i;
    for(i=0;i<m;i++) {
        alpha[0][i] = pi[i] * B[i][sequence[0]];
        C[0] += alpha[0][i];
    }
    C[0] = 1 / C[0];
    for(i=0;i<m;i++) {
        alpha[0][i] *= C[0];
    }
    // Récursivité
    for(int t=1;t<n;t++) {
        C[t] = 0;
        for(i=0;i<m;i++) {
            alpha[t][i] = 0;
            for(int j=0;j<m;j++) {
                alpha[t][i] += A[j][i] * alpha[t-1][j];
            }
        }
        for(i=0;i<m;i++) {
            alpha[t][i] *= B[i][sequence[t]];
            C[t] += alpha[t][i];
        }
        C[t] = 1/C[t];
        for(i=0;i<m;i++) {
            alpha[t][i] *= C[t];
        }
    }
}

void computeDelta(long double** delta, long double** A, long double** B, long double* pi, int* sequence, int An, int Am, int DeltaN, int DeltaM) {
    long double max = 0;
    long double c;
    int ind = -1;
    int i;
    for(i=0;i<DeltaM;i++) {
        delta[0][i] = pi[i] * B[i][sequence[0]];
        if(delta[0][i] > max) {
            max = delta[0][i];
            ind = i;
        }
    }
    for(int t=1;t<DeltaN;t++) {
        for(i=0;i<An;i++) {
            max = 0;
            ind = -1;
            for(int j=0;j<Am;j++) {
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

int* computeBestLikely(long double** A, long double** delta, int Sn, int An, int Am, int Dn, int Dm) {
    int *best = new int[Sn];
    int l = Dn - 1;
    long double max = 0;
    int ind = -1;
    for(int j=0;j<Dm;j++) {
        if(delta[l][j] > max) {
            max = delta[l][j];
            ind = j;
        }
    }
    int index = Sn-1;
    best[index] = ind;
    int b = ind;
    for(int i=Dn-1;i>=1;i--) {
        index--;
        max = 0;
        ind = -1;
        for(int j=0;j<Dm;j++) {
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

long double** getMatrice(int n, int m) {
    long double** matrice = new long double*[n];
    long double nb;
    for(int i=0;i<n;i++) {
        matrice[i] = new long double[m];
        for(int j=0;j<m;j++) {
            cin >> nb;
            matrice[i][j] = nb;
        }
    }
    return matrice;
}

long double** initMatrice(int n, int m, long double value) {
    long double** matrice = new long double*[n];
    for(int i=0;i<n;i++) {
        matrice[i] = new long double[m];
        for(int j=0;j<m;j++) {
            matrice[i][j] = value;
        }
    }
    return matrice;
}

int* initIntVector(int n, int value) {
    int* vect = new int[n];
    for(int i=0;i<n;i++) {
        vect[i] = value;
    }
    return vect;
}

long double* initDoubleVector(int n, long double value) {
    long double *vect = new long double[n];
    for(int i=0;i<n;i++) {
        vect[i] = value;
    }
    return vect;
}

long double* getDoubleVector(int n) {
    long double * vect = new long double[n];
    long double nb;
    for(int j=0;j<n;j++) {
        cin >> nb;
        vect[j] = nb;
    }
    return vect;
}

int* getIntVector(int n) {
    int * vect = new int[n];
    int nb;
    for(int j=0;j<n;j++) {
        cin >> nb;
        vect[j] = nb;
    }
    return vect;
}

void printMatrice(double long** matrice, int N, int M) {
    for(int i=0;i<N;i++) {
        for(int j=0;j<M;j++) {
            cout << matrice[i][j] << " ";
        }
        cout << endl;
    }
}

void printDoubleVector(double long*  vect, int size) {
    for(int i=0;i<size;i++) {
        cout << vect[i] << " ";
    }
    cout << endl;
}

void printIntVector(int* vect, int size) {
    for(int i=0;i<size;i++) {
        cout << vect[i] << " ";
    }
    cout << endl;
}
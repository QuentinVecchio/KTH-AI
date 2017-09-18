#include "Hmm.hpp"
#include <ctime>
#include <random>
#include <cstdlib>
#include <cstdio>

using namespace std;

// Publics Methods
HMM::HMM(double** A, double** B, double* PI, const unsigned nbState, const unsigned nbObservation) {
    this->nbState = nbState;
    this->nbObservation = nbObservation;
    this->A = A;
    this->B = B;
    this->PI = PI;
    srand (time(NULL));
}

HMM::HMM(const unsigned nbState, const unsigned nbObservation) {
    this->nbState = nbState;
    this->nbObservation = nbObservation;
    this->A = this->initRandomMatrice(nbState, nbState);
    this->B = this->initRandomMatrice(nbState, nbObservation);
    //this->B = this->initBMatrice(nbState, nbObservation);
    this->PI = this->initRandomVector(nbState);
}

HMM::~HMM() {
    // TODO (delete[])
}

void HMM::learn(int* sequence, const unsigned nbEmission, const unsigned maxIteration) {
    double oldLogProb = -1000000, logProb = 0;
    unsigned i = 0, j;

    // -------- Init all matrices ----------
    double** alpha = this->initMatrice(nbEmission, this->nbState, 0);
    double** beta = this->initMatrice(nbEmission, this->nbState, 0);
    double*** diGamma = new double**[nbEmission];
    for(i=0;i<nbEmission;++i) {
        diGamma[i] = initMatrice(this->nbState, this->nbState, 0);
    }
    double** gamma = initMatrice(nbEmission, this->nbState, 0);
    double* C = initVector(nbEmission, 0);

    // -------- Learning  ----------
    i = 0;
    do {
        if(i > 0)
            oldLogProb = logProb;
        logProb = this->computeAlpha(alpha, sequence, C, nbEmission);
        if(logProb != logProb)        
            break;

        this->computeBeta(beta, C, sequence, nbEmission);
        this->computeDiGammaAndGamma(diGamma, gamma, sequence, alpha, beta, nbEmission);
        PI = gamma[0];
        this->transitionEstimates(diGamma, gamma, alpha, nbEmission);
        this->emissionEstimates(diGamma, gamma, sequence, nbEmission);

        /*cerr << "Alpha ---------" << endl;
        this->printMatrice(alpha, nbEmission, nbState);
        cerr << "Beta ---------" << endl;
        this->printMatrice(beta, nbEmission, nbState);
        cerr << "Gamma ---------" << endl;
        this->printMatrice(gamma, nbEmission, nbState);
        cerr << "A ---------" << endl;
        this->printA();
        cerr << "B ---------" << endl;
        this->printB();
        cerr << "PI ---------" << endl;
        this->printPI();
        cerr << "---------" << endl;*/
        ++i;
    } while(logProb > oldLogProb && i < maxIteration);
    //cerr << "Iteration " << i << endl;
}

double HMM::emissionProbability(int *sequence, const unsigned nbEmission) {
    double** alpha = this->initMatrice(nbEmission, this->nbState);    
    return this->computeAlpha2(alpha, sequence, nbEmission);
}

int* HMM::estimateStatesSequence(int* sequence, const unsigned nbEmission) {
    double** delta = this->initMatrice(nbEmission, this->nbState);
    this->computeDelta(delta, sequence, nbEmission);
    return this->computeBestLikely(delta, nbEmission);
}

double** HMM::getA() const {
    return this->A;
}

double** HMM::getB() const {
    return this->B;
}

double* HMM::getPI() const {
    return this->PI;
}

int HMM::getNbState() const {
    return this->nbState;
}

int HMM::getNbObservation() const {
    return this->nbObservation;
}

void HMM::getA(double** A) {
    this->A = A;
}

void HMM::getB(double** B) {
    this->B = B;
}

void HMM::getPI(double* PI) {
    this->PI = PI;
}

void HMM::printA() const {
    this->printMatrice(this->A, this->nbState, this->nbState);
}

void HMM::printB() const {
    this->printMatrice(this->B, this->nbState, this->nbObservation);
}

void HMM::printPI() const {
    this->printDoubleVector(this->PI, this->nbState);
}

// Privates Methods

void HMM::printMatrice(double** matrice, const unsigned n, const unsigned m) const {
    for(unsigned i=0;i<n;++i) {
        for(unsigned j=0;j<m;++j) {
            cerr << matrice[i][j] << " ";
        }
        cerr << endl;
    }
}

void HMM::printDoubleVector(double* vect, const unsigned n) const {
    for(unsigned i=0;i<n;++i) {
        cerr << vect[i] << " ";
    }
    cerr << endl;
}

void HMM::printIntVector(int* vect, const unsigned n) const {
    for(unsigned i=0;i<n;++i) {
        cerr << vect[i] << " ";
    }
    cerr << endl;
}

double HMM::computeAlpha(double** alpha, int* sequence, double* C, const unsigned nbEmission) {
    C[0] = 0;
    unsigned i;
    for(i=0;i<this->nbState;++i) {
        alpha[0][i] = this->PI[i] * this->B[i][sequence[0]];
        C[0] += alpha[0][i];
    }
    /*if(C[0] <= 0) {
        C[0] = 1;
    }*/
    C[0] = 1 / C[0];
    double logProb = log(C[0]);
    for(i=0;i<this->nbState;++i) {
        alpha[0][i] *= C[0];
    }
    // Récursivité
    for(unsigned t=1;t<nbEmission;++t) {
        C[t] = 0;
        for(i=0;i<this->nbState;++i) {
            alpha[t][i] = 0;
            for(unsigned j=0;j<this->nbState;++j) {
                alpha[t][i] += this->A[j][i] * alpha[t-1][j];
            }
        }
        for(i=0;i<this->nbState;++i) {
            alpha[t][i] *= this->B[i][sequence[t]];
            C[t] += alpha[t][i];
        }
        /*if(C[t] <= 0) {
            C[t] = 1;
        }*/
        C[t] = 1/C[t];
        logProb += log(C[t]);
        for(i=0;i<this->nbState;++i) {
            alpha[t][i] *= C[t];  
        }
    }
    return -logProb;
}

double HMM::computeAlpha2(double** alpha, int* sequence, const unsigned nbEmission) {
    unsigned i;
    for(i=0;i<this->nbState;++i) {
        alpha[0][i] = this->PI[i] * this->B[i][sequence[0]];
    }

    double p = 0;
    // Récursivité
    for(unsigned t=1;t<nbEmission;++t) {
        for(i=0;i<this->nbState;++i) {
            alpha[t][i] = 0;
            for(unsigned j=0;j<this->nbState;++j) {
                alpha[t][i] += this->A[j][i] * alpha[t-1][j];
            }
        }
        for(i=0;i<this->nbState;++i) {
            alpha[t][i] *= this->B[i][sequence[t]];
        }
    }

    for(i=0;i<this->nbState;++i) {
        p += alpha[nbEmission-1][i];
    }

    return p;
}

void HMM::computeBeta(double** beta, double* C, int* sequence, const unsigned nbEmission) {
    unsigned i;
    for(i=0;i<this->nbState;++i) {
        beta[nbEmission-1][i] = C[nbEmission-1];
    }
    for(int t = nbEmission-2;t>=0;t--) {
        for(i=0;i<this->nbState;++i) {
            beta[t][i] = 0;
            for(unsigned j=0;j<this->nbState;++j) {
                beta[t][i] +=  beta[t+1][j] * this->A[i][j] * this->B[j][sequence[t+1]];
            }
            beta[t][i] *= C[t];
        }
    }
}

void HMM::computeDelta(double** delta, int* sequence, const unsigned nbEmission) {
    double max = 0;
    double c;
    int ind = -1;
    unsigned i;
    for(i=0;i<this->nbState;++i) {
        delta[0][i] = this->PI[i] * this->B[i][sequence[0]];
        if(delta[0][i] > max) {
            max = delta[0][i];
            ind = i;
        }
    }
    for(unsigned t=1;t<nbEmission;++t) {
        for(i=0;i<this->nbState;++i) {
            max = 0;
            ind = -1;
            for(unsigned j=0;j<this->nbState;++j) {
                c = this->A[j][i] * delta[t-1][j] * this->B[i][sequence[t]];
                if(c > max) {
                    max = c;
                    ind = j;
                }
            }
            delta[t][i] = max;
        }
    }
}

void HMM::computeDiGammaAndGamma(double*** diGamma, double** gamma, int* sequence, double** alpha, double** beta, const unsigned nbEmission) {
    double bottom = 0;
    double top;
    unsigned i, j;
    for(unsigned t=0;t<nbEmission-1;++t) {
        bottom = 0;
        for(i=0;i<this->nbState;++i) {
            for(j=0;j<this->nbState;++j) {
                bottom += alpha[t][i] * this->A[i][j] * this->B[j][sequence[t+1]] * beta[t+1][j];
            }
        }
        for(i=0;i<this->nbState;++i) {
            gamma[t][i] = 0;
            for(j=0;j<this->nbState;++j) {
                top = alpha[t][i] * this->A[i][j] * this->B[j][sequence[t+1]] * beta[t+1][j];
                if(bottom > 0)
                    diGamma[t][i][j] = top/bottom;
                gamma[t][i] += diGamma[t][i][j];
            }
        }
    }
    bottom = 0;
    for(i=0;i<this->nbState;++i) {
        bottom += alpha[nbEmission-1][i];
    }
    if(bottom > 0) {
        for(i=0;i<this->nbState;++i) {
            gamma[nbEmission-1][i] = alpha[nbEmission-1][i]/bottom;
        }
    }
}

void HMM::transitionEstimates(double*** diGamma, double** gamma, double** alpha, const unsigned nbEmission) {
    double top = 0;
    double bottom = 0;
    for(int i=0;i<this->nbState;i++) {
        for(int j=0;j<this->nbState;j++) {
            top = bottom = 0;
            for(int t=0;t<nbEmission-1;t++) {
                top += diGamma[t][i][j];
                bottom += gamma[t][i];
            }
            if(bottom > 0)
                A[i][j] = top/bottom;
        }
    }
}

void HMM::emissionEstimates(double*** diGamma, double** gamma, int* sequence, const unsigned nbEmission) {
    double top = 0;
    double bottom = 0;
    for(int i=0;i<this->nbState;i++) {
        for(int k=0;k<this->nbObservation;k++) {
            top = bottom = 0;
            for(int t=0;t<nbEmission;t++) {
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

void HMM::transitionAndEmissionEstimates(double*** diGamma, double** gamma, int* sequence, const unsigned nbEmission) {
    double topA = 0;
    double topB = 0;
    double bottom = 0;
    for(unsigned i=0;i<this->nbState;++i) {
        for(unsigned j=0;j<this->nbState||j<this->nbObservation;++j) {
            topA = topB = bottom = 0;
            for(unsigned t=0;t<nbEmission-1;++t) {
                if(j<nbState)
                    topA += diGamma[t][i][j];
                bottom += gamma[t][i];
                if(sequence[t] == j) {
                    topB += gamma[t][i];
                }
            }
            if(bottom > 0) {
                if(j<nbObservation)
                    this->B[i][j] = topB/bottom;
                if(j<nbState)
                    this->A[i][j] = topA/bottom;
            } 
        }
    }
}

double HMM::computeLog(double * C, const unsigned nbEmission) {
    double logProb = 0;
    for(unsigned i=0;i<nbEmission;i++) {
        logProb += log(C[i]);
    }
    return -logProb;
}

int* HMM::computeBestLikely(double** delta, const unsigned nbEmission) {
    int *best = new int[nbEmission];
    int l = nbEmission - 1;
    double max = 0;
    int ind = -1;
    for(unsigned j=0;j<this->nbState;++j) {
        if(delta[l][j] > max) {
            max = delta[l][j];
            ind = j;
        }
    }

    if(ind == -1)
        return nullptr;

    int index = nbEmission-1;
    best[index] = ind;
    unsigned b = ind;
    for(unsigned i=nbEmission-1;i>=1;i--) {
        index--;
        max = 0;
        ind = -1;
        for(unsigned j=0;j<this->nbState;++j) {
            if(delta[i-1][j] * this->A[j][b] > max) {
                max = delta[i-1][j] * this->A[j][b];
                ind = j;
            }
        }
        b = ind;
        best[index] = ind;
    }
    return best;
}

double HMM::computeProbability(double** alpha, double* C, const unsigned nbEmission) {
    double p = 0;
    double t = 0;
    for(unsigned i=0;i<this->nbState;++i) {
        t = alpha[nbEmission-1][i];
        for(unsigned i=0;i<nbEmission;++i)
            t /= C[i];
        p += t;
    }
    return p;
} 

double** HMM::initRandomMatrice(const unsigned n, const unsigned m) {
    double** matrice = new double*[n];
    for(unsigned i=0;i<n;++i) {
        matrice[i] = this->initRandomVector(m);
    }
    return matrice;
}

double* HMM::initRandomVector(const unsigned n) {
    uniform_real_distribution<double> dist(0.01, 1.0/n);  //(min, max)
    //Mersenne Twister: Good quality random number generator
    std::mt19937 rng; 
    //Initialize with non-deterministic seeds
    //rng.seed(std::random_device{}());

    double *vect = new double[n];
    for(unsigned i=0;i<n;++i) {
        vect[i] = 1.0/n;
    }
    for(unsigned i=0;i<rand() % 10+1;++i) {
        unsigned i1;
        unsigned i2;
        do {
            i1 = rand() % (n-1);
            i2 = rand() % (n-1);
        } while(i1 == i2);
        double epsilon = dist(rng);
        if(vect[i1] - epsilon > 0 && vect[i2] + epsilon < 1) {
            vect[i1] -= epsilon;
            vect[i2] += epsilon; 
        }

    }
    return vect;
}

double** HMM::initMatrice(const unsigned n, const unsigned m, double value) {
    double** matrice = new double*[n];
    for(unsigned i=0;i<n;++i) {
        matrice[i] = new double[m];
        for(unsigned j=0;j<m;++j) {
            matrice[i][j] = value;
        }
    }
    return matrice;
}

double* HMM::initVector(const unsigned n, double value) {
    double *vect = new double[n];
    for(unsigned i=0;i<n;++i) {
        vect[i] = value;
    }
    return vect;
}


double** HMM::initBMatrice(const unsigned n, const unsigned m) {
    double** matrice = new double*[n];
    for(unsigned i=0;i<n;++i) {
        matrice[i] = new double[m];
    }
    //             UP_L  UP   UP_R   L    S      R    D_L    D  D_R
    
    //Migrating left
    //UP_L
    matrice[0][0] = 1.0/4;
    //UP
    matrice[0][1] = 0.0;
    // UP_R
    matrice[0][2] = 0.0;
    //L
    matrice[0][3] = 1.0/4;
    //S
    matrice[0][4] = 1.0/4;
    //R
    matrice[0][5] = 0.0;
    //D_L
    matrice[0][6] = 1.0/4;
    //D
    matrice[0][7] = 0.0;
    //D_R
    matrice[0][8] = 0.0;

    //Migrating right
    //UP_L
    matrice[1][0] = 0.0;
    //UP
    matrice[1][1] = 0.0;
    // UP_R
    matrice[1][2] = 1.0/4;
    //L
    matrice[1][3] = 0.0;
    //S
    matrice[1][4] = 1.0/4;
    //R
    matrice[1][5] = 1.0/4;
    //D_L
    matrice[1][6] = 0.0;
    //D
    matrice[1][7] = 0.0;
    //D_R
    matrice[1][8] = 1.0/4;


    //Circling
    //UP_L
    matrice[2][0] = 1.0/3;
    //UP
    matrice[2][1] = 0.0;
    // UP_R
    matrice[2][2] = 1.0/3;
    //L
    matrice[2][3] = 0.0;
    //S
    matrice[2][4] = 1.0/4;
    //R
    matrice[2][5] = 0.0;
    //D_L
    matrice[2][6] = 0.0;
    //D
    matrice[2][7] = 0.0;
    //D_R
    matrice[2][8] = 0.0;

    //Hunting
    //UP_L
    matrice[3][0] = 0.0;
    //UP
    matrice[3][1] = 0.0;
    // UP_R
    matrice[3][2] = 0.0;
    //L
    matrice[3][3] = 1.0/6;
    //S
    matrice[3][4] = 1.0/6;
    //R
    matrice[3][5] = 1.0/6;
    //D_L
    matrice[3][6] = 1.0/6;
    //D
    matrice[3][7] = 1.0/6;
    //D_R
    matrice[3][8] = 1.0/6;


    //Drilling random for the moment
    matrice[4] = initRandomVector(9);

    //Zig zag random
    matrice[5] = initRandomVector(9);

    return matrice;
}
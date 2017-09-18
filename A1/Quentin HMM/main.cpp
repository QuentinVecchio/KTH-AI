#include <iostream>
#include <cmath>
#include "Hmm.hpp"

double** getMatrice(const unsigned n, const unsigned m);
double* getDoubleVector(const unsigned n);
int* getIntVector(const unsigned n);

using namespace std;

int main() {
    unsigned nbState, nbObservation, nbSuiteObservation, fuck;
    double **A, **B;
    double *PI;
    int *Sequence;
    int i,j;
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
    HMM* hmm = new HMM(A, B, PI, nbState, nbObservation);
    // HMM 1
    //cout << hmm->emissionProbability(Sequence, nbSuiteObservation) << endl;

    // HMM 2
    /*int* best = hmm->estimateStatesSequence(Sequence, nbSuiteObservation);
    for(int i=0;i<nbSuiteObservation;i++) {
        cout << best[i] << " ";
    }
    cout << endl;*/
    
    // HMM 4
    hmm->learn(Sequence, nbSuiteObservation, 30);
    A = hmm->getA();
    B = hmm->getB();
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
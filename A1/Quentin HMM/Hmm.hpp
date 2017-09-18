#ifndef _HMM_HPP_
#define _HMM_HPP_


#include <iostream>
#include <cmath>

class HMM {
private:
    double **A, **B, *PI;
    int nbState, nbObservation;
    void printMatrice(double** matrice, const unsigned n, const unsigned m) const;
    void printDoubleVector(double* vect, const unsigned n) const;
    void printIntVector(int* vect, const unsigned n) const;
    double computeAlpha(double** alpha, int* sequence, double* C, const unsigned nbEmission);
    void computeBeta(double** beta, double* C, int* sequence, const unsigned nbEmission);
    void computeDelta(double** delta, int* sequence, const unsigned nbEmission);
    void computeDiGammaAndGamma(double*** diGamma, double** gamma, int* sequence, double** alpha, double** beta, const unsigned nbEmission);
    void transitionAndEmissionEstimates(double*** diGamma, double** gamma, int* sequence, const unsigned nbEmission);
    double computeLog(double * C, const unsigned nbEmission);
    int* computeBestLikely(double** delta, const unsigned nbEmission);
    double** initRandomMatrice(const unsigned n, const unsigned m);
    double* initRandomVector(const unsigned n);
    double** initMatrice(const unsigned n, const unsigned m, double value = 0);
    double** initBMatrice(const unsigned n, const unsigned m);
    double* initVector(const unsigned n, double value = 0);
    double computeProbability(double** alpha, double* C, const unsigned nbEmission);
    void transitionEstimates(double*** diGamma, double** gamma, double** alpha, const unsigned nbEmission);
    void emissionEstimates(double*** diGamma, double** gamma, int* sequence, const unsigned nbEmission);
public:
    HMM(double** A, double** B, double* PI, const unsigned nbState, const unsigned nbObservation);
    HMM(const unsigned nbState, const unsigned nbObservation);
    ~HMM();

    // Methods
    void learn(int* sequence, const unsigned nbEmission, const unsigned maxIteration);
    double emissionProbability(int *sequence, const unsigned nbEmission);
    int* estimateStatesSequence(int* sequence, const unsigned nbEmission);

    // Getters
    double** getA() const;
    double** getB() const;
    double* getPI() const;
    int getNbState() const;
    int getNbObservation() const;

    // Setters
    void getA(double** A);
    void getB(double** B);
    void getPI(double* PI);

    // Printers
    void printA() const;
    void printB() const;
    void printPI() const;
};

#endif

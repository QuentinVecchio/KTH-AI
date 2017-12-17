#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

#define NBSTATES 3
#define NBOBSERVATIONS 9
#define NBBIRDS 1

#define PROBA_MIN_SHOOT 0.75
#define OFFSET_SHOOT 8

#define PROBA_MIN_GUESS 0.2
#define OFFSET_GUESS 6

#define MAX_ITER 300

#define MIN_TURNS 85

namespace ducks
{

Player::Player()
{ 
    this->nbTurns = 0;
    this->nbShoot = 0;
    this->nbHit = 0;
    this->nbGuess = 0;
    this->nbGoodGuess = 0;
    this->nbBirds = 0;
    this->learningTab = new bool[this->nbModels];
    this->guessByModels = new int[this->nbModels];
    this->goodGuessByModels = new int[this->nbModels];
    this->nbSpecies = new int[this->nbModels+1];
    this->modelDistribution = new double*[this->nbModels];
    this->groupsMovements =  std::vector<std::vector<std::vector<int> > >(this->nbModels);
    this->groupsMovementsDistribution = std::vector<std::vector<std::vector<double> > >(this->nbModels);
    this->groupsMovementsDistributionNormalize = std::vector<std::vector<std::vector<double> > >(this->nbModels);
    this->beardInGroups = std::vector<int>(this->nbModels,0);
    this->bestCentroid = new unsigned[this->nbModels];
    this->bestDistance = new double[this->nbModels];
    for(unsigned i=0;i<this->nbModels;i++) {
        this->learningTab[i] = false;
        this->guessByModels[i] = 0;
        this->goodGuessByModels[i] = 0;
        this->nbSpecies[i] = 0;
    }
    this->brain = new int*[10];
    this->shootBirdModel = new int*[10];
    this->nbModels = COUNT_SPECIES;

    this->hmm = std::vector<std::vector<HMM*> >(this->nbModels);
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    double pMax = 0.0;
    int movMax = -1;
    int birdMax = -1;
    int hmmMax = -1;

    ++nbTurns;

    if(pState.getRound() == 0)
        return cDontShoot;
    /*
    * ROUND 0
    */

    /*
    * SHOOT TURNS
    */
    if(this->nbTurns > MIN_TURNS) {
        for(unsigned i = 0; i < pState.getNumBirds(); ++i) {
            if(pState.getBird(i).isAlive()) {
                int* movements = new int[pState.getBird(i).getSeqLength()];
                for(unsigned j = 0; j < pState.getBird(i).getSeqLength(); ++j) {
                    movements[j] = pState.getBird(i).getObservation(j);
                }
                if(this->guessSpecie(movements, pState.getBird(i).getSeqLength()) == SPECIES_BLACK_STORK) {
                    std::cerr << "The bird is black" << std::endl;
                } else {

                    for(int m=0;m<this->nbModels-1;m++) {
                        for(int m1=0;m1<this->hmm[m].size();m1++) {
                            int* states = this->hmm[m][m1]->estimateStatesSequence(movements, pState.getBird(i).getSeqLength());
                            if(states == nullptr)
                                continue;
                            for(int k = 0; k < COUNT_MOVE; ++k) {
                                double tmp = 0;
                                for(unsigned j = 0; j < this->hmm[m][m1]->getNbState(); ++j) {
                                    tmp += this->hmm[m][m1]->getA()[states[pState.getBird(i).getSeqLength()-1]][j] * this->hmm[m][m1]->getB()[j][k];
                                }
                                if(tmp > pMax) {
                                    pMax = tmp;
                                    movMax = k;
                                    birdMax = i;
                                    hmmMax = m;
                                }
                            }
                            delete[] states;
                        }
                    }
                }
            }
        }
    }

        
    if(pMax > PROBA_MIN_SHOOT && birdMax != -1 && movMax != -1) {
        std::cerr << "######### SHOOT ON Bird N° " << birdMax << " at movement " << movMax << " with hmm " << hmmMax <<  " and p " << pMax << std::endl;
        this->nbShoot++;
        return Action(birdMax, EMovement(movMax));
    }
    else {
        if(this->nbTurns > MIN_TURNS)
            std::cerr << "######### BEST is bird : " << birdMax << " at movement " << movMax << " with hmm " << hmmMax <<  " and p " << pMax << std::endl;
    }


    //std::cerr << "Round " << pState.getNumNewTurns() << " " << pState.getRound() << std::endl;
    return cDontShoot;
}

double Player::isBlack(int* movements, int size) const {
    double pMin = -1;
    // Get the movements for the bird until his death

    for(unsigned m=0;m<this->hmm[this->nbModels-1].size();m++) {
        double p = this->hmm[this->nbModels-1][m]->emissionProbability(movements, size);
        if(p>pMin) {
            pMin = p;
        }
    }

    return pMin;
}

int Player::guessSpecie(int* movements, int size) const {
    int species = -1;
    double pMin = -1;
    for(unsigned k = 0; k < this->nbModels; ++k) { 
        for(unsigned m=0;m<this->hmm[k].size();m++) {
            double p = this->hmm[k][m]->emissionProbability(movements, size);
            if(p>pMin) {
                pMin = p;
                species = ESpecies(k);
            }
        }
    }
    return species;
}


void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */
    std::cerr << "######################## BIRD " << pBird << " KILLED" << std::endl;
    this->nbHit++;
}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    // Statistics
    this->guesses = new int[pState.getNumBirds()];
    this->nbBirds += pState.getNumBirds();

    // First round we do nothing
    if(pState.getRound() == 0) {
        std::vector<ESpecies> lGuesses(pState.getNumBirds(), ESpecies(0)); // Guess randomly
        return lGuesses;
    } else {
        // Search the max to normalize
        std::cerr << "================ TIME TO GUESS ================" << std::endl;
        std::vector<ESpecies> lGuesses(pState.getNumBirds(), ESpecies(-1));
        
        // Try to guess each birds
        for(unsigned i=0;i<pState.getNumBirds();++i) {
            double pMin = -1;
            // Get the movements for the bird until his death
            int* movements = new int[pState.getBird(i).getSeqLength()];
            unsigned nbMovements = 0;
            for(unsigned j = 0; j < pState.getBird(i).getSeqLength(); ++j) {
                nbMovements++;
                if(pState.getBird(i).getObservation(j) == -1) {
                    --nbMovements;
                    break;
                }
                movements[j] = pState.getBird(i).getObservation(j);
            }

            for(unsigned k = 0; k < this->nbModels; ++k) { 
                for(unsigned m=0;m<this->hmm[k].size();m++) {
                    double p = this->hmm[k][m]->emissionProbability(movements, nbMovements);
                    if(p>pMin) {
                        pMin = p;
                        lGuesses[i] = ESpecies(k);
                    }
                }
            }

            // Statistic
            if(lGuesses[i] != -1) {
                this->nbGuess++;
                this->guessByModels[lGuesses[i]]++;
            }
        }
        // Statistic
        for(unsigned i=0;i<lGuesses.size();i++) {
            this->guesses[i] = lGuesses[i];
            std::cerr << lGuesses[i] << " ";
        }
        std::cerr << std::endl;
        return lGuesses;
    }
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    this->nbTurns = 0;
    std::vector<std::vector<int> > groups = std::vector<std::vector<int> >(this->nbModels);
    for(unsigned i = 0;i<pSpecies.size();++i) {
        /*
        *   statistiques
        */
        if(pSpecies[i] != -1)
            this->nbSpecies[pSpecies[i]]++;
        else
            this->nbSpecies[this->nbModels]++;

        if(pState.getRound() > 0 && pSpecies[i] == this->guesses[i] && pSpecies[i] != -1) {
            this->nbGoodGuess++;
            this->goodGuessByModels[pSpecies[i]]++;
        }

        /*
        *   Concatenation mouvements
        */
        if(pSpecies[i] >= 0) {
            for(unsigned j = 0; j < pState.getBird(i).getSeqLength(); ++j) {
                if(pState.getBird(i).getObservation(j) >= 0)
                    groups[pSpecies[i]].push_back(pState.getBird(i).getObservation(j));
            }
        }
    }

    for(unsigned i=0;i<this->nbModels;i++) {

        if(groups[i].size() >= 70) {
            int* movements = new int[groups[i].size()];
            for(unsigned j=0;j<groups[i].size();j++) {
                movements[j] = groups[i][j];
            }

            HMM* nHmm = new HMM(NBSTATES, NBOBSERVATIONS);

            nHmm->learn(movements, groups[i].size(), 300);

            this->hmm[i].push_back(nHmm);
        }
    }

    this->printStatics();
}

void Player::printStatics() const {
    std::cerr << "LEARNING STATS ===========================" << std::endl;
        for(unsigned i=0;i<this->nbModels;i++) {
            std::cerr << this->learningTab[i] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "ACTIONS STATS ===========================" << std::endl;
        std::cerr << "Nb Shoot : " << this->nbShoot << std::endl;
        std::cerr << "Nb Hit : " << this->nbHit << std::endl;
        std::cerr << "Ratio : " << double(this->nbHit)/this->nbShoot << std::endl;
        std::cerr << "Nb Birds : " << this->nbBirds << std::endl;
        std::cerr << "Ratio : " << double(this->nbHit)/this->nbBirds << std::endl;
        std::cerr << "Nb guess : " << this->nbGuess << std::endl;
        std::cerr << "Nb good guess : " << this->nbGoodGuess << std::endl;
        std::cerr << "Ratio : " <<  double(this->nbGoodGuess)/this->nbGuess << std::endl;
        std::cerr << "Ratio GuessDid : " <<  double(this->nbGuess)/this->nbBirds << std::endl;
        std::cerr << "Ratio GoodGuessDid : " <<  double(this->nbGoodGuess)/this->nbBirds << std::endl;
        std::cerr << "Stats species :  \t\t";
        std::cerr.precision(3);
        double nbPointLost = 0;
        for(unsigned x=0;x<=this->nbModels;x++)
            std::cerr << this->nbSpecies[x] << "  \t\t";
        std::cerr << std::endl;
        std::cerr << "Stats guess models :  \t\t";
        for(unsigned x=0;x<this->nbModels;x++)
            std::cerr << this->guessByModels[x] << "  \t\t";
        std::cerr << std::endl;
        std::cerr << "Stats good guess models : \t";
        for(unsigned x=0;x<this->nbModels;x++) {
            std::cerr << this->goodGuessByModels[x] << " \t\t";
            nbPointLost +=  this->guessByModels[x] - this->goodGuessByModels[x];
        }
        std::cerr << std::endl;
        std::cerr << "Ratio : \t \t \t";
        for(unsigned x=0;x<this->nbModels;x++)
            std::cerr << double(this->goodGuessByModels[x])/this->guessByModels[x] << "  \t\t";
        std::cerr << std::endl << "Nb points perdus : " << nbPointLost;
        std::cerr << std::endl << "================================================" << std::endl;
}

} /*namespace ducks*/
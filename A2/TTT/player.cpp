#include "player.hpp"
#include <cstdlib>
#include <limits>
#include <iostream>
#include <stdint.h>

#define LARGEUR 4
#define HAUTEUR 4

#define DEPTH 4

#define MINIMUM std::numeric_limits<int>::min()
#define MAXIMUM std::numeric_limits<int>::max()

namespace TICTACTOE
{

Player::Player() {
    this->states[CELL_X] = std::map<std::string, int>();
    this->states[CELL_O] = std::map<std::string, int>();
}

GameState Player::play(const GameState &pState,const Deadline &pDue)
{
    this->states[CELL_X].clear();
    this->states[CELL_O].clear();
    this->count = 0;
    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

    if (lNextStates.size() == 0) return GameState(pState, Move());

    unsigned indexBestMove = 0;
    int bestEvaluation = MINIMUM;

    std::cerr << "Analyze possibilities" << std::endl;
    for(unsigned i=0;i<lNextStates.size();++i) {
        //int v = this->minMax(lNextStates[i], pState.getNextPlayer(), 3);
        int v = this->alphaBeta(lNextStates[i], pState.getNextPlayer(), 15, bestEvaluation, MAXIMUM);
        if(v>bestEvaluation) {
            bestEvaluation = v;
            indexBestMove = i;
        }
    }

    std::cerr << "Best movements find : " << indexBestMove << " with a score of " << bestEvaluation << std::endl;
    std::cerr << "Nb de node visités " << this->count << std::endl;
    return lNextStates[indexBestMove];
}

int Player::gain(int nbCell, bool currentPlayer) const {
    int coef = 1;
    if(currentPlayer && nbCell == 4)
        coef = -1;
    else if(currentPlayer)
        coef = 0;
    switch(nbCell) {
        case 1:
            return coef;
        case 2:
            return coef*10;
        case 3:
            return coef*100;
        case 4:
            return coef*1000;
        default:
            return 0;
    }
}

int Player::eval(GameState state, uint8_t player) const {
    int sum = 0;
    Cell valueMine = Cell(player);
    Cell valueOpp = CELL_O;
    if(valueMine == CELL_O)
        valueOpp = CELL_X;
    int nbMine;
    int nbOpp;
    //rows
    for(unsigned i=0;i<4;i++) {
        nbMine = 0;
        nbOpp = 0;
        for(unsigned j=0;j<4;j++) {
            if(state.at(i, j) == valueMine) {
                nbMine++;
            } else if(state.at(i, j) == valueOpp) {
                nbOpp++;
            }
        }
        sum += this->gain(nbMine, player == state.getNextPlayer());
        sum += this->gain(nbOpp, player == state.getNextPlayer());
    }
    //columns
    for(unsigned i=0;i<4;i++) {
        nbMine = 0;
        nbOpp = 0;
        for(unsigned j=0;j<4;j++) {
            if(state.at(j, i) == valueMine) {
                nbMine++;
            } else if(state.at(j, i) == valueOpp) {
                nbOpp++;
            }
        }
        sum += this->gain(nbMine, player == state.getNextPlayer());
        sum += this->gain(nbOpp, player == state.getNextPlayer());
    }
    ///diagonals left-top to rigth-bottom
    for(unsigned i=0;i<4;i++) {
        nbMine = 0;
        nbOpp = 0;
        if(state.at(i, i) == valueMine) {
            nbMine++;
        } else if(state.at(i, i) == valueOpp) {
            nbOpp++;
        }
        sum += this->gain(nbMine, player == state.getNextPlayer());
        sum += this->gain(nbOpp, player == state.getNextPlayer());
    }
    //diagonals rigth-top to left-bottom
    for(unsigned i=0;i<4;i++) {
        nbMine = 0;
        nbOpp = 0;
        if(state.at(i, 3-i) == valueMine) {
            nbMine++;
        } else if(state.at(i, 3-i) == valueOpp) {
            nbOpp++;
        }
        sum += this->gain(nbMine, player == state.getNextPlayer());
        sum += this->gain(nbOpp, player == state.getNextPlayer());
    }
    return sum;
}

int Player::minMax(GameState state, uint8_t player, int depth) {
    uint8_t nextPlayer = state.getNextPlayer();
    std::vector<GameState> possibleMoves;
    state.findPossibleMoves(possibleMoves);
    if(state.isEOG() || depth == 0) {
        return this->eval(state, player);
    } else {
        this->count++;
        int bestPossible;
        if(player != nextPlayer) {
            bestPossible = MINIMUM;
            for(unsigned i=0;i<possibleMoves.size();++i) {
                int tmp;
                if(!this->stateExist(player, possibleMoves[i])) {
                    tmp = minMax(possibleMoves[i], player, depth-1);
                    this->createState(player, possibleMoves[i], tmp);
                } else {
                    tmp = this->getGainForState(player, possibleMoves[i]);
                }
                if(tmp>bestPossible) {
                    bestPossible = tmp;
                }
            }
            return bestPossible;
        } else {
            bestPossible = MAXIMUM;
            for(unsigned i=0;i<possibleMoves.size();++i) {
                int tmp;
                if(!this->stateExist(player, possibleMoves[i])) {
                    tmp = minMax(possibleMoves[i], player, depth-1);
                    this->createState(player, possibleMoves[i], tmp);
                } else {
                    tmp = this->getGainForState(player, possibleMoves[i]);
                }
                if(tmp<bestPossible) {
                    bestPossible = tmp;
                }
            }
        }
        return bestPossible;
    }
}

int Player::alphaBeta(GameState state, uint8_t player, unsigned depth, int alpha, int beta) {
    uint8_t nextPlayer = state.getNextPlayer();
    int nextDepth = depth;
    nextDepth--;
    /*if(factoriel(16-dept) > 3000) {
        
    }*/
    std::vector<GameState> possibleMoves;
    state.findPossibleMoves(possibleMoves);
    this->count++;
    if(depth == 0 || state.isEOG()) {
        return this->eval(state, player);
    } else if(player == nextPlayer) {
        int max = MINIMUM;
        for(unsigned i=0;i<possibleMoves.size();++i) {
            int tmp;
            if(!this->stateExist(player, possibleMoves[i])) {
                tmp = alphaBeta(possibleMoves[i], player, nextDepth, alpha, beta);
                this->createState(player, possibleMoves[i], tmp);
            } else {
                tmp = this->getGainForState(player, possibleMoves[i]);
            }
            if(tmp > max) {
                max = tmp;
            }
            if(max > alpha) {
                alpha = max;
            }
            if(beta<=alpha) {
                break;
            }
        }
        return max;
    } else {
        int min = MAXIMUM;
        for(unsigned i=0;i<possibleMoves.size();++i) {
            int tmp;
            if(!this->stateExist(player, possibleMoves[i])) {
                tmp = alphaBeta(possibleMoves[i], player, nextDepth, alpha, beta);
                this->createState(player, possibleMoves[i], tmp);
            } else {
                tmp = this->getGainForState(player, possibleMoves[i]);
            }
            if(tmp<min) {
                min = tmp;
            }
            if(min<beta) {
                beta = min;
            }
            if(beta<=alpha) {
                break;
            }
        }
        return min;
    }
}

bool Player::stateExist(unsigned player, GameState state) const {
    // Test the first matrice    
    int **matrice = new int*[4];
    for(unsigned x=0;x<4;x++) {
        matrice[x] = new int[4];
        for(unsigned y=0;y<4;y++) {
            matrice[x][y] = state.at(x, y);
        }
    }
    std::stringstream ss;
    for(unsigned x=0;x<4;x++) {
        for(unsigned y=0;y<4;y++) {
            ss << MESSAGE_SYMBOLS[matrice[x][y]];
        }
    }
    if(this->states.find(player)->second.find(ss.str()) != this->states.find(player)->second.end()) {
        //std::cerr << "Matrice exists " << std::endl;
        return true;
    }
    // Test for the rotated matrices
    for(unsigned i=0;i<3;i++) {
        int** matriceRotated = this->matriceRotation(matrice, 4, 4);
        ss.str(""); // clean the string stream
        for(unsigned x=0;x<4;x++) {
            for(unsigned y=0;y<4;y++) {
                ss << MESSAGE_SYMBOLS[matriceRotated[x][y]];
            }
        }
        if(this->states.find(player)->second.find(ss.str()) != this->states.find(player)->second.end()) {
            //std::cerr << "Matrice rotated exists " << std::endl;
            return true;
        } else {
            matrice = matriceRotated;
        }
    }
    return false;
}

int Player::getGainForState(unsigned player, GameState state) const {
    std::stringstream ss;
    for(int i=0;i<state.cSquares;i++)
        ss << MESSAGE_SYMBOLS[state.at(i)];
    return this->states.find(player)->second.find(ss.str())->second;
}

void Player::createState(unsigned player, GameState state, int value) {
    std::stringstream ss;
    for(int i=0;i<state.cSquares;i++)
        ss << MESSAGE_SYMBOLS[state.at(i)];
    //std::cerr << "Analyse " << ss.str() << std::endl;
    this->states.find(player)->second.insert(std::pair<std::string, int>(ss.str(), value));
}

int** Player::matriceRotation(int** matrice, int n, int m) const {
    int** newMatrice = new int*[m];
    for(unsigned i =0;i<n;i++) {
        newMatrice[i] = new int[n];
    }
    for(unsigned i=0;i<m;i++) {
        for(unsigned j=0;j<n;j++) {
            newMatrice[i][j] = matrice[n-j-1][i];
        }
    }
    return newMatrice;
}

int Player::factoriel(int n) {
    int f = 1;
    for(int i=n;i>0;i--) {
        f *= i;
    }
    return f;
}
/*namespace TICTACTOE*/ 
}

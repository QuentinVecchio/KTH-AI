#include "player.hpp"
#include <cstdlib>
#include <climits>

namespace TICTACTOE
{

GameState Player::play(const GameState &pState,const Deadline &pDue)
{
    this->states.clear();
    this->count = 0;
    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

    if (lNextStates.size() == 0) return GameState(pState, Move());

    unsigned indexBestMove = 0;
    unsigned bestEvaluation = 0;

    std::cerr << "Analyze possibilities" << std::endl;
    for(unsigned i=0;i<lNextStates.size();++i) {
        //int v = this->minMax(lNextStates[i], 1, 4);
        int v = this->alphaBeta(lNextStates[i], 5, INT_MIN, INT_MAX, 1);
        if(v>bestEvaluation) {
            bestEvaluation = v;
            indexBestMove = i;
        }
    }

    std::cerr << "Best movements find : " << indexBestMove << " with a score of " << bestEvaluation << std::endl;
    std::cerr << "Nb de node visités " << this->count << std::endl;
    return lNextStates[indexBestMove];
}

int Player::gain(int nbCell, int player) const {
    int coef = 1;
    if(player == 1)
        coef = -1;
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

int Player::eval(GameState state, int player) const {
    int sum = 0;
    Cell valueMine = CELL_X;
    Cell valueOpp = CELL_O;
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
        sum += this->gain(nbMine, 0);
        sum += this->gain(nbOpp, 1);
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
        sum += this->gain(nbMine, 0);
        sum += this->gain(nbOpp, 1);
    }
    //diagonals left-top to rigth-bottom
    for(unsigned i=0;i<4;i++) {
        nbMine = 0;
        nbOpp = 0;
        if(state.at(i, i) == valueMine) {
            nbMine++;
        } else if(state.at(i, i) == valueOpp) {
            nbOpp++;
        }
        sum += this->gain(nbMine, 0);
        sum += this->gain(nbOpp, 1);
    }
    //diagonals rigth-top to left-bottom
    for(unsigned i=0;i<4;i++) {
        nbMine = 0;
        nbOpp = 0;
        if(state.at(i, 4-i) == valueMine) {
            nbMine++;
        } else if(state.at(i, 4-i) == valueOpp) {
            nbOpp++;
        }
        sum += this->gain(nbMine, 0);
        sum += this->gain(nbOpp, 1);
    }
    return sum;
}

int Player::minMax(GameState state, int player, int depth) {
    if(state.isEOG() || depth == 0) {
        return this->eval(state, player);
    } else {
        this->count++;
        int bestPossible;
        if(player == 0) {
            bestPossible = INT_MIN;
            std::vector<GameState> possibleMoves;
            state.findPossibleMoves(possibleMoves);
            for(unsigned i=0;i<possibleMoves.size();++i) {
                int v;
                if(!this->stateExist(possibleMoves[i])) {
                    v = minMax(possibleMoves[i], 0, depth-1);
                    this->createState(state, v);
                } else {
                    v = this->getGainForState(possibleMoves[i]);
                }
                if(v>bestPossible) {
                    bestPossible = v;
                }
            }
            return bestPossible;
        } else {
            bestPossible = INT_MAX;
            std::vector<GameState> possibleMoves;
            state.findPossibleMoves(possibleMoves);
            for(unsigned i=0;i<possibleMoves.size();++i) {
                int v;
                if(!this->stateExist(possibleMoves[i])) {
                    v = minMax(possibleMoves[i], 0, depth-1);
                    this->createState(state, v);
                } else {
                    v = this->getGainForState(possibleMoves[i]);
                }
                if(v<bestPossible) {
                    bestPossible = v;
                }
            }
        }
        return bestPossible;
    }
}

int Player::alphaBeta(GameState state, int depth, int alpha, int beta, int player) {
    //return 0;
    int v;
    if(depth == 0 || state.isEOG()) {
        return this->eval(state, player);
    } else if(player == 0) {
        v = INT_MIN;
        this->count++;
        std::vector<GameState> possibleMoves;
        state.findPossibleMoves(possibleMoves);
        for(unsigned i=0;i<possibleMoves.size();++i) {
            if(!this->stateExist(possibleMoves[i])) {
                int val;
                if(!this->stateExist(possibleMoves[i])) {
                    val = alphaBeta(possibleMoves[i], depth-1, alpha, beta, 1);
                    this->createState(state, val);
                } else {
                    val = this->getGainForState(possibleMoves[i]);
                }
                if(val>v) {
                    v = val;
                }
                if(v>alpha) {
                    alpha = v;
                }
                if(beta<=alpha) {
                    return v;
                }
            }
        }
    } else if(player == 1) {
        v = INT_MAX;
        std::vector<GameState> possibleMoves;
        state.findPossibleMoves(possibleMoves);
        for(unsigned i=0;i<possibleMoves.size();++i) {
            if(!this->stateExist(possibleMoves[i])) {
                int val;
                if(!this->stateExist(possibleMoves[i])) {
                    val = alphaBeta(possibleMoves[i], depth-1, alpha, beta, 0);
                    this->createState(state, val);
                } else {
                    val = this->getGainForState(possibleMoves[i]);
                }
                if(val<v) {
                    v = val;
                }
                if(v<beta) {
                    beta = v;
                }
                if(beta<=alpha) {
                    return v;
                }
            }
        }
    }
    return v;
}

bool Player::stateExist(GameState state) const {
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
    if(this->states.find(ss.str()) != this->states.end()) {
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
        if(this->states.find(ss.str()) != this->states.end()) {
            //std::cerr << "Matrice rotated exists " << std::endl;
            return true;
        } else {
            matrice = matriceRotated;
        }
    }
    return false;
}

int Player::getGainForState(GameState state) const {
    std::stringstream ss;
    for(int i=0;i<state.cSquares;i++)
        ss << MESSAGE_SYMBOLS[state.at(i)];
    return this->states.find(ss.str())->second;
}

void Player::createState(GameState state, int value) {
    std::stringstream ss;
    for(int i=0;i<state.cSquares;i++)
        ss << MESSAGE_SYMBOLS[state.at(i)];
    //std::cerr << "Analyse " << ss.str() << std::endl;
    this->states.insert(std::pair<std::string, int>(ss.str(), value));
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


/*namespace TICTACTOE*/ 
}

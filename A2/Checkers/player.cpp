#include "player.hpp"
#include <cstdlib>
#include <limits>
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <algorithm>

#define MINIMUM -100000
#define MAXIMUM std::numeric_limits<double>::max()

namespace checkers
{

unsigned Player::countPieces(const GameState &pState) {
    unsigned sum = 0;
    std::string grid = pState.toMessage();
    for(unsigned i = 0; i < grid.size(); ++i) {
        if(grid[i] == 'r' || grid[i] == 'R' || grid[i] == 'w' || grid[i] == 'W')
            ++sum;
    }
    return sum;
}

std::string Player::boardToString(const GameState & pState) const {
    std::stringstream ss;
    for(unsigned i=0;i<pState.cSquares;i++) {
        ss << SIMPLE_TEXT[pState.at(i)];
    }

    return ss.str();
}

std::string Player::boardToMessage(const GameState & pState) const {
    std::stringstream ss;
    for(unsigned i=0;i<pState.cSquares;i++) {
        ss << SIMPLE_TEXT[pState.at(i)];
    }
    ss << pState.getNextPlayer();

    return ss.str();
}



double Player::computeHeuristic(const GameState &pState, const uint8_t & player, bool currentPlayer) {
    if(pState.isDraw())
        return -1000;
    if(pState.isRedWin() || pState.isWhiteWin()) {
        if(currentPlayer) {
            return 1000;
        }
        else {
            return -1000;
        }
    }

    double nbRed = 0;
    double nbWhite = 0;

    if(currentPlayer) {
        for(int i=0;i<8;i++) {
            for(int j=0;j<8;j++) {
                if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'r')
                    nbRed += i;
                else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')
                    nbRed += i+ 4;
                else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'w')
                    nbWhite += 7 - i;
                else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'W')
                    nbWhite += 7 - i + 4;

                if(i>0 && j>0 && MESSAGE_SYMBOLS[pState.at(i-1, j-1)] == '.') {
                    if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'r' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')
                        nbRed++;
                    else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'w' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'W')
                        nbWhite++;
                }
                if(i>0 && j<7 && MESSAGE_SYMBOLS[pState.at(i-1, j+1)] == '.') {
                    if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'r' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')
                        nbRed++;
                    else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'w' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'W')
                        nbWhite++;
                }
                if(i<7 && j>0 && MESSAGE_SYMBOLS[pState.at(i+1, j-1)] == '.') {
                    if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'r' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')
                        nbRed++;
                    else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'w' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'W')
                        nbWhite++;
                }
                if(i<7 && j<7 && MESSAGE_SYMBOLS[pState.at(i+1, j+1)] == '.') {
                    if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'r' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')
                        nbRed++;
                    else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'w' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'W')
                        nbWhite++;
                }

                if((MESSAGE_SYMBOLS[pState.at(i, j)] == 'r' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'R') && (j==0 || j == 7)){
                    nbRed++;
                } else if((MESSAGE_SYMBOLS[pState.at(i, j)] == 'w' || MESSAGE_SYMBOLS[pState.at(i, j)] == 'W') && (j==0 || j == 7)) {
                    nbWhite++;
                }

                /*
                // Defends
                if(i<=1 && (MESSAGE_SYMBOLS[pState.at(i, j)] == 'r'||MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')) {
                    nbRed++;
                } else if(i>=6 && (MESSAGE_SYMBOLS[pState.at(i, j)] == 'w'||MESSAGE_SYMBOLS[pState.at(i, j)] == 'W')) {
                    nbWhite++;
                }

                // Attack
                if(i<=2 && (MESSAGE_SYMBOLS[pState.at(i, j)] == 'w'||MESSAGE_SYMBOLS[pState.at(i, i)] == 'W')) {
                    nbWhite--;
                } else if(i>=5 && (MESSAGE_SYMBOLS[pState.at(i, j)] == 'r'||MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')) {
                    nbRed--;
                }
                */

                // Diagonal
                
                /*if(i<=1 && MESSAGE_SYMBOLS[pState.at(i, j)] == '.') {
                    nbRed++;
                } else if(i>=6 && MESSAGE_SYMBOLS[pState.at(i, j)] == '.') {
                    nbWhite++;
                }*/
            }
            if(MESSAGE_SYMBOLS[pState.at(i, i)] == 'r' || MESSAGE_SYMBOLS[pState.at(i, i)] == 'R'){
                nbRed++;
            } else if(MESSAGE_SYMBOLS[pState.at(i, i)] == 'w' || MESSAGE_SYMBOLS[pState.at(i, i)] == 'W') {
                nbWhite++;
            }
            if(MESSAGE_SYMBOLS[pState.at(i, 7-i)] == 'r' || MESSAGE_SYMBOLS[pState.at(i, 7-i)] == 'R'){
                nbRed++;
            } else if(MESSAGE_SYMBOLS[pState.at(i, 7-i)] == 'w' || MESSAGE_SYMBOLS[pState.at(i, 7-i)] == 'W') {
                nbWhite++;
            }
        }

        // Oreo
        if((MESSAGE_SYMBOLS[pState.at(6, 3)] == 'w' || MESSAGE_SYMBOLS[pState.at(6, 3)] == 'W') && 
            (MESSAGE_SYMBOLS[pState.at(7, 2)] == 'w' || MESSAGE_SYMBOLS[pState.at(7, 2)] == 'W') &&
            (MESSAGE_SYMBOLS[pState.at(7, 4)] == 'w' || MESSAGE_SYMBOLS[pState.at(7, 4)] == 'W')) {
            nbWhite++;
        }
        if((MESSAGE_SYMBOLS[pState.at(0, 3)] == 'r' || MESSAGE_SYMBOLS[pState.at(0, 3)] == 'R') && 
            (MESSAGE_SYMBOLS[pState.at(0, 5)] == 'r' || MESSAGE_SYMBOLS[pState.at(0, 5)] == 'R') &&
            (MESSAGE_SYMBOLS[pState.at(1, 4)] == 'r' || MESSAGE_SYMBOLS[pState.at(1, 4)] == 'R')) {
            nbRed++;
        }

        // Triangle
        if((MESSAGE_SYMBOLS[pState.at(6, 5)] == 'w' || MESSAGE_SYMBOLS[pState.at(6, 5)] == 'W') && 
            (MESSAGE_SYMBOLS[pState.at(7, 4)] == 'w' || MESSAGE_SYMBOLS[pState.at(7, 4)] == 'W') &&
            (MESSAGE_SYMBOLS[pState.at(7, 6)] == 'w' || MESSAGE_SYMBOLS[pState.at(7, 6)] == 'W')) {
            nbWhite++;
        }
        if((MESSAGE_SYMBOLS[pState.at(0, 1)] == 'r' || MESSAGE_SYMBOLS[pState.at(0, 1)] == 'R') && 
            (MESSAGE_SYMBOLS[pState.at(0, 3)] == 'r' || MESSAGE_SYMBOLS[pState.at(0, 3)] == 'R') &&
            (MESSAGE_SYMBOLS[pState.at(1, 2)] == 'r' || MESSAGE_SYMBOLS[pState.at(1, 2)] == 'R')) {
            nbRed++;
        } 
    }
    else {
        for(int i=0;i<8;i++) {
            for(int j=0;j<8;j++) {
                if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'r')
                    nbRed += 1;
                else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'R')
                    nbRed += 4;
                else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'w')
                    nbWhite += 1;
                else if(MESSAGE_SYMBOLS[pState.at(i, j)] == 'W')
                    nbWhite += 4;
            }
        }
    }

    

    
    if(player == CELL_RED)
        return (nbRed-nbWhite)/this->countPieces(pState);
    else
        return (nbWhite-nbRed)/this->countPieces(pState);
}




double Player::isInHistory(const GameState & pState, unsigned depth) {
    std::string message = this->boardToMessage(pState);

    if(this->history[depth].find(message) != this->history[depth].end()) {
        return this->history[depth][message];
    }

    return MINIMUM;
}


double Player::alphaBeta(const GameState &pState, const uint8_t & player, unsigned depth, double alpha, double beta, bool currentPlayer) {
    //IF NO DEPTH ONLY COMPUTE SIMPLE HEURISTIC
    if(depth == 0) {
        return computeHeuristic(pState, player, currentPlayer);
    }

    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

    //IF NO POSSIBLE MOVE STOP
    if (lNextStates.size() == 0)
        return computeHeuristic(pState, player, currentPlayer);

    //VAR TO ONLY TAKE FIRST sizeSolution CHILDREN
    unsigned sizeSolution = lNextStates.size();


    // VECTORS USED FOR SORTING CHILDREN
    std::vector<double> heuristics = std::vector<double> (lNextStates.size());
    std::vector<double> heuristicsIdx = std::vector<double> (lNextStates.size());

    double tmp;
    //IF PLAYER TURN
    if(pState.getNextPlayer() == player) {
        //COMPUTE SIMPLE HEURISTIC OF EACH CHILD
        for(unsigned i = 0; i < lNextStates.size(); ++i) {
            heuristics[i] = computeHeuristic(lNextStates[i], player, currentPlayer);
            heuristicsIdx[i] = i;
        }

        //SORT CHILDREN BY HEURISTIC
        auto comparator = [&heuristics](double a, double b){ return heuristics[a] > heuristics[b]; };
        std::sort(heuristicsIdx.begin(), heuristicsIdx.end(), comparator);

        //ALPHA
        double max = MINIMUM;
        for(unsigned i = 0; i < sizeSolution; ++i) {
            tmp = isInHistory(lNextStates[heuristicsIdx[i]], depth-1);
            if(tmp == MINIMUM) {
                tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
                this->history[depth-1][this->boardToMessage(lNextStates[heuristicsIdx[i]])] = tmp;
            }
            
            if(tmp > max) {
                max = tmp;
                if(max > alpha)
                    alpha = max;
            }
            if(beta <= alpha) {
                break;
            }
        }
        //std::cerr << "HEURISTIC : " << max << std::endl;
        return max;
    }
    //IF NOT PLAYER TURN
    else {
        //COMPUTE SIMPLE HEURISTIC OF EACH CHILD
        for(unsigned i = 0; i < lNextStates.size(); ++i) {
            heuristics[i] = computeHeuristic(lNextStates[i], player, currentPlayer);
            heuristicsIdx[i] = i;
        }

        //SORT CHILDREN BY HEURISTIC
        auto comparator = [&heuristics](double a, double b){ return heuristics[a] < heuristics[b]; };
        std::sort(heuristicsIdx.begin(), heuristicsIdx.end(), comparator);

        //BETA
        double min = MAXIMUM;
        for(unsigned i = 0; i < sizeSolution; ++i) {
            tmp = isInHistory(lNextStates[heuristicsIdx[i]], depth-1);
            if(tmp == MINIMUM) {
                tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
                this->history[depth-1][this->boardToMessage(lNextStates[heuristicsIdx[i]])] = tmp;
            }

            if(tmp < min) {
                min = tmp;
                if(min < beta)
                    beta = min;
            }
            if(beta <= alpha) {
                break;
            }
        }
        //std::cerr << "HEURISTIC : " << min << std::endl;
        return min;
    }
}

GameState Player::play(const GameState &pState,const Deadline &pDue)
{    
    unsigned depth = 8;

    //POSSIBLE MOVES
    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);


    if (lNextStates.size() == 0) {
        std::cerr << "NO MOVE" << std::endl;
        return GameState(pState, Move());
    }

    double alpha = double(MINIMUM);
    double beta = double(MAXIMUM);

    double max = MINIMUM;
    GameState bestMove;
    for(unsigned i = 0; i < lNextStates.size(); ++i) {
        double tmp = alphaBeta(lNextStates[i], pState.getNextPlayer(), depth, alpha, beta, false);
        if(tmp > max) {
            max = tmp;
            bestMove = lNextStates[i];
            alpha = max;
        }
    }


    return bestMove;
    //return lNextStates[rand() % lNextStates.size()];
}

/*namespace checkers*/ }

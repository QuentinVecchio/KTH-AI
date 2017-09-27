#include "player.hpp"
#include <cstdlib>
#include <limits>
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <algorithm>

#define MINIMUM std::numeric_limits<double>::min()
#define MAXIMUM std::numeric_limits<double>::max()

namespace checkers
{

Player::Player() {
    this->nbTurns = 0;
    this->history = std::vector< std::map<std::string, double> > (10);
    for(unsigned i = 0; i < this->history.size(); ++i) {
        this->history[i] = std::map<std::string, double>();
    }
    this->states[CELL_RED] = std::map<std::string, double>();
    this->states[CELL_WHITE] = std::map<std::string, double>();
}

unsigned Player::countPieces(const GameState &pState) {
    unsigned sum = 0;
    std::string grid = this->boardToString(pState);
    for(unsigned i = 0; i < grid.size(); ++i) {
        if(grid[i] == 'r' || grid[i] == 'R' || grid[i] == 'w' || grid[i] == 'W')
            ++sum;
    }
    return sum;
}

double Player::computeHeuristic(const GameState &pState, const uint8_t & player, bool currentPlayer) {
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
    std::string grid = this->boardToString(pState);

    if(!currentPlayer) {
        
        for(unsigned i = 0; i < grid.size(); ++i) {
            if(grid[i] == 'r' || grid[i] == 'R')
                ++nbRed;
            else if(grid[i] == 'w' || grid[i] == 'W')
                ++nbWhite;
        }
        if(player == CELL_RED)
            return nbRed/nbWhite;
        else
            return nbWhite/nbRed;
    } else {
        for(unsigned i = 0; i < grid.size(); ++i) {
            if(grid[i] == 'r')
                ++nbRed;
            else if(grid[i] == 'R')
                nbRed += 2;
            else if(grid[i] == 'w')
                ++nbWhite;
            else if(grid[i] == 'W')
                nbWhite += 2;
        }
        if(player == CELL_RED)
            return nbRed/nbWhite;
        else
            return nbWhite/nbRed;
    }
}

std::string Player::boardToString(const GameState & pState) const {
    std::stringstream ss;
    for(unsigned i=0;i<pState.cSquares;i++) {
        ss << SIMPLE_TEXT[pState.at(i)];
    }

    return ss.str();
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
        for(unsigned i = 0; i < heuristicsIdx.size(); ++i) {
            tmp = MINIMUM;
            /*tmp = isInHistory(lNextStates[heuristicsIdx[i]], depth-1);
            if(tmp == MINIMUM) {
                
                this->history[depth-1][this->boardToString(lNextStates[heuristicsIdx[i]])] = tmp;
            }*/
            tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
            /*if(!this->stateExist(player, lNextStates[heuristicsIdx[i]])) {
                tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
                this->createState(player, lNextStates[heuristicsIdx[i]], tmp);
            } else {
                tmp = this->getGainForState(player, lNextStates[heuristicsIdx[i]]);
            }*/
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
            tmp = MAXIMUM;
            /*tmp = isInHistory(lNextStates[heuristicsIdx[i]], depth-1);
            if(tmp == MINIMUM) {
                tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
                this->history[depth-1][this->boardToString(lNextStates[heuristicsIdx[i]])] = tmp;
            }*/
            tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
            /*if(!this->stateExist(pState.getNextPlayer(), lNextStates[heuristicsIdx[i]])) {
                tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
                this->createState(pState.getNextPlayer(), lNextStates[heuristicsIdx[i]], tmp);
            } else {
                tmp = this->getGainForState(pState.getNextPlayer(), lNextStates[heuristicsIdx[i]]);
            }*/
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

bool Player::stateExist(uint8_t player, GameState state) const {
    std::string grid = this->boardToString(state);
    return this->states.find(player)->second.find(grid) != this->states.find(player)->second.end();
}

double Player::getGainForState(uint8_t player, GameState state) const {
    std::string grid = this->boardToString(state);
    return this->states.find(player)->second.find(grid)->second;
}

void Player::createState(uint8_t player, GameState state, double value) {
    std::string grid = this->boardToString(state);
    this->states.find(player)->second.insert(std::pair<std::string, double>(grid, value));
}

/*GameState Player::play(const GameState &pState,const Deadline &pDue)
{    
    this->states[CELL_RED].clear();
    this->states[CELL_WHITE].clear();

    unsigned depth = 6;

    if(countPieces(pState) <= 4) {
        depth = 8;
    }

    //POSSIBLE MOVES
    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);


    if (lNextStates.size() == 0) {
        std::cerr << "NO MOVE" << std::endl;
        return GameState(pState, Move());
    }

    double alpha = double(MINIMUM);
    double beta = double(MAXIMUM);

    GameState bestMove;
    for(unsigned i = 0; i < lNextStates.size(); ++i) {
        double tmp = alphaBeta(lNextStates[i], pState.getNextPlayer(), depth, alpha, beta, false);
        if(tmp > alpha) {
            bestMove = lNextStates[i];
            alpha = tmp;
        }
    }


    return bestMove;
    //return lNextStates[rand() % lNextStates.size()];
}*/

GameState Player::play(const GameState &pState,const Deadline &pDue)
{    
    this->states[CELL_RED].clear();
    this->states[CELL_WHITE].clear();

    unsigned depth = 6;

    if(countPieces(pState) <= 4) {
        depth = 8;
    }

    //POSSIBLE MOVES
    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);


    if (lNextStates.size() == 0) {
        std::cerr << "NO MOVE" << std::endl;
        return GameState(pState, Move());
    }

    double alpha = double(MINIMUM);
    double beta = double(MAXIMUM);
    GameState bestMove;
    for(unsigned i = 0; i < lNextStates.size(); ++i) {
        double tmp = alphaBeta(lNextStates[i], pState.getNextPlayer(), depth, alpha, beta, false);
        if(tmp > alpha) {
            bestMove = lNextStates[i];
            alpha = tmp;
        }
    }


    return bestMove;
    //return lNextStates[rand() % lNextStates.size()];
}

/*namespace checkers*/ }

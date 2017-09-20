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

int Player::minMax(const GameState &pState, uint8_t player, unsigned depth) {
	if(depth == 0) {
		return computeHeuristic(pState, player);
	}

	uint8_t nextPlayer = pState.getNextPlayer();

	std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

    if (lNextStates.size() == 0)
    	return computeHeuristic(pState, player);

    if(nextPlayer == player) {
    	int max = MINIMUM;
	    for(unsigned i = 0; i < lNextStates.size(); ++i) {
	    	int tmp = minMax(lNextStates[i], player, depth-1);
	    	if(tmp > max) {
	    		max = tmp;
	    	}
	    }
	    return max;
    }
    else {
    	int min = MAXIMUM;
	    for(unsigned i = 0; i < lNextStates.size(); ++i) {
	    	int tmp = minMax(lNextStates[i], player, depth-1);
	    	if(tmp < min) {
	    		min = tmp;
	    	}
	    }
	    return min;
    }
}

int Player::alphaBeta(const GameState &pState, uint8_t player, unsigned depth, int alpha, int beta) {
	if(depth == 0) {
		return computeHeuristic(pState, player);
	}

	uint8_t nextPlayer = pState.getNextPlayer();

	std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

    if (lNextStates.size() == 0)
    	return computeHeuristic(pState, player);

    if(nextPlayer == player) {
    	int max = MINIMUM;
	    for(unsigned i = 0; i < lNextStates.size(); ++i) {
	    	int tmp = alphaBeta(lNextStates[i], player, depth-1, alpha, beta);
	    	if(tmp > max) {
	    		max = tmp;
	    		if(max > alpha)
	    			alpha = max;
	    	}
	    	if(beta <= alpha)
	    		break;
	    }
	    return max;
    }
    else {
    	int min = MAXIMUM;
	    for(unsigned i = 0; i < lNextStates.size(); ++i) {
	    	int tmp = alphaBeta(lNextStates[i], player, depth-1, alpha, beta);
	    	if(tmp < min) {
	    		min = tmp;
	    		if(min < beta)
	    			beta = min;
	    	}
	    	if(beta <= alpha)
	    		break;
	    }
	    return min;
    }
}


int Player::computeHeuristic(const GameState &pState, uint8_t player) {
	int tmpSum = 0;
	int i, j;
	// COUNT OUR POINTS
	int sum = 0;

	if(pState.isXWin() || pState.isOWin()) {
		if(player != pState.getNextPlayer()) {
			sum = 100;
			return sum;
		}
		else {
			sum = -100;
			return sum;
		}
	}

	// Verify each line
	for(unsigned i = 0; i < HAUTEUR; ++i) {
		int tmpSum = 0;
		for(unsigned j = 0; j < LARGEUR; ++j) {
			if(pState.at(i,j) != player && (pState.at(i,j) == CELL_X || pState.at(i,j) == CELL_O)) {
				tmpSum = 0;
				break;
			}
			if(pState.at(i,j) == player) {
				tmpSum++;
			}
		}
		sum += tmpSum;
	}

	// Verify each column
	for(unsigned i = 0; i < LARGEUR; ++i) {
		int tmpSum = 0;
		for(unsigned j = 0; j < HAUTEUR; ++j) {
			if(pState.at(j,i) != player && (pState.at(j,i) == CELL_X || pState.at(j,i) == CELL_O)) {
				tmpSum = 0;
				break;
			}
			if(pState.at(j,i) == player) {
				tmpSum++;
			}
		}
		sum += tmpSum;
	}

	// Verify diagonals
	for(i = 0, j = 0; i < HAUTEUR && j < LARGEUR; ++i, ++j) {
		if(pState.at(i,j) != player && (pState.at(i,j) == CELL_X || pState.at(i,j) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(i,j) == player) {
			tmpSum++;
		}
	}
	sum += tmpSum;

	tmpSum = 0;
	for(i = HAUTEUR-1, j = LARGEUR-1; i >= 0 && j >= 0; --i, --j) {
		if(pState.at(i,j) != player && (pState.at(i,j) == CELL_X || pState.at(i,j) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(i,j) == player) {
			tmpSum++;
		}
	}
	sum += tmpSum;

	return sum;
}

GameState Player::play(const GameState &pState,const Deadline &pDue)
{
    //std::cerr << "Processing " << pState.toMessage() << std::endl;

    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

	/*
	 * Here you should write your clever algorithms to get the best next move, ie the best
	 * next state. This skeleton returns a random move instead.
	 */
	if (lNextStates.size() == 0) {
		std::cerr << "NO MOVE" << std::endl;
		return GameState(pState, Move());
	}


    int max = MINIMUM;
    GameState bestMove;
    for(unsigned i = 0; i < lNextStates.size(); ++i) {
    	int tmp = alphaBeta(lNextStates[i], pState.getNextPlayer(), DEPTH, MINIMUM, MAXIMUM);
    	if(tmp > max) {
    		max = tmp;
    		bestMove = lNextStates[i];
    	}
    }


    return bestMove;
    //return lNextStates[rand() % lNextStates.size()];
}

/*namespace TICTACTOE*/ }

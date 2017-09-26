#include "player.hpp"
#include <cstdlib>
#include <limits>
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <algorithm>

#define MINIMUM std::numeric_limits<int>::min()
#define MAXIMUM std::numeric_limits<int>::max()

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

int Player::computeHeuristic(const GameState &pState, const uint8_t & player, bool currentPlayer) {
	if(pState.isRedWin() || pState.isWhiteWin()) {
		if(currentPlayer) {
			return 100;
		}
		else {
			return -100;
		}
	}

	unsigned nbRed = 0;
	unsigned nbWhite = 0;

	std::string grid = pState.toMessage();
	for(unsigned i = 0; i < grid.size(); ++i) {
		if(grid[i] == 'r')
			++nbRed;
		else if (grid[i] == 'R')
			nbRed += 2;
		else if(grid[i] == 'w')
			++nbWhite;
		else if (grid[i] == 'W')
			nbWhite += 2;
	}

	if(player == CELL_RED)
		return nbRed-nbWhite;
	else
		return nbWhite-nbRed;
}


int Player::isInHistory(const GameState & pState, unsigned depth) {
	std::string message = std::string(pState.toMessage(), 0, 34);

	if(this->history[depth].find(message) != this->history[depth].end()) {
		return this->history[depth][message];
	}

	return MINIMUM;
}


int Player::alphaBeta(const GameState &pState, const uint8_t & player, unsigned depth, int alpha, int beta, bool currentPlayer) {
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
    std::vector<int> heuristics = std::vector<int> (lNextStates.size());
    std::vector<int> heuristicsIdx = std::vector<int> (lNextStates.size());

    int tmp;
    //IF PLAYER TURN
    if(pState.getNextPlayer() == player) {
    	//COMPUTE SIMPLE HEURISTIC OF EACH CHILD
    	for(unsigned i = 0; i < lNextStates.size(); ++i) {
    		heuristics[i] = computeHeuristic(lNextStates[i], player, currentPlayer);
    		heuristicsIdx[i] = i;
    	}

    	//SORT CHILDREN BY HEURISTIC
    	auto comparator = [&heuristics](int a, int b){ return heuristics[a] > heuristics[b]; };
  		std::sort(heuristicsIdx.begin(), heuristicsIdx.end(), comparator);

  		//ALPHA
    	int max = MINIMUM;
	    for(unsigned i = 0; i < sizeSolution; ++i) {
	    	tmp = isInHistory(lNextStates[heuristicsIdx[i]], depth-1);
	    	if(tmp == MINIMUM) {
	    		tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
	    		this->history[depth-1][std::string(lNextStates[heuristicsIdx[i]].toMessage(), 0, 34)] = tmp;
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
    	auto comparator = [&heuristics](int a, int b){ return heuristics[a] < heuristics[b]; };
  		std::sort(heuristicsIdx.begin(), heuristicsIdx.end(), comparator);

  		//BETA
    	int min = MAXIMUM;
	    for(unsigned i = 0; i < sizeSolution; ++i) {
    		tmp = isInHistory(lNextStates[heuristicsIdx[i]], depth-1);
	    	if(tmp == MINIMUM) {
	    		tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
	    		this->history[depth-1][std::string(lNextStates[heuristicsIdx[i]].toMessage(), 0, 34)] = tmp;
	    	}

    		tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, !currentPlayer);
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

	if(countPieces(pState) <= 16) {
		depth = 6;
	}

	//POSSIBLE MOVES
    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);


	if (lNextStates.size() == 0) {
		std::cerr << "NO MOVE" << std::endl;
		return GameState(pState, Move());
	}

	int alpha = int(MINIMUM);
	int beta = int(MAXIMUM);

    int max = MINIMUM;
    GameState bestMove;
    for(unsigned i = 0; i < lNextStates.size(); ++i) {
    	int tmp = alphaBeta(lNextStates[i], pState.getNextPlayer(), depth, alpha, beta, false);
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

#include "player.hpp"
#include <cstdlib>
#include <limits>
#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <algorithm>

#define DIM 4

#define DEPTH 2

#define MINIMUM std::numeric_limits<int>::min()
#define MAXIMUM std::numeric_limits<int>::max()



namespace TICTACTOE3D
{

int evalMatrix2D(int* matrice, const uint8_t & player) {
	int tmpSum = 0;
	int i, j;
	// COUNT OUR POINTS
	int sum = 0;

	// Verify each line
	for(unsigned i = 0; i < DIM; ++i) {
		int tmpSum = 0;
		for(unsigned j = 0; j < DIM; ++j) {
			if(matrice[i*DIM+j] != player && (matrice[i*DIM+j]== CELL_X || matrice[i*DIM+j] == CELL_O)) {
				tmpSum = 0;
				break;
			}
			if(matrice[i*DIM+j] == player) {
				tmpSum++;
			}
		}
		sum += (tmpSum*10);
	}

	// Verify each column
	for(unsigned i = 0; i < DIM; ++i) {
		int tmpSum = 0;
		for(unsigned j = 0; j < DIM; ++j) {
			if(matrice[j*DIM+i] != player && (matrice[j*DIM+i] == CELL_X || matrice[j*DIM+i] == CELL_O)) {
				tmpSum = 0;
				break;
			}
			if(matrice[j*DIM+i] == player) {
				tmpSum++;
			}
		}
		sum += (tmpSum*10);
	}

	// Verify diagonals
	for(i = 0, j = 0; i < DIM && j < DIM; ++i, ++j) {
		if(matrice[i*DIM+j] != player && (matrice[i*DIM+j] == CELL_X || matrice[i*DIM+j] == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(matrice[i*DIM+j] == player) {
			tmpSum++;
		}
	}
	sum += (tmpSum*10);

	tmpSum = 0;
	for(i = DIM-1, j = DIM-1; i >= 0 && j >= 0; --i, --j) {
		if(matrice[i*DIM+j] != player && (matrice[i*DIM+j] == CELL_X || matrice[i*DIM+j] == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(matrice[i*DIM+j] == player) {
			tmpSum++;
		}
	}
	sum += (tmpSum*10);

	return sum;
}

int Player::isInHistory(const GameState & pState) {
	std::string message = std::string(pState.toMessage());

	if(this->history.find(message) != this->history.end())
		return this->history[message];

	return MINIMUM;

}

int Player::computeHeuristic(const GameState &pState, const uint8_t & player) {
	int i,j,k;

	// COUNT OUR POINTS
	int sum = 0;

	uint8_t otherPlayer = player ^ (CELL_X | CELL_O);

	if(pState.isXWin() || pState.isOWin()) {
		if(player != pState.getNextPlayer()) {
			sum = 100000;
			return sum;
		}
		else {
			sum = -100000;
			return sum;
		}
	}

    //Verify top to bottom
    for(i = 0; i < DIM; ++i) {
    	for(j = 0; j < DIM; ++j) {
    		for(k = 0; k < DIM; ++k) {
    			this->matrice[j*DIM+k] = pState.at(j,k,i);
    		}
    	}
    	sum += evalMatrix2D(this->matrice, player);
    	sum -= evalMatrix2D(this->matrice, otherPlayer);
    }


    //Verify behind to front   
	for(j = 0; j < DIM; ++j) {
		for(i = 0; i < DIM; ++i) {
    		for(k = 0; k < DIM; ++k) {
    			this->matrice[i*DIM+k] = pState.at(j,k,i);
    		}
    	}
    	sum += evalMatrix2D(this->matrice, player);
    	sum -= evalMatrix2D(this->matrice, otherPlayer);
    }


    //Verify left to right
    for(k = 0; k < DIM; ++k) {
		for(i = 0; i < DIM; ++i) {
			for(j = 0; j < DIM; ++j) {
				this->matrice[i*DIM+j] = pState.at(j,k,i);
			}
		}
    	sum += evalMatrix2D(this->matrice, player);
    	sum -= evalMatrix2D(this->matrice, otherPlayer);
    }

    //4 diagonals 3D
    unsigned tmpSum = 0;

    for(i = 0, j = 0, k = 0; i < DIM && j < DIM && k < DIM; ++i, ++j, ++k) {
		if(pState.at(j,k,i) != player && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == player) {
			tmpSum++;
		}
    }
    sum += (tmpSum*10);

    tmpSum = 0;
    for(i = 0, j = DIM, k = DIM; i < DIM && j >= 0 && k >= 0; ++i, --j, --k) {
		if(pState.at(j,k,i) != player && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == player) {
			tmpSum++;
		}
    }
    sum += (tmpSum*10);


    tmpSum = 0;
    for(i = 0, j = DIM, k = 0; i < DIM && j >= 0 && k < DIM; ++i, --j, ++k) {
		if(pState.at(j,k,i) != player && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == player) {
			tmpSum++;
		}
    }
    sum += (tmpSum*10);


    tmpSum = 0;
    for(i = 0, j = 0, k = DIM; i < DIM && j < DIM && k >= 0; ++i, ++j, --k) {
		if(pState.at(j,k,i) != player && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == player) {
			tmpSum++;
		}
    }
    sum += (tmpSum*10);


    //4 diagonals 3D OTHERPLAYER 
    tmpSum = 0;
    for(i = 0, j = 0, k = 0; i < DIM && j < DIM && k < DIM; ++i, ++j, ++k) {
		if(pState.at(j,k,i) != otherPlayer && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == otherPlayer) {
			tmpSum++;
		}
    }
    sum -= (tmpSum*10);

    tmpSum = 0;
    for(i = 0, j = DIM, k = DIM; i < DIM && j >= 0 && k >= 0; ++i, --j, --k) {
		if(pState.at(j,k,i) != otherPlayer && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == otherPlayer) {
			tmpSum++;
		}
    }
    sum -= (tmpSum*10);


    tmpSum = 0;
    for(i = 0, j = DIM, k = 0; i < DIM && j >= 0 && k < DIM; ++i, --j, ++k) {
		if(pState.at(j,k,i) != otherPlayer && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == otherPlayer) {
			tmpSum++;
		}
    }
    sum -= (tmpSum*10);


    tmpSum = 0;
    for(i = 0, j = 0, k = DIM; i < DIM && j < DIM && k >= 0; ++i, ++j, --k) {
		if(pState.at(j,k,i) != otherPlayer && (pState.at(j,k,i) == CELL_X || pState.at(j,k,i) == CELL_O)) {
			tmpSum = 0;
			break;
		}
		if(pState.at(j,k,i) == otherPlayer) {
			tmpSum++;
		}
    }
    sum -= (tmpSum*10);

	return sum;
}



int Player::alphaBeta(const GameState &pState, const uint8_t & player, unsigned depth, int alpha, int beta, const Deadline &deadline) {
	//IF NO DEPTH ONLY COMPUTE SIMPLE HEURISTIC
	if(depth == 0) {
		return computeHeuristic(pState, player);
	}

	std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

    //IF NO POSSIBLE MOVE STOP
    if (lNextStates.size() == 0 || deadline - TICTACTOE3D::Deadline::now() <= 0.2)
    	return computeHeuristic(pState, player);

    //VAR TO ONLY TAKE FIRST sizeSolution CHILDREN
    unsigned sizeSolution = lNextStates.size();
	if(sizeSolution > 40)
    	sizeSolution = 4;
    else if(sizeSolution > 30)
    	sizeSolution = 14;
    // VECTORS USED FOR SORTING CHILDREN
    std::vector<int> heuristics = std::vector<int> (lNextStates.size());
    std::vector<int> heuristicsIdx = std::vector<int> (lNextStates.size());

    int tmp;
    //IF PLAYER TURN
    if(pState.getNextPlayer() == player) {
    	//COMPUTE SIMPLE HEURISTIC OF EACH CHILD
    	for(unsigned i = 0; i < lNextStates.size(); ++i) {
    		heuristics[i] = computeHeuristic(lNextStates[i], player);
    		heuristicsIdx[i] = i;
    	}

    	//SORT CHILDREN BY HEURISTIC
    	auto comparator = [&heuristics](int a, int b){ return heuristics[a] > heuristics[b]; };
  		std::sort(heuristicsIdx.begin(), heuristicsIdx.end(), comparator);

  		//ALPHA
    	int max = MINIMUM;
	    for(unsigned i = 0; i < sizeSolution; ++i) {
	    	tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, deadline);
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
    		heuristics[i] = computeHeuristic(lNextStates[i], player);
    		heuristicsIdx[i] = i;
    	}

    	//SORT CHILDREN BY HEURISTIC
    	auto comparator = [&heuristics](int a, int b){ return heuristics[a] < heuristics[b]; };
  		std::sort(heuristicsIdx.begin(), heuristicsIdx.end(), comparator);

  		//BETA
    	int min = MAXIMUM;
	    for(unsigned i = 0; i < sizeSolution; ++i) {
    		tmp = alphaBeta(lNextStates[heuristicsIdx[i]], player, depth-1, alpha, beta, deadline);
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
	//DEPTH IS ONE FOR 6*2 = 12 TURN
	//THEN IT'S TWO
	unsigned depth = 2;
	if(this->nbTurn > 2) {
		depth = 3;
	}


	this->nbTurn++;
    

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
    	int tmp = alphaBeta(lNextStates[i], pState.getNextPlayer(), depth, alpha, beta, pDue);
    	if(tmp > max) {
    		max = tmp;
    		bestMove = lNextStates[i];
    	}
    }


    return bestMove;
    //return lNextStates[rand() % lNextStates.size()];
}

/*namespace TICTACTOE3D*/ }

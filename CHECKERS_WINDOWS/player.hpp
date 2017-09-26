#ifndef _CHECKERS_PLAYER_HPP_
#define _CHECKERS_PLAYER_HPP_

#include "constants.hpp"
#include "deadline.hpp"
#include "move.hpp"
#include "gamestate.hpp"
#include <vector>
#include <map>

namespace checkers
{

class Player
{
public:
	unsigned nbTurns;
	std::vector< std::map<std::string, int> > history;

	Player() {
        this->nbTurns = 0;
        this->history = std::vector< std::map<std::string, int> > (10);
        for(unsigned i = 0; i < this->history.size(); ++i) {
        	this->history[i] = std::map<std::string, int>();
        }
	}


	unsigned countPieces(const GameState &pState);
	int computeHeuristic(const GameState &pState, const uint8_t & player, bool currentPlayer);
	int isInHistory(const GameState & pState, unsigned depth);
	int alphaBeta(const GameState &pState, const uint8_t & player, unsigned depth, int alpha, int beta, bool currentPlayer);
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move
    GameState play(const GameState &pState, const Deadline &pDue);
};

/*namespace checkers*/ }

#endif

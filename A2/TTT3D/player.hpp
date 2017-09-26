#ifndef _TICTACTOE3D_PLAYER_HPP_
#define _TICTACTOE3D_PLAYER_HPP_

#include "constants.hpp"
#include "deadline.hpp"
#include "move.hpp"
#include "gamestate.hpp"
#include <vector>
#include <map>

#define DIM 4

namespace TICTACTOE3D
{

class Player
{
private:
	std::map<std::string, int> history;
public:
	Player() {
		this->history = std::map<std::string, int>();
        this->nbTurn = 0;
        this->matrice = new int[DIM*DIM];
	}

    unsigned nbTurn;
    int* matrice;
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move
    GameState play(const GameState &pState, const Deadline &pDue);
    int alphaBeta(const GameState &pState, const uint8_t & player, unsigned depth, int alpha, int beta, const Deadline &deadline);
    int computeHeuristic(const GameState &pState, const uint8_t & player);
	int isInHistory(const GameState & pState);
};

/*namespace TICTACTOE3D*/ }

#endif

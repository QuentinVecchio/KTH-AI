#ifndef _TICTACTOE_PLAYER_HPP_
#define _TICTACTOE_PLAYER_HPP_

#include "constants.hpp"
#include "deadline.hpp"
#include "move.hpp"
#include "gamestate.hpp"
#include <vector>
#include <stdint.h>

namespace TICTACTOE
{

class Player
{
public:
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move
    GameState play(const GameState &pState, const Deadline &pDue);

private:
	int alphaBeta(const GameState &pState, uint8_t player, unsigned depth, int alpha, int beta);
    int minMax(const GameState &pState, uint8_t player, unsigned depth);
    int computeHeuristic(const GameState &pState, uint8_t player);


};

/*namespace TICTACTOE*/ }

#endif

#ifndef _TICTACTOE_PLAYER_HPP_
#define _TICTACTOE_PLAYER_HPP_

#include "constants.hpp"
#include "deadline.hpp"
#include "move.hpp"
#include "gamestate.hpp"
#include <vector>
#include <map>

namespace TICTACTOE
{

class Player
{
private:
    std::map<std::string, int> states;
    int count;
public:
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move
    GameState play(const GameState &pState, const Deadline &pDue);

private:
    int minMax(GameState state, int player, int depth);
    int alphaBeta(GameState state, int depth, int alpha, int beta, int player);
    int eval(GameState state, int player)  const;
    bool stateExist(GameState state) const;
    void createState(GameState state, int v);
    int** matriceRotation(int** matrice, int n, int m) const;
    int getGainForState(GameState state) const;
    int gain(int nbCell, int player) const;
};

/*namespace TICTACTOE*/ }

#endif

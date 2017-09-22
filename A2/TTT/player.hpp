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
    std::map<int, std::map<std::string, int> > states;
    int count;

public:
    Player();
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move
    GameState play(const GameState &pState, const Deadline &pDue);

private:
    int minMax(GameState state, uint8_t player, int depth);
    //int alphaBeta(GameState state, int depth, int alpha, int beta, uint8_t player);
    int alphaBeta(GameState state, uint8_t player, unsigned depth, int alpha, int beta);
    int eval(GameState state, uint8_t player)  const;
    bool stateExist(unsigned player, GameState state) const;
    void createState(unsigned player, GameState state, int v);
    int** matriceRotation(int** matrice, int n, int m) const;
    int getGainForState(unsigned player, GameState state) const;
    int gain(int nbCell, bool currentPlayer) const;
    int factoriel(int n);
};

/*namespace TICTACTOE*/ }

#endif

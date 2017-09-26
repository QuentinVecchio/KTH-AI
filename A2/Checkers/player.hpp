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
private:
    unsigned nbTurns;
    std::vector< std::map<std::string, double> > history;
    std::map<int, std::map<std::string, double> > states;

    bool stateExist(uint8_t player, GameState state) const;
    void createState(uint8_t player, GameState state, double v);
    double getGainForState(uint8_t player, GameState state) const;
    std::string boardToString(const GameState & pState) const;
    unsigned countPieces(const GameState &pState);
    double computeHeuristic(const GameState &pState, const uint8_t & player, bool currentPlayer);
    double alphaBeta(const GameState &pState, const uint8_t & player, unsigned depth, double alpha, double beta, bool currentPlayer);

public:
	Player();
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move
    GameState play(const GameState &pState, const Deadline &pDue);
};

/*namespace checkers*/ }

#endif

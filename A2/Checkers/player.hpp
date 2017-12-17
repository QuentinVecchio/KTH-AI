#ifndef _CHECKERS_PLAYER_HPP_
#define _CHECKERS_PLAYER_HPP_

#include "constants.hpp"
#include "deadline.hpp"
#include "move.hpp"
#include "gamestate.hpp"
#include <vector>
#include <map>
#include <time.h>

namespace checkers
{

class Player
{
public:
    unsigned nbTurns;
    std::vector< std::map<std::string, double> > history;

    Player() {
        this->nbTurns = 0;
        this->history = std::vector< std::map<std::string, double> > (10);
        for(unsigned i = 0; i < this->history.size(); ++i) {
            this->history[i] = std::map<std::string, double>();
        }
        srand(time(NULL));
    }


    unsigned countPieces(const GameState &pState);
    std::string boardToString(const GameState & pState) const;
    std::string boardToMessage(const GameState & pState) const;
    double computeHeuristic(const GameState &pState, const uint8_t & player, bool currentPlayer);
    double isInHistory(const GameState & pState, unsigned depth);
    double alphaBeta(const GameState &pState, const uint8_t & player, unsigned depth, double alpha, double beta, bool currentPlayer);
    ///perform a move
    ///\param pState the current state of the board
    ///\param pDue time before which we must have returned
    ///\return the next state the board is in after our move
    GameState play(const GameState &pState, const Deadline &pDue);
};

/*namespace checkers*/ }

#endif

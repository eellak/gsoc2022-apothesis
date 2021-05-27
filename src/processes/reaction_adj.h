//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposotion processes.
//    Copyright (C) 2019  Nikolaos (Nikos) Cheimarios
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//============================================================================

#ifndef REACTION_ADJ_H
#define REACTION_ADJ_H

#include <iostream>
#include <string>
#include <vector>

#include "process.h"
#include "reaction.h"

using namespace std;

namespace MicroProcesses
{

class ReactionAdj: public Reaction
{
public:
    ReactionAdj();
    ~ReactionAdj();

    bool rules ( Site* );
    void perform( Site* );
    
private:
    REGISTER_PROCESS(ReactionAdj)
};    
}

#endif // REACTION_H

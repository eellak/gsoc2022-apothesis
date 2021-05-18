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

#ifndef REACTION_H
#define REACTION_H

#include <iostream>
#include <string>
#include <vector>

#include "process.h"
#include "species_new.h"

using namespace std;

namespace MicroProcesses
{

class Reaction: public Process
{
public:
    Reaction();
    ~Reaction();

    // Getters for reactants and products
    inline vector< pair < int, species_new* > > getReactants(){ return m_vpReactants; }
    inline vector< pair < int, species_new* > > getProducts(){ return m_vpProducts; }

    // Set reactants and products
    void addReactants( const double coeff, species_new* species );
    void addProducts( const double coeff, species_new* species );

    double getProbability();
    bool rules ( Site* );
    void perform( Site* );
    void print();

private:
    REGISTER_PROCESS(Reaction)
};    
}

#endif // REACTION_H

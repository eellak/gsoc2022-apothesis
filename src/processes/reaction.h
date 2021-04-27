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
    inline vector< pair < double, species_new* > > getReactants(){ return m_vpReactants; }
    inline vector< pair < double, species_new* > > getProducts(){ return m_vpProducts; }

    // Set reactants and products
    void addReactants( const double coeff, species_new* species );
    void addProducts( const double coeff, species_new* species );

    // Getters and setters for activation energy and pre-exponential factors
    inline void setActivationEnergy( double Ea) { m_dEa = Ea; }
    inline double getActivationEnergy() { return m_dEa; }
    inline void setPreExpFactor( double k0 ){ m_dK0 = k0; }
    inline double getPreExpFactor(){ return m_dK0; }

    double getProbability();
    bool rules ( Site* );
    void perform( Site* );
    void print();

private:
    /// The reactants participating in this reaction
    vector< pair< double, species_new* > > m_vpReactants;

    /// The products formed by this reaction
    vector< pair< double, species_new* > > m_vpProducts;

    /// The activation energy of this reaction
    double m_dEa;

    /// The pre-exponential factors for this reaction
    double m_dK0;

    REGISTER_PROCESS(Reaction)
};    
}

#endif // REACTION_H

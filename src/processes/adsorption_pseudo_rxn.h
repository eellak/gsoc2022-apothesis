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
#ifndef ADSORPTIONPSEUDORXN_H
#define ADSORPTIONPSEUDORXN_H

#include "process.h"
#include "species_new.h"

namespace MicroProcesses
{

class AdsorptionPseudoRxn: public Process
{
public:
    AdsorptionPseudoRxn();
    ~AdsorptionPseudoRxn() override;

    double getProbability() override;
    bool rules( Site* ) override;
    void perform( Site* ) override;

    void setAdsorptionSpecies(species_new* species) { m_aSpecies = species; }
    species_new* getAdsorptionSpecies() { return m_aSpecies; }

    void setDesorptionSpecies(species_new* species) { m_dSpecies = species; }
    species_new* getDesorptionSpecies() { return m_dSpecies; }

private:

    bool mf_isInLowerStep( Site* s );
    bool mf_isInHigherStep( Site* s );


    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The mole fraction of the AdsorptionPseudoRxn process
    double m_dMolFrac;

    ///The species that must adsopt
    species_new* m_aSpecies;

    ///The species that must desorb
    species_new* m_dSpecies;

    /// A member function to calculate the neighbors of a given site
    int mf_calculateNeighbors(Site*);

    REGISTER_PROCESS(AdsorptionPseudoRxn)
};
}

#endif // AdsorptionPseudoRxn_H

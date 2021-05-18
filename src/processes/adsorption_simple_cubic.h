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
#ifndef ADSORPTIONSIMPLECUBIC_H
#define ADSORPTIONSIMPLECUBIC_H

#include "process.h"

namespace MicroProcesses
{

class AdsorptionSimpleCubic: public Process
{
public:
    AdsorptionSimpleCubic();
    ~AdsorptionSimpleCubic() override;

    double getProbability() override;
    bool rules( Site* ) override;
    void perform( Site* ) override;

private:

    bool mf_isInLowerStep( Site* s );
    bool mf_isInHigherStep( Site* s );

    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    ///The mole fraction of the AdsorptionSimpleCubic process
    double m_dMolFrac;

    ///The species that must adsopt
//    species_new* m_Species;

    /// A member function to calculate the neighbors of a given site
    int mf_calculateNeighbors(Site*);

    REGISTER_PROCESS(AdsorptionSimpleCubic)
};
}

#endif // AdsorptionSimpleCubic_H

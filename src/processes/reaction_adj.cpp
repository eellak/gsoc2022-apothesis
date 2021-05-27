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

#include "reaction_adj.h"

namespace MicroProcesses
{

    REGISTER_PROCESS_IMPL ( ReactionAdj )

    ReactionAdj::ReactionAdj() {}
    ReactionAdj::~ReactionAdj() {}

    bool ReactionAdj::rules(Site* site)
    {
        for (pair<int, species_new*> reactant:m_vpReactants)
        {
            // If the species map contains enough of the reactant
            if (site->getSpeciesMap()[reactant.second->getID()] >= reactant.first)
            {
                continue;
            }
            else
            {
                return false;
            }
        }
        return true;
    }

    void ReactionAdj::perform(vector<Site*> sites) // This will perform when a vector is provided (rather than single site)
    {
        for (pair<int, species_new*> reactant:m_vpReactants)
        {
            // Keep track of the number of species needed to be removed from each site.
            // Assume each site can only hold 1 of an individual species
            int stoichCoeff = reactant.first;
            for (Site* site:sites)
            {
                // If the species map contains enough of the reactant
                if (site->getSpeciesMap()[reactant.first] > 1 && stoichCoeff > 0)
                {
                    site->removeSpecies(reactant.second, reactant.first);
                    // decrement stoichiometric coefficient
                    stoichCoeff--;
                }
            }
        }
        
        if (getApothesis()->getCaseStudy() == 1)
        {
            for (Site* site:sites)
            {
                site->increaseHeight(1);
            }
        }
        else if (getApothesis()->getCaseStudy() == 2)
        {
            for (Site* site:sites)
            {
                for (pair<int, species_new*> product:m_vpProducts)
                {
                    site->addSpecies(product.second, product.first);   
                }
            }
        }       
    }

}

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
            int stoichiometry = reactant.first;
            // if site has at least one
            if (site->getSpeciesMap()[reactant.second->getID()] == 1)
            {
                stoichiometry--;
                // what are the levels?
                for (Site* s:site->get1stNeighbors()[0])
                {
                    if (s->getSpeciesMap()[reactant.second->getID()] > 0)
                        stoichiometry--;
                    if (stoichiometry < 1)
                        return true;
                }
                if (stoichiometry > 1)
                    return false;
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
        // Add products if we have any
        if (m_vpProducts.size() > 0)
        {
            for (Site* site:sites)
            {
                for (pair<int, species_new*> product:m_vpProducts)
                {
                    site->addSpecies(product.second, product.first);   
                }
            }
        }
        else // Else simply increase the height
        {
            
            for (Site* site:sites)
            {
                site->increaseHeight(1);
            }
        }
        
        for (Site* site:sites)
        {
            m_seAffectedSites.insert(site);
            for ( Site* neigh:site->getNeighs() ) 
            {
                mf_calculateNeighbors( neigh );
                m_seAffectedSites.insert( neigh );
            }

        } 
    }

    int ReactionAdj::mf_calculateNeighbors(Site* s)
    {
        int neighs = 1;
        for ( Site* neigh:s->getNeighs() ) {
            if ( s->isLowerStep() && neigh->isHigherStep() ){
                if ( neigh->getHeight() >= s->getHeight() + m_pLattice->getStepDiff() + 1 )
                    neighs++;
            }
            else if ( neigh->isLowerStep() && s->isHigherStep() ){
                if ( neigh->getHeight() >= s->getHeight() - m_pLattice->getStepDiff() + 1 )
                    neighs++;
            }
            else {
                if ( neigh->getHeight() >= s->getHeight() )
                    neighs++;
            }
        }

        s->setNeighsNum( neighs );
        return neighs;
    }

}

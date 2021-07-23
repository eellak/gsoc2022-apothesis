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
#include "desorption_pseudo_rxn.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL(DesorptionPseudoRxn);

DesorptionPseudoRxn::DesorptionPseudoRxn():m_iNeigh(0){}
DesorptionPseudoRxn::~DesorptionPseudoRxn(){}


bool DesorptionPseudoRxn::rules( Site* s)
{
    if (s->getSpeciesVec().size() == 1)
    {
        if (s->getSpeciesVec()[0]->getChemFormula() != getSpecies()->getChemFormula())
        {
            return false;
        }
        else if ( mf_calculateNeighbors( s ) == any_cast<int>(m_mParams["neigh"] ) )
        {
            return true;
        }
    }
    
    return false;
}

void DesorptionPseudoRxn::perform( Site* s)
{
    //For PVD results
    s->removeSpecies(getSpecies(), 1);
    
    mf_calculateNeighbors( s ) ;
    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        mf_calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( mf_calculateNeighbors( firstNeigh ) );
            m_seAffectedSites.insert( firstNeigh );
        }
    }
}

species_new* DesorptionPseudoRxn::getLatticeSpecies(Site* site)
{
    return site->getLattice()->getLatticeSpecies();
}

int DesorptionPseudoRxn::mf_calculateNeighbors(Site* s)
{

    //We do not need to count the neighbours here!!!
    //We need it only in the rules!
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

    //For flat surfaces
/*    int neighs = 1;
    if ( mf_isInLowerStep( s ) ){
        for ( Site* neigh:s->getNeighs() ) {
            if ( mf_isInHigherStep( neigh ) ){
                if ( neigh->getHeight() >= s->getHeight() + m_pLattice->getStepDiff() + 1 )
                    neighs++;
            }
            else{
                if ( neigh->getHeight() >= s->getHeight() )
                    neighs++;
            }
        }
    }
    else if ( mf_isInHigherStep( s ) ){
        for ( Site* neigh:s->getNeighs() ) {
            if ( mf_isInLowerStep( neigh ) ){
                if ( neigh->getHeight() >= s->getHeight() - (m_pLattice->getStepDiff() + 1 ) )
                    neighs++;
            }
            else{
                if ( neigh->getHeight() >= s->getHeight() )
                    neighs++;
            }
        }
    }
    else {
        for ( Site* neigh:s->getNeighs() ) {
            if ( neigh->getHeight() >= s->getHeight() )
                neighs++;
        }
    }

    s->setNeighsNum( neighs );

    return neighs; */

  //For flat surfaces
/*    int neighs = 1;
    for ( Site* neigh:s->getNeighs() ) {
        if ( neigh->getHeight() >= s->getHeight() )
            neighs++;
    }
    return neighs;*/
}

bool DesorptionPseudoRxn::mf_isInLowerStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, 0 )->getID() )
            return true;

    return false;
}

bool DesorptionPseudoRxn::mf_isInHigherStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++){
   //     cout<< m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() << endl;
        if ( s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() ){
            return true;
        }
    }

    return false;
}

double DesorptionPseudoRxn::getProbability(){

    //These must trenafered in the global definitions
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double E_d = any_cast<double>(m_mParams["E_d"]);			// [j]
    double E = any_cast<double>(m_mParams["E_m"]); //any_cast<double>(m_mParams["E"]); //71128/Na;   //(7.14e+4)/Na;			// [j] -> 17 kcal
    double v0 = any_cast<double>(m_mParams["Freq"]);				// [s^-1]
    /*--------------------------------------------------*/

    double probability = v0*exp(-(double)any_cast<int>(m_mParams["neigh"])*E/(k*T));
    return probability;			//DesorptionPseudoRxn 1 neigh
}

}

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

#include "BCC.h"
#include "read.h"

#include <map>

BCC::BCC(Apothesis *apothesis) : Lattice(apothesis), m_iMinNeigs(1)
{
	;
}

BCC::BCC(Apothesis *apothesis, bool step, vector<int> stepInfo) : Lattice(apothesis),
																  m_hasSteps(step),
																  m_stepInfo(stepInfo)
{
	;
}

void BCC::setInitialHeight(int height) { m_iHeight = height; }

void BCC::build()
{
	if (m_Type == NONE)
	{
		cout << "Not supported lattice type" << endl;
		EXIT;
	}

	if (m_iSizeX == 0 || m_iSizeY == 0)
	{
		m_errorHandler->error_simple_msg("The lattice size cannot be zero in either dimension.");
		EXIT;
	}

	if (m_iHeight < 5)
	{
		m_errorHandler->warningSimple_msg("The lattice initial height is too small.Consider revising.");
	}

	// The sites of the lattice.
	m_vSites.resize(getSize());
	for (int i = 0; i < m_vSites.size(); i++)
		m_vSites[i] = new Site();

	//  m_pSites = new Site[ m_iSizeX*m_iSizeY];

	for (int i = 0; i < m_iSizeX; i++)
	{
		for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
		{
			m_vSites[j]->setID(j);
			m_vSites[j]->setHeight(m_iHeight - 1);
		}
	}

	for (int i = 0; i < m_iSizeX; i++)
	{
		for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
			cout << m_vSites[j]->getHeight() << " ";

		cout << " " << endl;
	}

	for (int i = 0; i < m_iSizeX; i++)
	{
		for (int j = i * m_iSizeY; j < (m_iSizeY + i * m_iSizeY); j++)
			cout << m_vSites[j]->getID() << " ";

		cout << " " << endl;
	}

	if (m_bHasSteps)
		mf_buildSteps();

	mf_neigh();
}

BCC::~BCC()
{
	for (int i = 0; i < getSize(); i++)
		delete m_vSites[i];
}

void BCC::setSteps(bool hasSteps)
{
	m_bHasSteps = hasSteps;
}

void BCC::setStepInfo(int sizeX, int sizeY, int sizeZ)
{
	m_iStepX = sizeX;
	m_iStepY = sizeY;
	m_iStepZ = sizeZ;
}

void BCC::mf_buildSteps()
{
	// Pick dimension of step
	// Largest value is the direction of stepping
	// i.e. lattice is [40 x 20], having a step array of [40, 2, 1] 
	// steps of height 1 in the z direction, and 2 consecutive sites 
	// in the y direction are the same height

	// 1. Find the initial height from arbitrary site
	int initialHeight = m_vSites[0]->getHeight();
	vector<int> currentDimensions{m_iSizeX, m_iSizeY, initialHeight};

	int iteration = 0;
	int stepDimension = 0, stepSoFar = 0, growthDimension = 0, growthSoFar = 0, latentDimension = 0;
	for (auto &dim : m_stepInfo)
	{
		// 2. If the step information is the same value as lattice dim, this will not be the step/growth dimension
		if (dim != currentDimensions[iteration])
		{
			// The step dimension will be the largest value
			if (dim > stepSoFar)
			{
				stepDimension = iteration;
				stepSoFar = dim;
			}
			else
			{
				growthDimension = iteration;
			}
		}
		else
		{
			latentDimension = iteration;
		}

		iteration++;
	}

	for (unsigned int firstDim = 0; firstDim < currentDimensions[latentDimension]; ++firstDim)
	{
		for (unsigned int secondDim = 0; secondDim < currentDimensions[stepDimension]; ++secondDim)
		{
			// Calculate how much we increase the step by.
			// Calculation is split up to ensure we have proper integer division in the first step.
			int growth = secondDim / m_stepInfo[stepDimension];
			growth *= m_stepInfo[growthDimension];
			int index = firstDim * currentDimensions[stepDimension] + secondDim;
			m_vSites[index]->increaseHeight(growth);
		}
	}
}

void BCC::mf_neigh()
{
	/* All except the boundaries */
	for (int i = 1; i < m_iSizeY - 1; i++)
	{
		for (int j = 1; j < m_iSizeX - 1; j++)
		{
			m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i - 1) * m_iSizeX + j]);
			m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[(i + 1) * m_iSizeX + j]);
			m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j + 1]);
			m_vSites[i * m_iSizeX + j]->setNeigh(m_vSites[i * m_iSizeX + j - 1]);
		}
	}

	int iFirstCorner = 0;
	int iSecondCorner = m_iSizeX - 1;
	int iThirdCorner = m_iSizeX * m_iSizeY - m_iSizeX;
	int iForthCorner = m_iSizeX * m_iSizeY - 1;

	/*First row */
	for (int j = iFirstCorner; j <= iSecondCorner; j++)
	{
		if (j != 0 && j != m_iSizeX - 1)
		{
			m_vSites[j]->setNeigh(m_vSites[j - 1]);
			m_vSites[j]->setNeigh(m_vSites[j + 1]);
			m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
			m_vSites[j]->setNeigh(m_vSites[iThirdCorner + j]);
		}
		else if (j == iFirstCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
			m_vSites[j]->setNeigh(m_vSites[1]);
			m_vSites[j]->setNeigh(m_vSites[iSecondCorner + 1]);
			m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
		}
		else if (j == iSecondCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[j - 1]);
			m_vSites[j]->setNeigh(m_vSites[0]);
			m_vSites[j]->setNeigh(m_vSites[2 * m_iSizeX - 1]);
			m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
		}
	}

	/*Last row */
	int iPos = 1;
	for (int j = iThirdCorner; j <= iForthCorner; j++)
	{
		if (j != iThirdCorner && j != iForthCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[j - 1]);
			m_vSites[j]->setNeigh(m_vSites[j + 1]);
			m_vSites[j]->setNeigh(m_vSites[iFirstCorner + iPos]);
			m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
			iPos++;
		}
		else if (j == iThirdCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[iForthCorner]);
			m_vSites[j]->setNeigh(m_vSites[iThirdCorner + 1]);
			m_vSites[j]->setNeigh(m_vSites[iFirstCorner]);
			m_vSites[j]->setNeigh(m_vSites[iThirdCorner - m_iSizeX]);
		}
		else if (j == iForthCorner)
		{
			m_vSites[j]->setNeigh(m_vSites[iForthCorner - 1]);
			m_vSites[j]->setNeigh(m_vSites[iThirdCorner]);
			m_vSites[j]->setNeigh(m_vSites[iSecondCorner]);
			m_vSites[j]->setNeigh(m_vSites[iThirdCorner - 1]);
		}
	}

	/* First column */
	for (int j = iFirstCorner + m_iSizeX; j < iThirdCorner; j += m_iSizeX)
	{
		m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX - 1]);
		m_vSites[j]->setNeigh(m_vSites[j + 1]);
		m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
		m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
	}

	/* Last column */
	for (int j = iSecondCorner + m_iSizeX; j < iForthCorner; j += m_iSizeX)
	{
		m_vSites[j]->setNeigh(m_vSites[j - 1]);
		m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX + 1]);
		m_vSites[j]->setNeigh(m_vSites[j + m_iSizeX]);
		m_vSites[j]->setNeigh(m_vSites[j - m_iSizeX]);
	}

	/*	int iCount = 0;
	int pos = 0;
	while (iCount < 100) {
		cout << "Enter pos to print neighbours: ";
		cin >> pos;
		cout << m_vSites[pos]->getID() << ": " << endl;
//		for (int i = 0; i < 4; i++) {
			cout << "WEST: " << m_vSites[pos]->getNeighPosition( Site::WEST )->getID()  << endl;
			cout << "EAST: " << m_vSites[pos]->getNeighPosition(Site::EAST)->getID() << endl;
			cout << "NORTH: " << m_vSites[pos]->getNeighPosition(Site::NORTH)->getID() << endl;
			cout << "SOUTH: " << m_vSites[pos]->getNeighPosition(Site::SOUTH)->getID() << endl;
			//		}
	}*/
}

//Site* BCC::getSite(int id) { return m_vSites[ id ]; }

void BCC::check()
{
	cout << "Checking lattice..." << endl;

	int test = 2;
	cout << test << ": ";
	cout << "W:" << getSite(test)->getNeighPosition(Site::WEST)->getID() << " ";
	cout << "E:" << getSite(test)->getNeighPosition(Site::EAST)->getID() << " ";
	cout << "N:" << getSite(test)->getNeighPosition(Site::NORTH)->getID() << " ";
	cout << "S:" << getSite(test)->getNeighPosition(Site::SOUTH)->getID() << endl;
}

int BCC::calculateNeighNum(int id)
{
	int neighs = 1;
	for (Site *s : m_vSites[id]->getNeighs())
	{
		if (s->getHeight() >= m_vSites[id]->getHeight())
			neighs++;
	}

	return neighs;
}

// We have to define the processes e.g. Adsorption Simple, Asorption MultipleSites, Diffusion 1s Neighbors,
// Desorption 1st Neighbors, Reaction Decomposition etc...

void BCC::adsorp(int siteID, species_new *chemSpec)
{
	// >--------  For Lam & Vlachos (2000) ------------------------------------//
	//Remove site and its neihbors from its previous position in diffusion and desorption classes
	for (auto &p : *m_pProcMap)
	{
		string adsorption = "Adsorption_" + chemSpec->getChemFormula();
		if (!IO::contains(p.first, "Adsorption_" + adsorption))
		{
			p.second.erase(siteID);

			for (Site *s : m_vSites[siteID]->getNeighs())
			{
				p.second.erase(s->getID());
				for (Site *firstNeigh : s->getNeighs())
					p.second.erase(firstNeigh->getID());
			}
		}
	}

	//For PVD results
	m_vSites[siteID]->increaseHeight();
	m_iSiteNeighsNum = calculateNeighNum(siteID);

	
	string strProc = "Desorption_" + chemSpec->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N" ;
	map<string, set<int>>:: iterator itr = (*m_pProcMap).find(strProc);
	if (itr != (*m_pProcMap).end())
	{
		m_pProcMap->at(strProc).insert(siteID);
	}
	
	strProc = "Diffusion_" + chemSpec->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
	itr = (*m_pProcMap).find(strProc);
	if (itr != (*m_pProcMap).end())
	{
		m_pProcMap->at(strProc).insert(siteID);
	}

	for (Site *s : m_vSites[siteID]->getNeighs())
	{
		m_iSiteNeighsNum = calculateNeighNum(s->getID());
		string strProc = "Desorption_" + chemSpec->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
		map<string, set<int>>:: iterator itr = (*m_pProcMap).find(strProc);
		if (itr != (*m_pProcMap).end())
		{
			m_pProcMap->at(strProc).insert(s->getID());
		}
		strProc = "Diffusion_" + chemSpec->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
		itr = (*m_pProcMap).find(strProc);
		if (itr != (*m_pProcMap).end())
		{
			m_pProcMap->at(strProc).insert(s->getID());
		}

		for (Site *firstNeigh : s->getNeighs())
		{
			m_iSiteNeighsNum = calculateNeighNum(firstNeigh->getID());
			string strProc = "Desorption_" + chemSpec->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
			itr = (*m_pProcMap).find(strProc);
			if (itr != (*m_pProcMap).end())
			{
				m_pProcMap->at(strProc).insert(firstNeigh->getID());
			}
			strProc = "Diffusion_" + chemSpec->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
			itr = (*m_pProcMap).find(strProc);
			if (itr != (*m_pProcMap).end())
			{
				m_pProcMap->at(strProc).insert(firstNeigh->getID());
			}
		}
	}
	// ---------  For Lam & Vlachos (2000) ------------------------------------<//

	//1.Add to this site species the new adroped species
	//m_vSites[ siteID ]->getReactSpecies()[ chemSpec->getID() ]++;

	//2.Check if this site can react (according to the reactio matrix) and put it in the appropriate process map classes.
	//  2a. If a reaction with the max coeff has been completed. Remove it from adsoprtion and wait to react.

	//3.Check if this site (according to the reactio matrix) can accept diffused species from the neighbour sites or be a donor.
	//  3a. Do that for the neighbour sites and remove/add them accordingly in the process map.
}

void BCC::desorp(int siteID, species_new *chemSpecies)
{
	// <--------  For Lam & Vlachos (2000) ------------------------------------//
	//Remove site and its neihbors from its previous position in diffusion and desorption classes
	for (auto &p : *m_pProcMap)
	{
		if (!IO::contains(p.first, "Adsorption"))
		{
			p.second.erase(siteID);

			for (Site *s : m_vSites[siteID]->getNeighs())
			{
				p.second.erase(s->getID());
				for (Site *firstNeigh : s->getNeighs())
					p.second.erase(firstNeigh->getID());
			}
		}
	}

	//For PVD results
	m_vSites[siteID]->decreaseHeight();
	m_iSiteNeighsNum = calculateNeighNum(siteID);

	string strProc = "Desorption_" + chemSpecies->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
	m_pProcMap->at(strProc).insert(siteID);

	strProc = "Diffusion_" + chemSpecies->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
	m_pProcMap->at(strProc).insert(siteID);

	for (Site *s : m_vSites[siteID]->getNeighs())
	{
		m_iSiteNeighsNum = calculateNeighNum(s->getID());
		string strProc = "Desorption_" + chemSpecies->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
		m_pProcMap->at(strProc).insert(s->getID());
		strProc = "Diffusion_" + chemSpecies->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
		m_pProcMap->at(strProc).insert(s->getID());

		for (Site *firstNeigh : s->getNeighs())
		{
			m_iSiteNeighsNum = calculateNeighNum(firstNeigh->getID());
			string strProc = "Desorption_" + chemSpecies->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
			m_pProcMap->at(strProc).insert(firstNeigh->getID());
			strProc = "Diffusion_" + chemSpecies->getChemFormula() + " " + to_string(m_iSiteNeighsNum) + "N";
			m_pProcMap->at(strProc).insert(firstNeigh->getID());
		}
	}
	// ---------  For Lam & Vlachos (2000) ------------------------------------>//

	// <--------  For Lam & Vlachos (2000) ------------------------------------//
	//Remove site from its previous positin in diffusion and desorption classes
	/*    for ( auto &p:*m_pProcMap ) {
        if ( !IO::contains( p.first, "Adsoprtion" ) )
            p.second.erase( siteID );
    }

    this->print();

    //For PVD results
    m_vSites[ siteID ]->decreaseHeight();

    //This is for the basic site
    int neighs = 1;
    for ( Site* s:m_vSites[ siteID ]->getNeighs() ){
        if ( s->getHeight() <= m_vSites[ siteID ]->getHeight() )
            m_iMinNeigs++;
        else {
            for ( auto &p:*m_pProcMap ) {
                if ( !IO::contains( p.first, "Adsoprtion" ) )
                    p.second.erase( s->getID() );
            }

            s->decreaseNeighsNum();

            string strProcN = "Desorption " + to_string( s->getNeighboursNum() ) + "N";
            m_pProcMap->at( strProcN ).insert( s->getID() );

            strProcN = "Diffusion " + to_string( s->getNeighboursNum() ) + "N";
            m_pProcMap->at( strProcN ).insert( s->getID() );
        }
    }
    m_vSites[ siteID ]->setNeighboursNum( m_iMinNeigs );

    string strProc = "Desorption " + to_string( m_iMinNeigs ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );

    strProc = "Diffusion " + to_string( m_iMinNeigs ) + "N";
    m_pProcMap->at( strProc ).insert( siteID );*/

	// ---------  For Lam & Vlachos (2000) ------------------------------------>//
}

void BCC::react(int siteID)
{
	//1.React means increase the sites height by one.

	//2.Remove every species from this site (that means that the site is occupied by a lattice site).

	//3.Update neighbours and process map according to the new height as we would do in PVD.
}

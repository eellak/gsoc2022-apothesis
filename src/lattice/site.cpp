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

#ifndef SITE_CPP
#define SITE_CPP

#include "site.h"

namespace SurfaceTiles
{

  Site::Site():m_phantom(false),m_isLowerStep(false), m_isHigherStep(false)
  {
      vector<Site* > vec;
      m_m1stNeighs = { {-1, vec}, { 0, vec }, {1, vec }, };
  }

  Site::~Site() {}

  void Site::addSpecies(species_new *s, int stoich)
  {
    for (int num = 0; num < stoich; ++num)
    {
      m_vSpecies.push_back(s);
    }
    m_mapSpecies[s->getID()]+=stoich;
  }

  void Site::removeSpecies(species_new *s, int stoich)
  {
    // TODO: Is this better to do by species pointer or by string?. I think pointer is easier
    // Search through species list until we've found the specific one to remove

    // Store previous state to ensure that we remove the desired species
    int numIter = 0;
    int m_species_prevSize = m_vSpecies.size();

    for (int num = 0; num < stoich; ++num)
    {
      int numOfSpecies = m_mapSpecies[s->getID()];

      if (numOfSpecies > 0)
      {
        for (vector<species_new *>::iterator itr = m_vSpecies.begin(); itr != m_vSpecies.end(); ++itr)
        {
          if (*itr == s)
          {
            m_vSpecies.erase(itr);
            break;
          }
          ++numIter;
        }

      }
      
      // Output warning message if we didn't remove anything
      if (numOfSpecies < 1)
      {
        cout << "Warning: did not find an instance of " << s->getChemFormula() << "in site " << getID() << endl;
      }
    }
    // Decrement number of said species
    m_mapSpecies[s->getID()]-=stoich;
  }

  vector<string> Site::getSpeciesName()
  {
    vector<string> names;
    for (vector<species_new *>::iterator itr = m_vSpecies.begin(); itr != m_vSpecies.end(); ++itr)
    {
      names.push_back((*itr)->getChemFormula());
    }
    return names;
  }

} // namespace SurfaceTiles

#endif

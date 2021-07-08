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

  void Site::addSpecies(Species *s)
  {
    m_species.push_back(s);
    m_mapSpecies[s->getId()]++;
  }

  void Site::removeSpecies(Species *s)
  {
    // TODO: Is this better to do by species pointer or by string?. I think pointer is easier
    // Search through species list until we've found the specific one to remove

    // Store previous state to ensure that we remove the desired species
    int numIter = 0;
    int m_species_prevSize = m_species.size();

    int numOfSpecies = m_mapSpecies[s->getId()];

    if (numOfSpecies > 0)
    {
      for (vector<Species *>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
      {
        if (*itr == s)
        {
          m_species.erase(itr);
          break;
        }
        ++numIter;
      }
      // Decrement number of said species
      m_mapSpecies[s->getId()]--;
    }

    // Output warning message if we didn't remove anything
    if (numOfSpecies < 1)
    {
      cout << "Warning: did not find an instance of " << s->getName() << "in site " << getID() << endl;
    }
  }

  vector<Species *> Site::getSpecies()
  {
    return m_species;
  }

  vector<string> Site::getSpeciesName()
  {
    vector<string> names;
    for (vector<Species *>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
    {
      names.push_back((*itr)->getName());
    }
    return names;
  }

} // namespace SurfaceTiles

#endif

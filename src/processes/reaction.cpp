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

#include "reaction.h"
//#include "species_new.cpp"

namespace MicroProcesses
{

    REGISTER_PROCESS_IMPL ( Reaction )

    Reaction::Reaction() {}
    Reaction::~Reaction() {}

    // move to process.h class
    void Reaction::addReactants(const double coeff, species_new *species)
    {
        pair<int, species_new *> data;
        data.first = coeff;
        data.second = species;
        m_vpReactants.push_back(data);
    }

    void Reaction::addProducts(const double coeff, species_new *species)
    {
        pair<int, species_new *> data;
        data.first = coeff;
        data.second = species;
        m_vpProducts.push_back(data);
    }

    double Reaction::getProbability()
    {
        // Initialize class once
        // Rules for case # 1 and # 2 are the same
        // For all species in the reaction
        //auto v_reactants = (vector<pair<int, species_new*>>) m_mParams["reactants"];
        double T = any_cast<double>(getParameter("T"));
        double R = any_cast<double>(getParameter("R"));
        return m_dK0*exp(-m_dActNrg/T/R); // Need the number of sites that can perform a given reaction
    }

    bool Reaction::rules(Site* site)
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

    void Reaction::perform(Site* s)
    {
        cout<<"Reacting"<<endl;
        for (pair<int, species_new*> reactant:m_vpReactants)
        {
            // If the species map contains enough of the reactant
            s->removeSpecies(reactant.second, reactant.first);
        }
        
        if (getApothesis()->getCaseStudy() == 1)
        {
            s->increaseHeight(1);
        }
        else if (getApothesis()->getCaseStudy() == 2)
        {
            for (pair<int, species_new*> product:m_vpProducts)
            {
                s->addSpecies(product.second, product.first);   
            }
        }       
    }

    void Reaction::print()
    {
        int iCount = 0;
        for (pair<int, species_new *> &p : m_vpReactants)
        {
            if (iCount != m_vpReactants.size() - 1)
                if (p.first != 1)
                    cout << p.first << " " << p.second->getChemFormula() << " + ";
                else
                    cout << p.second->getChemFormula() << " + ";
            else if (p.first != 1)
                cout << p.first << " " << p.second->getChemFormula();
            else
                cout << p.second->getChemFormula();

            iCount++;
        }

        cout << " = ";

        iCount = 0;
        for (pair<int, species_new *> &p : m_vpProducts)
        {
            if (iCount != m_vpProducts.size() - 1)
                if (p.first != 1)
                    cout << p.first << " " << p.second->getChemFormula() << " + ";
                else
                    cout << p.second->getChemFormula() << " + ";
            else if (p.first != 1)
                cout << p.first << " " << p.second->getChemFormula();
            else
                cout << p.second->getChemFormula();

            iCount++;
        }

        cout << endl;
    }

}

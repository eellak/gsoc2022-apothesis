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

#include "adsorption.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( Adsorption )

Adsorption::Adsorption():m_Species(0){}

Adsorption::~Adsorption(){}

bool Adsorption::rules( Site* s )
{
    //You can always adsorb as is
    return true;
}

void Adsorption::perform( Site* s )
{
//    m_pLattice->adsorp( siteID, m_Species );
}

//--------------------- Transitions probabilities ----------------------------------------//
/*        (*prob)[0] = pa*Nx*Ny;									//Adsorption
        (*prob)[1] = group[0].size()*v0*exp(-1.0e0*E/(k*T));			//Desorption 1 neigh
        (*prob)[3] = group[1].size()*v0*exp(-2.0e0*E/(k*T));        //Desorption 2 neigh
        (*prob)[5] = group[2].size()*v0*exp(-3.0e0*E/(k*T));		//Desorption 3 neigh
        (*prob)[7] = group[3].size()*v0*exp(-4.0e0*E/(k*T));		//Desorption 4 neigh
        (*prob)[9] = group[4].size()*v0*exp(-5.0e0*E/(k*T));		//Desorption 5 neigh
        (*prob)[2] = A*(*prob)[1];	  							//Diffusion  1 neisgh
        (*prob)[4] = A*(*prob)[3];								//Diffusion  2 neisgh
        (*prob)[6] = A*(*prob)[5];								//Diffusion  3 neisgh
        (*prob)[8] = A*(*prob)[7];								//Diffusion  4 neisgh
        (*prob)[10] = A*(*prob)[9];								//Diffusion  5 neisgh */
//----------------------------------------------------------------------------------------//

/*        (*p_tot) = 0;
        for (unsigned int i=0; i<nof_trans_prob; i++)
                 *p_tot += (*prob)[i];

        for (unsigned int i=0; i<nof_trans_prob; i++)
                cout<< "NO"<< "\t" <<(*prob)[i] << endl;

        cout<<"***************"<< endl;*/

//---------Canonical-form--------------//
//       for (unsigned int i=0; i<nof_trans_prob; i++)
//               (*prob)[i] /= (*p_tot);
//-------------------------------------//

/*	for (unsigned int i=0; i<nof_trans_prob; i++)
                cout<<"CAN"<< "\t" <<(*prob)[i] << endl;
        cout<<"***************"<< endl;
        cout<<"PROBS"<<endl;
        system("pause");
--------------------------------------------------*/

double Adsorption::getProbability(){

    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = 500;						// [K]
    double k = 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = 0.1;
    double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
    double E_d = (7.14e+4)/Na;			// [j]
    double E = 71128/Na;   //(7.14e+4)/Na;			// [j]
    double m = 32e-3/Na;				// [kg]
    double E_m = (4.28e+4)/Na;			// [j]
    double k_d = 1.0e+13;				// [s^-1]
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafer

    return s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}

}

#include "adsorption_new.h"

Adsorption_new::Adsorption_new() : m_Species(0) {}

Adsorption_new::Adsorption_new(double actNrg, double molFrac, double sticking,
                               double siteDensity, double mass, species_new *species,
                               Parameters *parameter)
    : m_dActNrg(actNrg),
      m_dMolFrac(molFrac),
      m_dStick(sticking),
      m_siteDensity(siteDensity),
      m_mass(mass),
      m_Species(species),
      pParameters(parameter)
{
}

Adsorption_new::~Adsorption_new() {}

void Adsorption_new::perform(int siteID)
{
        m_pLattice->adsorp(siteID, m_Species);
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

double Adsorption_new::getProbability()
{

        //Constants
        double Na = 6.0221417930e+23; // Avogadro's number [1/mol]
        double k = 1.3806503e-23;     // Boltzmann's constant [j/K]

        //Case-dependent values read from input file
        double P = pParameters->getPressure();    // [Pa]
        double T = pParameters->getTemperature(); // [K]
        double s0 = getSticking();
        double C_tot = getSiteDensity(); // [sites/m^2] Vlachos code says [moles sites/m^2]
        double m = getMass() / Na;           // [kg]
        double y = getMolFrac();         // Mole fraction of the precursor on the wafer

        return s0 * y * P / (C_tot * sqrt(2.0e0 * 3.14159265 * m * k * T));
}

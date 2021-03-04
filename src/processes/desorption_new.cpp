#include "desorption_new.h"

Desorption_new::Desorption_new() : m_iNeigh(0) {}

Desorption_new::Desorption_new(int numNeigh, double activationEnergy, double freq, species_new *species, Parameters *parameter)
    : m_iNeigh(numNeigh),
      m_dActNrg(activationEnergy),
      m_freq(freq),
      m_Species(species),
      pParameters(parameter)
{
}

Desorption_new::~Desorption_new() {}

void Desorption_new::perform(int siteID)
{
    m_pLattice->desorp(siteID, m_Species);
}

double Desorption_new::getProbability()
{
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    // Constants
    double Na = 6.0221417930e+23; // Avogadro's number [1/mol]
    double k = 1.3806503e-23;     // Boltzmann's constant [j/K]

    //Case-dependent values read from the input file
    double P = pParameters->getPressure();    // [Pa]
    double T = pParameters->getTemperature(); // [K]

    double E = getActivationEnergy() / Na; //(7.14e+4)/Na;			// [j] -> 17 kcal
    double k_d = getFrequency();             // [s^-1]
    /*--------------------------------------------------*/

    double v0 = k_d;  //*exp(-E/(k*T));

    return v0 * exp(-(double)m_iNeigh * E / (k * T)); //Desorption 1 neigh
}

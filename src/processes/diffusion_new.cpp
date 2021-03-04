#include "diffusion_new.h"

Diffusion_new::Diffusion_new() : m_iNeighNum(0) {}

Diffusion_new::Diffusion_new(int numNeigh, double E_d, double E_m, double freq, species_new *species, Parameters *parameter)
    : m_iNeighNum(numNeigh),
      m_dActNrg(E_d),
      m_Em(E_m),
      m_freq(freq),
      m_Species(species),
      pParameters(parameter) {}

Diffusion_new::~Diffusion_new() {}

void Diffusion_new::perform(int siteID)
{
    m_pLattice->desorp(siteID, m_Species);

    int targetID = m_pLattice->getSite(siteID)->getNeighs().at(rand() % 4)->getID();

    //From this site get the neighobours id and peformt it.
    m_pLattice->adsorp(targetID, m_Species);
}

double Diffusion_new::getProbability()
{

    //These must trenafered in the global definitions
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    // Constants
    double Na = 6.0221417930e+23; // Avogadro's number [1/mol]
    double k = 1.3806503e-23;     // Boltzmann's constant [j/K]

    double P = pParameters->getPressure();    // [Pa]
    double T = pParameters->getTemperature(); // [K]
    double E_d = getActivationEnergy() / Na;              // [j]
    double E_m = getDesorptionEnergy() / Na;              // [j]
    double k_d = getFrequency();                     // [s^-1]
    /*--------------------------------------------------*/

    double v0 = k_d; //*exp(-E/(k*T));
    double A = exp((E_d - E_m) / (k * T));

    //--------------------- Transitions probability ----------------------------------------//
    return A*v0*exp( -m_iNeighNum*E_d/(k*T) );
    //----------------------------------------------------------------------------------------//
}

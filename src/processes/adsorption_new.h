#ifndef ADSORPTION_NEW_H
#define ADSORPTION_NEW_H

#include <iostream>
#include <string>

#include "process_new.h"
#include "site.h"
#include "species_new.h"
#include "parameters.h"

using namespace std;

class Adsorption_new: public Process_new
{
public:
    Adsorption_new();
    Adsorption_new( double actNrg, double molFrac, double sticking, double siteDensity, double mass, species_new* species, Parameters* parameter );
    virtual ~Adsorption_new();

    inline void setActivationEnergy( double nrg ){ m_dActNrg = nrg; }
    inline double getActivationEnergy(){ return m_dActNrg; }

    inline void setMolFrac( double val ){ m_dMolFrac = val; }
    inline double getMolFrac(){ return m_dMolFrac; }

    inline void setTargetSite( Site* site ){ m_Site = site;}
    inline Site* getTargetSite(){ return m_Site; }

    inline void setSpecies( species_new* s ){ m_Species = s; }
    inline species_new* getSpecies(){ return m_Species; }

    inline void setSticking( double sticking ) { m_dStick = sticking; }
    inline double getSticking() { return m_dStick; }

    inline void setSiteDensity( double sitedensity ) { m_siteDensity = sitedensity; }
    inline double getSiteDensity () { return m_siteDensity; }

    inline void setMass( double mass ) { m_mass = mass; }
    inline double getMass() { return m_mass; }

    double getProbability();

    void rules(set<string, std::any>) override {}

    void perform( int siteID  ) override;

private:
    ///The activation energy of the adsorption process
    double m_dActNrg;

    ///The sticking coefficient of the adsorption process
    double m_dStick;

    ///The mole fraction of the adsorption process
    double m_dMolFrac;

    /// Site density [sites/m^2] 
    /// Vlachos code says [moles sites/m^2]
    double m_siteDensity;

    /// mass of system [kg]
    double m_mass;

    ///The site that adsorption will be performed
    Site* m_Site;

    /// The species that must adsopt
    species_new* m_Species;

    /// Pointer to parameters
    Parameters* pParameters;
};

#endif // ADSORPTION_NEW_H

#ifndef DESORPTION_NEW_H
#define DESORPTION_NEW_H

#include "process_new.h"
#include "species_new.h"
#include "parameters.h"

class Desorption_new : public Process_new
{
public:
    Desorption_new();
    Desorption_new(int numNeigh, double activationenergy, double freq, species_new *species, Parameters* parameter );
    ~Desorption_new();

    inline void setActivationEnergy(double nrg) { m_dActNrg = nrg; }
    inline double getActivationEnergy() { return m_dActNrg; }

    inline void setFrequency(double freq) { m_freq = freq; }
    inline double getFrequency() { return m_freq; }

    inline void setTargetSite(Site *site) { m_Site = site; }
    inline Site *getTargetSite() { return m_Site; }

    inline void setSpecies(species_new *s) { m_Species = s; }
    inline species_new *getSpecies() { return m_Species; }

    inline void setNumNeigh(int n) { m_iNeigh = n; }
    inline int getNumNeigh() { return m_iNeigh; }

    double getProbability() override;

    void rules(set<string, std::any>) override {}
    void perform(int siteID) override;

private:
    ///The activation energy of the adsoprtion process
    double m_dActNrg;

    /// Frequency
    double m_freq;

    ///The site that adsorption will be performed
    Site *m_Site;

    ///The species that must be removed from the site
    species_new *m_Species;

    /// Pointer to parameters
    Parameters* pParameters;

    ///The neighbours of this diffusion process
    int m_iNeigh;
};

#endif // DESORPTION_NEW_H

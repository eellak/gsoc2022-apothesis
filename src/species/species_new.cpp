#include "species_new.h"
#include "reaction_new.h"

species_new::species_new() {}

species_new::species_new(string name, double mw, int id) : m_sChemForm(name),
                                                           m_mw(mw),
                                                           m_iID(id){};

species_new::~species_new() {}

void species_new::addReaction(reaction_new *reaction)
{
    m_vReactions.push_back(reaction);
}

vector<reaction_new *> species_new::getReactions()
{
    return m_vReactions;
}

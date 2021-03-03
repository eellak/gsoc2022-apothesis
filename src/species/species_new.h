#ifndef SPECIES_NEW_H
#define SPECIES_NEW_H

#include <iostream>
#include <string>
#include <vector>
#include <string>

using namespace std;

class reaction_new;

class species_new
{
public:
    species_new();
    species_new( string name, double  mw, int id );
    ~species_new();

    inline void setChemFormula( string chemForm ){ m_sChemForm = chemForm; }
    inline string getChemFormula(){ return m_sChemForm; }

    inline void setID( int id ){ m_iID = id; }
    inline int getID(){ return m_iID; }

    inline void setMW( int mw ){ m_mw = mw; }
    inline int getMW(){ return m_mw; }

    void addReaction( reaction_new* reaction );
    vector< reaction_new* > getReactions();

    inline void setMaxReacCoreff( double maxCoeff ){ m_dMaxCoeff = maxCoeff; }
    inline double getMaxReacCoreff(){ return m_dMaxCoeff; }

private:
    /// The chemical formula of this species
    string m_sChemForm;

    /// The id of this species
    int m_iID;

    /// Molecular weight
    double m_mw;

    /// The reactions that species participates in
    vector< reaction_new* > m_vReactions;

    /// This is the maximun number of the coefficient
    /// that this species participates in a reaction
    double m_dMaxCoeff;

};

#endif // SPECIES_NEW_H

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

#ifndef KMC_H
#define KMC_H

#include <map>
#include <list>
#include <vector>
#include <string>
#include <functional>
#include <set>
#include <valarray>
#include "species_new.h"
#include "species/species.h"


#define EXIT { printf("Apothesis terminated. \n"); exit( EXIT_FAILURE ); }

using namespace std;

/** The basic class of the kinetic monte carlo code. */

namespace Utils{ class ErrorHandler; class Parameters;}
namespace SurfaceTiles{ class Site; }
namespace MicroProcesses { class Process; class Adsorption; class Desorption; class Diffusion; class SurfaceReaction;}
namespace RandomGen { class RandomGenerator; }

class Lattice;
class IO;
class Read;

class Apothesis
{
public:
    Apothesis( int argc, char* argv[] );
    virtual ~Apothesis();

    /// Pointers to the classes that will share the common space i.e. the "pointer"
    /// Pointer to the input/output class
    IO* pIO;

    /// Ponter to the read class
    Read* pRead;

    /// Pointer to the lattice class
    Lattice* pLattice;

    /// Pointer to the error class
    Utils::ErrorHandler* pErrorHandler;

    /// Pointer to the paramters class
    Utils::Parameters* pParameters;

    /// Random generator
    RandomGen::RandomGenerator *pRandomGen;

    /// Holds chemical and processes data that are needed
//    Utils::pChemAndProcData;

    /// Intialization of the KMC method. For example here the processes to be performed
    /// as these are written in the input file are constcucted through the factory method
    void init();

    /// Perform the KMC iteratios
    void exec();

    /// Add a process
    void addProcess(string process);

    /// Function to log to output file whether a parameter is properly read
    void logSuccessfulRead(bool read, string parameter);

    /// Return access to list of species
    map<string, Species*> getAllSpecies();

    // Return species
    Species* getSpecies(string species);

    /// Return normalized probabilities of each process
    vector<double> calculateProbabilities(vector<MicroProcesses::Process*>);

    MicroProcesses::Process* getProcessAt(int index, vector<MicroProcesses::Process*> pProcesses);

    MicroProcesses::Process* pickProcess(vector<double> probabilities, double random, vector<MicroProcesses::Process*> pProcesses);

    MicroProcesses::Adsorption* findAdsorption(string species);

    MicroProcesses::Desorption* findDesorption(string species);

    /// Return access to IO pointer
    IO* getIOPointer();

    void setDebugMode(bool);

    bool getDebugMode();

    void setLatticePointer(Lattice* pLattice);

    /// Return access to adsorption
    vector<MicroProcesses::Adsorption*> getAdsorptionPointers();

    /// Return access to adsorption
    vector<MicroProcesses::SurfaceReaction*> getReactionPointers();

private:
    /// The process map which holds all the processes and the sites that each can be performed.
    //Not to handy. Re-think... I have found another way... Implement it
    map< MicroProcesses::Process*, list< SurfaceTiles::Site* >* > m_processMap;

    //This holds the process and the list of sites that can be performed.
    //This should replace m_processMap.
    //Have to check if string is the name of the process or should we put
    //set is log(n) in insert and delete and log(1) for delete
    // "Adsorption:CuAMD" must be a name of a process because we want to be able to do
    // m_procMap["Adsorption:CuAMD"]->addSite(Site);
    map< string, set< int > > m_procMap;

    //The map holding all tbe species. The int is the same as the ID of the species.
    // e.g Assuming the user has procided CuAMD, H2, SiH4, SiH2
    // then m_speciesMAP[ 0 ]= > CuAMD
    //      m_speciesMAP[ 1 ]= > H2
    //      m_speciesMAP[ 2 ]= > SiH4
    //      m_speciesMAP[ 3 ]= > SiH2
    map< int, species_new* > m_speciesMap;

    // Hold the name and the stichiometric coeficient of the reactants in the surface reactions
    // The elements of the vector can be accessed through thr id of the species.

    // e.g. 2CuAMD + 2H2 => ...
    // e.g. CuAMD + S => ...
    // CuAMD: 0
    // H2: 1
    // S: 2
    // e.g. Reaction_0: 2 2 0 0 0 0 0 0 etc.
    //      Reaction_1: 1 0 1 0 0 0 0 0 etc.
    // m_surfReactionsMap[ 0 ] -> will return the stoichiometric coeff for this reactions.
    map< string, valarray<int> > m_surfReacMap;

    // Same as m_surfReacMap holding the procudts. The string should be the same.
    map< string, valarray<int> > m_surfProdMap;

    /// Vector holding the processes to be performed.
    vector< MicroProcesses::Process*> m_vProcesses;

    vector< MicroProcesses::Adsorption*> m_vAdsorption;
    
    vector< MicroProcesses::Desorption*> m_vDesorption;

    vector< MicroProcesses::SurfaceReaction*> m_vSurfaceReaction;

    vector <reference_wrapper<MicroProcesses::SurfaceReaction>> m_refSurfaceReaction;

    /// Vector holding the name of the processes (string)
    vector<string> m_processes;

    /// The number of flags given by the user
    int m_iArgc;

    /// The flags given by the user
    char** m_vcArgv;

    // map of species
    map<string, Species*> m_species;

    // map of interactions
    vector< tuple<string, string> > m_interactions;

    // Set debug mode
    bool m_debugMode;

    /// number of species
    int m_nSpecies;
};

#endif // KMC_H

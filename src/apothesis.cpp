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

#include "pointers.h"
#include "apothesis.h"
#include "lattice/lattice.h"
#include "FCC.h"
#include "io.h"
#include "read.h"
#include "errorhandler.h"
#include "parameters.h"
#include "properties.h"
#include "process.h"
#include "species.h"
#include "string.h"
#include "adsorption.h"
#include "desorption.h"
#include "diffusion.h"
#include "SurfaceReaction.h"
#include "reaction_new.h"
#include "aux/random_generator.h"

#include "processpool.h"

#include "adsorption_new.h"
#include "desorption_new.h"
#include "diffusion_new.h"

#include <numeric>

using namespace MicroProcesses;
using namespace newDesign;

typedef rapidjson::Document Document;
typedef rapidjson::Value Value;
typedef rapidjson::SizeType SizeType;

//using namespace Utils;

Apothesis::Apothesis(int argc, char *argv[])
    : pLattice(0),
      pRead(0),
      m_dProcTime(0.0),
      m_dRTot(0.0),
      m_dProcRate(0.0),
      m_debugMode(false)
{
    m_iArgc = argc;
    m_vcArgv = argv;

    pParameters = new Utils::Parameters(this);
    pProperties = new Utils::Properties(this);
    pRandomGen = new RandomGen::RandomGenerator( this );

    /* This must be constructed before the input */

    // Create input instance
    pIO = new IO(this);
    pRead = new Read(this);

    vector<string> pName = pRead->getSpeciesNames();

    // Build the lattice. This should always follow the read input

    std::cout << "Building the lattice" << std::endl;
    pLattice->setOrientation("110");
    pLattice->build();
    std::cout << "Finished building the lattice" << std::endl;

    // initialize number of species
    m_nSpecies = 0;

    // Initialize Random generator.
    pRandomGen->init( 0 );

    // Initialize process pool
    m_procPool = new newDesign::ProcessPool();
}

Apothesis::~Apothesis()
{
    delete pIO;
    delete pRead;
    delete pLattice;

    // Delete the processes created by the factory method
    for (vector<Process *>::iterator it = m_vProcesses.begin();
         it != m_vProcesses.end(); it++)
    {
        delete *it;
    }
}

void Apothesis::init()
{
    cout << "Opening output file" << endl;
    /// The output file name will come from the user and will have the extenstion .log
    /// This would come as a parameter from the user from the args (also the input).
    /// Now both are hard copied.

    if (!pIO->outputOpen())
        pIO->openOutputFile("Output");

    //---------------------- Creation of the process map & initialization (must be transferred to init) ------------------------------>//
    Adsorption_new* adsorption = new Adsorption_new();
    adsorption->setActivationEnergy( 12.0 );
    adsorption->setName("Adsorption");
    adsorption->setID( 0 );

    pair<string, set<int> > p;
    p.first = adsorption->getName();
    set< int > ids;

    m_procMap.insert( p );
    m_procPool->addProcess( adsorption->getName(), adsorption );
    m_procPool->addProcess( adsorption->getID(),  adsorption);

    for (Site* s:pLattice->getSites() )
        m_procMap[ adsorption->getName() ].insert( s->getID() );

    Desorption_new* desorption_1N = new Desorption_new();
    desorption_1N->setName("Desorption 1N");
    desorption_1N->setNumNeigh( 1 );
    desorption_1N->setID( 1 );

    p.first = desorption_1N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( desorption_1N->getName(), desorption_1N );
    m_procPool->addProcess( desorption_1N->getID(), desorption_1N );

    Desorption_new* desorption_2N = new Desorption_new();
    desorption_2N->setName("Desorption 2N");
    desorption_2N->setNumNeigh( 2 );
    desorption_2N->setID( 2 );

    p.first = desorption_2N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( desorption_2N->getName(), desorption_2N );
    m_procPool->addProcess( desorption_2N->getID(), desorption_2N );

    Desorption_new* desorption_3N = new Desorption_new();
    desorption_3N->setName("Desorption 3N");
    desorption_3N->setNumNeigh( 3 );
    desorption_3N->setID( 3 );

    p.first = desorption_3N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( desorption_3N->getName(), desorption_3N );
    m_procPool->addProcess( desorption_3N->getID(), desorption_3N );

    Desorption_new* desorption_4N = new Desorption_new();
    desorption_4N->setName("Desorption 4N");
    desorption_4N->setNumNeigh( 4 );
    desorption_4N->setID( 4 );

    p.first = desorption_4N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( desorption_4N->getName(), desorption_4N );
    m_procPool->addProcess( desorption_4N->getID(), desorption_4N );

    Desorption_new* desorption_5N = new Desorption_new();
    desorption_5N->setName("Desorption 5N");
    desorption_5N->setNumNeigh( 5 );
    desorption_5N->setID( 5 );

    p.first = desorption_5N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( desorption_5N->getName(), desorption_5N );
    m_procPool->addProcess( desorption_5N->getID(), desorption_5N );

    for (Site* s:pLattice->getSites() )
        m_procMap[ desorption_5N->getName() ].insert( s->getID() );

    Diffusion_new* diffusion_1N = new Diffusion_new();
    diffusion_1N->setName("Diffusion 1N");
    diffusion_1N->setNeigh(1);
    diffusion_1N->setID(6);

    p.first = diffusion_1N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( diffusion_1N->getName(), diffusion_1N );
    m_procPool->addProcess( diffusion_1N->getID(), diffusion_1N );

    Diffusion_new* diffusion_2N = new Diffusion_new();
    diffusion_2N->setName("Diffusion 2N");
    diffusion_2N->setID(7);

    p.first = diffusion_2N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( diffusion_2N->getName(), diffusion_2N );
    m_procPool->addProcess( diffusion_2N->getID(), diffusion_2N );

    Diffusion_new* diffusion_3N = new Diffusion_new();
    diffusion_3N->setName("Diffusion 3N");
    diffusion_3N->setNeigh(3);
    diffusion_3N->setID(8);

    p.first = diffusion_3N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( diffusion_3N->getName(), diffusion_3N );
    m_procPool->addProcess( diffusion_3N->getID(), diffusion_3N );

    Diffusion_new* diffusion_4N = new Diffusion_new();
    diffusion_4N->setName("Diffusion 4N");
    diffusion_4N->setNeigh(4);
    diffusion_4N->setID(9);

    p.first = diffusion_4N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( diffusion_4N->getName(), diffusion_4N );
    m_procPool->addProcess( diffusion_4N->getID(), diffusion_4N );

    Diffusion_new* diffusion_5N = new Diffusion_new();
    diffusion_5N->setName("Diffusion 5N");
    diffusion_5N->setNeigh(5);
    diffusion_5N->setID(10);

    p.first = diffusion_5N->getName();
    m_procMap.insert( p );
    m_procPool->addProcess( diffusion_5N->getName(), diffusion_5N );
    m_procPool->addProcess( diffusion_5N->getID(), diffusion_5N );

    for (Site* s:pLattice->getSites() )
        m_procMap[ diffusion_5N->getName() ].insert( s->getID() );

    //Set reference to each process5
    for (pair<string, set<int> > p:m_procMap)
        m_procPool->getProcessByName( p.first )->setLattice( pLattice );
    //<---------------------- End creation of the process map & initialization  ------------------------------//



    // Processes in this case
    vector<string> pProc = m_processes;

    cout << pProc[0] << endl;

    Document &doc = pRead->getDoc();

    pIO->writeLogOutput("Initializing instances of species");
    if (std::find(pProc.begin(), pProc.end(), "Reaction") == pProc.end())
    {
        vector<double> mws = pRead->getMWs();
        vector<string> names = pRead->getSpeciesNames();

        for (int i = 0; i < mws.size(); ++i)
        {
            Species *s = new Species(names[i], mws[i], m_nSpecies);
            m_nSpecies++;
            m_species[names[i]] = s;
        }
    }

    // Initializing interactions between species
    pIO->writeLogOutput("Reading interactions between species");
    Value &speciesName = doc["Species"];
    for (Value::ConstMemberIterator itr = speciesName.MemberBegin(); itr != speciesName.MemberEnd(); ++itr)
    {
        const char *name = itr->name.GetString();
        Value &singleSpecies = speciesName[name];
        if (singleSpecies.HasMember("Interactions"))
        {
            Value &interactions = singleSpecies["Interactions"];

            pIO->writeLogOutput("Initializing interactions between species...");

            int numInters = interactions.Size();
            for (int i = 0; i < numInters; i++)
            {
                m_interactions.push_back(make_tuple(name, interactions[i].GetString()));
                pIO->writeLogOutput("(" + get<0>(m_interactions.back()) + ", " + get<1>(m_interactions.back()) + ")");

            }
        }
    }

    pIO->writeLogOutput("Initializing processes");

    if (std::find(pProc.begin(), pProc.end(), "Adsorption") != pProc.end())
    {
        pIO->writeLogOutput("Initializing Adsorption");

        // Read parameters for Adsorption
        Value &specie = doc["Process"]["Adsorption"]["Species"];
        Value &stick = doc["Process"]["Adsorption"]["Sticking"];
        Value &mFrac = doc["Process"]["Adsorption"]["MassFraction"];
        Value &ads = doc["Process"]["Adsorption"];

        // Verify presence of each parameter in input file
        logSuccessfulRead(specie.IsArray(), "Adsorption species");
        logSuccessfulRead(stick.IsArray(), "Adsorption sticking coefficients");
        logSuccessfulRead(mFrac.IsArray(), "Adsorption mass fraction");

        // Initialize vectors
        vector<string> species;
        vector<double> sticking;
        vector<double> massFraction;

        // Sum of mass fraction. Later used to normalize.
        double sum = 0;

        for (SizeType i = 0; i < specie.Size(); i++)
        {
            // Output possible errors
            if (!specie[i].IsString())
                pErrorHandler->error_simple_msg("Species format is not a string");
            if (!stick[i].IsNumber())
                pErrorHandler->error_simple_msg("Sticking coefficient format is not a double");
            if (!mFrac[i].IsNumber())
                pErrorHandler->error_simple_msg("Mass fraction format is not a double");

            // Push values to corresponding vectors
            species.push_back(specie[i].GetString());
            sticking.push_back(stick[i].GetDouble());
            massFraction.push_back(mFrac[i].GetDouble());
            sum += massFraction[i];
        }

        // Normalize the values of the mass fraction
        for (vector<double>::iterator itr = massFraction.begin(); itr != massFraction.end(); ++itr)
        {
            *itr = *itr / sum;
        }

        for (int i = 0; i < species.size(); ++i)
        {
            bool direct = false;
            if (ads.HasMember("Direct"))
            {
                string dir = ads["Direct"].GetString();
                if (dir.compare("Yes") == 0)
                {
                    direct = true;
                }
                else
                {
                    // Set as warning, in case an unexpected result is received. Default value is false.
                    pErrorHandler->warningSimple_msg("The adsorption of species " + species[i] + " has direct value " + dir);
                }

            }
            Adsorption* a = new Adsorption(this, species[i], m_species[species[i]], sticking[i], massFraction[i], direct);

            // Keep two separate vectors: one for all processes, one for adsorption processes only
            m_vAdsorption.push_back(a);
            m_vProcesses.push_back(a);
        }

        // Resolving Conflicts

        pIO->writeLogOutput("...Done initializing Adsorption process.");
    }
    if (std::find(pProc.begin(), pProc.end(), "Desorption") != pProc.end())
    {
        pIO->writeLogOutput("Initializing Desorption");

        // Read parameters for Desorption
        Value &vSpecie = doc["Process"]["Desorption"]["Species"];
        Value &vEnergy = doc["Process"]["Desorption"]["Energy"];
        Value &vFreq = doc["Process"]["Desorption"]["Frequency"];

        // Verify presence of each parameter in input file
        logSuccessfulRead(vSpecie.IsArray(), "Desorption species");
        logSuccessfulRead(vEnergy.IsArray(), "Desorption energy");
        logSuccessfulRead(vFreq.IsArray(), "Desorption frequency");

        // Initialize vectors
        vector<string> species;
        vector<double> energy;
        vector<double> frequency;

        for (SizeType i = 0; i < vSpecie.Size(); i++)
        {
            // Output possible errors
            if (!vSpecie[i].IsString())
                pErrorHandler->error_simple_msg("Species format is not a string");
            if (!vEnergy[i].IsNumber())
                pErrorHandler->error_simple_msg("Desorption energy format is not a number");
            if (!vFreq[i].IsNumber())
                pErrorHandler->error_simple_msg("Desorption frequency format is not a number");

            // Push values to corresponding vectors
            species.push_back(vSpecie[i].GetString());
            energy.push_back(vEnergy[i].GetDouble());
            frequency.push_back(vFreq[i].GetDouble());
        }

        // Add process to m_vProcesses
        for (int i = 0; i < species.size(); ++i)
        {
            // Call function to find adsorption class
            Adsorption *pAdsorption = findAdsorption(species[i]);

            // Create new instance of desorption class
            Desorption *d = new Desorption(this, species[i], m_species[species[i]], energy[i], frequency[i]);

            // Check if associated adsorption class is a valid (non-null) pointer. Associate desorption class with appropriate adsorption and vice versa
            if (pAdsorption)
            {
                d->setAdsorptionPointer(pAdsorption);
                pAdsorption->setDesorptionPointer(d);
                pAdsorption->setDesorption(true);
            }

            // if not, create associations in adsorption and desorption

            // else, output warning

            m_vDesorption.push_back(d);
            m_vProcesses.push_back(d);
        }
        pIO->writeLogOutput("...Done initializing desorption process.");
    }
    if (std::find(pProc.begin(), pProc.end(), "Diffusion") != pProc.end())
    {
        pIO->writeLogOutput("Initializing Diffusion");

        // Read parameters for Diffusion
        Value &vSpecie = doc["Process"]["Diffusion"]["Species"];
        Value &vEnergy = doc["Process"]["Diffusion"]["Energy"];
        Value &vFreq = doc["Process"]["Diffusion"]["Frequency"];

        // Verify presence of each parameter in input file
        logSuccessfulRead(vSpecie.IsArray(), "Diffusion species");
        logSuccessfulRead(vEnergy.IsArray(), "Diffusion energy");
        logSuccessfulRead(vFreq.IsArray(), "Diffusion frequency");

        // Initialize vectors
        vector<string> species;
        vector<double> energy;
        vector<double> frequency;

        for (SizeType i = 0; i < vSpecie.Size(); i++)
        {
            // Output possible errors
            if (!vSpecie[i].IsString())
                pErrorHandler->error_simple_msg("Species format is not a string");
            if (!vEnergy[i].IsNumber())
                pErrorHandler->error_simple_msg("Diffusion energy format is not a number");
            if (!vFreq[i].IsNumber())
                pErrorHandler->error_simple_msg("Diffusion frequency format is not a number");

            // Push values to corresponding vectors
            species.push_back(vSpecie[i].GetString());
            energy.push_back(vEnergy[i].GetDouble());
            frequency.push_back(vFreq[i].GetDouble());
        }

        // Add process to m_vProcesses
        for (int i = 0; i < species.size(); ++i)
        {
            // Call function to find adsorption class
            Adsorption *pAdsorption = findAdsorption(species[i]);
            Desorption *pDesorption = findDesorption(species[i]);

            Diffusion *diff = new Diffusion(this, species[i], energy[i], frequency[i]);
            diff->setAdsorptionPointer(pAdsorption);
            diff->setDesorptionPointer(pDesorption);

            pAdsorption->setDiffusion(true);
            pAdsorption->setDiffusionPointer(diff);

            pDesorption->setDiffusion(true);
            pDesorption->setDiffusionPointer(diff);

            m_vProcesses.push_back(diff);
        }
        pIO->writeLogOutput("...Done initializing diffusion process.");
    }
    if (std::find(pProc.begin(), pProc.end(), "Reaction") != pProc.end())
    {
        pIO->writeLogOutput("Initializing Reaction");

        // Read parameters for Reaction
        Value &pRxn = doc["Process"]["Reaction"];
        Value &vSpecie = doc["Process"]["Reaction"]["Species"];
        Value &vStoich = doc["Process"]["Reaction"]["Stoichiometry"];
        Value &vEnergy = doc["Process"]["Reaction"]["Energy"];
        Value &vPreExp = doc["Process"]["Reaction"]["PreExp"];

        // Verify presence of each parameter in input file
        logSuccessfulRead(vSpecie.IsArray(), "Reaction species");
        logSuccessfulRead(vStoich.IsArray(), "Reaction Stoichiometry");
        logSuccessfulRead(vEnergy.IsNumber(), "Enthalpy of reaction");
        logSuccessfulRead(vPreExp.IsNumber(), "Pre-exponential factor for the reaction");

        // Initialize vectors
        vector<Species *> species;
        vector<double> stoichiometry;
        double energy;
        double preexp;

        auto mws = pRead->getMWs();

        // Loop through species
        for (SizeType i = 0; i < vSpecie.Size(); i++)
        {
            // Output possible errors
            if (!vSpecie[i].IsString())
                pErrorHandler->error_simple_msg("Species format is not a string");
            if (!vStoich[i].IsNumber())
                pErrorHandler->error_simple_msg("Diffusion energy format is not a number");

            // Push values to corresponding vectors
            Species *s = getSpecies(vSpecie[i].GetString());
            stoichiometry.push_back(vStoich[i].GetDouble());

            string name = vSpecie[i].GetString();
            // If the species is not already defined, push new member onto maps
            if (m_species[name] == NULL)
            {
                Species *s = new Species(name, mws[i], vStoich[i].GetDouble(), m_nSpecies);
                m_nSpecies++;
                m_species[name] = s;
                species.push_back(s);
            }
            else
            {
                species.push_back(m_species[name]);
            }
        }

        // Store value for energy and pre-exponential factor
        energy = vEnergy.GetDouble();
        preexp = vPreExp.GetDouble();

        // Check mass balance on reaction
        double cumulativemass = 0;
        for (map<string, Species *>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
        {
            Species *s = itr->second;
            cumulativemass += s->getMW() * s->getStoicCoeff();
        }

        if (abs(cumulativemass) > 1e-10)
        {
            cout << "Warning! Mass balance of Reaction is not balanced" << endl;
        }

        // Read immobilization variable
        bool immobilized = true;
        if (pRxn.HasMember("Immobilize"))
        {
            immobilized = pRxn["Immobilize"].GetBool();
        }

        SurfaceReaction *s = new SurfaceReaction(this, species, stoichiometry, energy, preexp, immobilized);
        m_vProcesses.push_back(s);
        m_vSurfaceReaction.push_back(s);
        pIO->writeLogOutput("...Done initializing reaction.");
    }


    // Initialize interactions between adsorption species and classes
    vector<tuple<string, string>>::iterator itr = m_interactions.begin();
    // For each pair of interactions, add the possible interaction (2-way) within the appropriate pointers
    for (itr; itr != m_interactions.end(); ++itr)
    {
        tuple<string, string> temp = *itr;
        string spec1 = get<0>(temp);
        string spec2 = get<1>(temp);

        Species *s1 = m_species[spec1];
        Species *s2 = m_species[spec2];

        Adsorption *pAdsorption1 = findAdsorption(spec1);
        pAdsorption1->addInteraction(s2);

        Adsorption *pAdsorption2 = findAdsorption(spec2);
        pAdsorption2->addInteraction(s1);
    }
    vector<Adsorption*> adsorptionpointers = getAdsorptionPointers();
    for (vector<Adsorption *>::iterator iter = adsorptionpointers.begin(); iter != adsorptionpointers.end(); ++iter)
    {
        Adsorption *pAds = *iter;
        vector<Species *> possibleInteractions = pAds->getInteractions();
        bool found = false;
        for (int i = 0; i < possibleInteractions.size(); ++i)
        {
            cout<< possibleInteractions[i]->getName() << ", ";
        }
    }
    cout<<endl;
    cout<<endl;
    /// First the processes that participate in the simulation

    /// that were read from the file input and the I/O functionality
    //m_vProcesses[0]->setInstance( this );
    for (vector<Process *>::iterator itr = m_vProcesses.begin(); itr != m_vProcesses.end(); ++itr)
    {
        Process *p = *itr;
        p->activeSites(pLattice);
    }

    // Initialize species map in lattice
    vector<Site*> sites = pLattice->getSites();
    for (int site = 0; site < sites.size(); ++site)
    {
        Site* pSite = sites[site];
        pSite->initSpeciesMap(m_nSpecies);
    }

}

void Apothesis::exec()
{
    
    //--------------- Open files for writting ---------------------->
    pIO->openRoughnessFile( "testRough" );

    // Here we set the process map to the lattice in order for the lattice to be able to modified according to the structural properties of the lattice.
    //(must be transferred to init)
    pLattice->setProcMap( &m_procMap );

    //Perform the number of KMC steps read from the input.
    m_dEndTime = pParameters->getEndTime();

    //if (m_dEndTime == 0.0){
    //    pErrorHandler->error_simple_msg("Zero iterations found.");
    //    exit(0);
    //}
    //else
    //    pIO->writeLogOutput("Running " + to_string( m_dEndTime ) + " sec");

    //Calculate the total probability (R) --------------------------//
    m_dRTot = 0.0;
    for (pair<string, set<int> > p:m_procMap)
        m_dRTot += (double) m_procPool->getProcessByName( p.first )->getProbability()*p.second.size();

    m_dEndTime = 0.01;

    pProperties->calculateRoughness();
    pIO->writeRoughness( m_dProcTime, pProperties->getRoughness() );

    while ( m_dProcTime <= m_dEndTime ){
        //1. Get a random numbers
        m_dRandom = pRandomGen->getDoubleRandom();

        m_dSum = 0.0;
        for (pair<string, set<int> > p:m_procMap) {
            m_dProcRate = (double) m_procPool->getProcessByName( p.first )->getProbability()*p.second.size();
            m_dSum += m_dProcRate/m_dRTot;

            //2. Pick a process according to the rates
            if ( m_dRandom <= m_dSum ){

                //Get a random number which is the ID of the site where this process can performed
                m_iSiteNum = pRandomGen->getIntRandom(0, m_procMap[ p.first ].size() - 1 );

                //3. From this process pick the random site with id and perform it:
                m_procPool->getProcessByName( p.first )->perform( *next(m_procMap[ p.first ].begin(), m_iSiteNum) );
                m_procPool->getProcessByName( p.first )->eventHappened();

                //4. Re-compute the processes rates and re-compute Rtot (see ppt).
                m_dRTot = 0.0;
                for (pair<string, set<int> > p:m_procMap)
                    m_dRTot += (double) m_procPool->getProcessByName( p.first )->getProbability()*p.second.size();

                //5. Compute dt = -ln(ksi)/Rtot
                m_dt = -log( pRandomGen->getDoubleRandom()  )/m_dRTot;
                break;
            }
        }
        //6. advance time: time += dt;
        m_dProcTime += m_dt;

        //Calulate the roughness
        pProperties->calculateRoughness();

        //Then write it
        pIO->writeRoughness( m_dProcTime, pProperties->getRoughness() );
    }

    delete m_procPool;
}

void Apothesis::addProcess(string process)
{
    m_processes.push_back(process);
}

void Apothesis::logSuccessfulRead(bool read, string parameter)
{
    if (!pIO->outputOpen())
    {
        pIO->openOutputFile("Output");
    }

    read ? pIO->writeLogOutput("Reading " + parameter)
         : pErrorHandler->error_simple_msg("No " + parameter + " found in input file");
}

map<string, Species *> Apothesis::getAllSpecies()
{
    return m_species;
}

Species *Apothesis::getSpecies(string species)
{
    return m_species[species];
}

vector<double> Apothesis::calculateProbabilities(vector<Process *> pProcesses)
{
    /// Calculate probabilities for each process
    int numProcesses = pProcesses.size();

    /// Calculate probability of each
    vector<double> probability;
    double total = 0;
    vector<Process *>::iterator itr = pProcesses.begin();

    // Find the probability for each process. Push onto prob.
    for (; itr != pProcesses.end(); ++itr)
    {
        Process *process = *itr;
        double prob = process->getProbability();
        probability.push_back(prob + total);
        total += prob;
    }

    // Normalize all values
    for (vector<double>::iterator itr = probability.begin(); itr != probability.end(); ++itr)
    {
        *itr = *itr / total;
    }

    return probability;
}

// May be possible to delete (if pProcesses is a vector, can simply access by index)
Process *Apothesis::getProcessAt(int index, vector<Process *> pProcesses)
{
    vector<Process *>::iterator it = pProcesses.begin();
    std::advance(it, index);
    return *it;
}

Process *Apothesis::pickProcess(vector<double> probabilities, double random, vector<Process *> pProcesses)
{
    for (int index = 0; index < probabilities.size(); ++index)
    {
        if (random < probabilities[index])
        {
            return pProcesses[index];
        }
    }
    return pProcesses[probabilities.size() - 1];
}

Adsorption *Apothesis::findAdsorption(string species)
{
    for (vector<Adsorption *>::iterator itr = m_vAdsorption.begin(); itr != m_vAdsorption.end(); ++itr)
    {
        Adsorption *a = *itr;
        if (!species.compare(a->getSpeciesName()))
        {
            return a;
        }
        // find name of adsorption species
    }
    pErrorHandler->warningSimple_msg("Warning! Could not find instance of Adsorption class for desorbed species " + species);
}

Desorption *Apothesis::findDesorption(string species)
{
    for (vector<Desorption *>::iterator itr = m_vDesorption.begin(); itr != m_vDesorption.end(); ++itr)
    {
        Desorption *d = *itr;
        if (!species.compare(d->getSpeciesName()))
        {
            return d;
        }
        // find name of adsorption species
    }
    cout << "Warning! Could not find instance of Adsorption class for desorbed species " << species << endl;
}

IO *Apothesis::getIOPointer()
{
    return pIO;
}

void Apothesis::setDebugMode(bool ifDebug)
{
    m_debugMode = ifDebug;
}

bool Apothesis::getDebugMode()
{
    return m_debugMode;
}

void Apothesis::setLatticePointer(Lattice *lattice)
{
    pLattice = lattice;
}

vector<Adsorption *> Apothesis::getAdsorptionPointers()
{
    return m_vAdsorption;
}

vector<SurfaceReaction *> Apothesis::getReactionPointers()
{
    return m_vSurfaceReaction;
}

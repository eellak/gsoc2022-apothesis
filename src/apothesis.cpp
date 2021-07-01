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
#include "aux/random_generator.h"
#include <bits/stdc++.h>
#include "factory_process.h"
#include <numeric>
#include "adsorption_simple_cubic.h"

using namespace MicroProcesses;

typedef rapidjson::Document Document;
typedef rapidjson::Value Value;
typedef rapidjson::SizeType SizeType;

//using namespace Utils;

Apothesis::Apothesis(int argc, char *argv[])
    : pLattice(0),
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

    //For building with steps surface. Works only for simple cubic
    pLattice->buildSteps( 20, 1, 0);

    std::cout << "Finished building the lattice" << std::endl;

    // initialize number of species
    m_nSpecies = 0;
}

Apothesis::~Apothesis()
{
    delete pIO;
    delete pRead;
    delete pLattice;
}

void Apothesis::init()
{
                

    cout << "Opening output file" << endl;
    /// The output file name will come from the user and will have the extenstion .log
    /// This would come as a parameter from the user from the args (also the input).
    /// Now both are hard copied.


    if (!pIO->outputOpen())
        pIO->openOutputFile("Output");
    // Initialize Random generator.
    pRandomGen->init( 213212 );

    // Get access to document file
    Document &doc = pRead->getDoc();

    // --------------------------- set case number ---------------------------- //
    // Note to Nikos: 
    // Case 0 = Simple cubic, only CuAMD is adsorbed. Desorption of CuAMD. No reaction or diffusion. Low temperature regime.
    // Case 1 = Simple cubic, CuAMD is adsorbed, reacts to HAMD, which is then desorbed.
    // Case 2 = Simple cubic, CuAMD and H2 are both adsorbed, reacts, wait for HAMD desorption.
    
    setCaseStudy(0);

    // Set end time
    m_dEndTime = pParameters->getEndTime();

    map<string, any> params;
    params.insert( {"T", pParameters->getTemperature()} );
    params.insert( {"P", pParameters->getPressure()} );
    
      
    set< Site* > emptySet;

    pIO->writeLogOutput("Initializing instances of species");


    vector<double> mws = pRead->getMWs();
    vector<string> names = pRead->getSpeciesNames();
//
    for (int i = 0; i < mws.size(); ++i)
    {
      species_new *s = new species_new(); //Species(names[i], mws[i], m_nSpecies);
      s->setChemFormula(names[i]);
      s->setID(i);
      m_nSpecies++;
      m_speciesMap[i] = s;
    }
  

    // Read parameters for Reaction
    Value& pProc = doc["Process"];

    for (Value::ConstMemberIterator itr = pProc.MemberBegin(); itr != pProc.MemberEnd(); ++itr)
    {
        string rxnName = itr->name.GetString();
        if (rxnName.compare("Reaction") == 0)
        {
            Value& pRxn = doc["Process"]["Reaction"];
            int specCounter = 0;
            for (Value::ConstMemberIterator itr = pRxn.MemberBegin(); itr != pRxn.MemberEnd(); ++itr)
            {
                const char *rxnName = itr->name.GetString();
                Value &vSpecies = pRxn[rxnName]["Species"];
                Value &vStoich = pRxn[rxnName]["Stoichiometry"];
                Value &vEnergy = pRxn[rxnName]["Energy"];
                Value &vPreExp = pRxn[rxnName]["PreExp"];

                auto pos = m_processMap.insert( { FactoryProcess::createProcess("Reaction"), emptySet } );
                pos.first->first->setName(rxnName);
                pos.first->first->setActivationEnergy(vEnergy.GetDouble());
                pos.first->first->setPreExpFactor(vPreExp.GetDouble());
                vector<pair<int, species_new*>> reactants;
                vector<pair<int, species_new*>> products;

                for (SizeType counter = 0; counter < vSpecies.Size(); ++counter)
                {
                    species_new* spec = new species_new();
                    spec->setChemFormula(vSpecies[counter].GetString());
                    spec->setID(specCounter);
                    spec->setMaxReacCoreff(vStoich[counter].GetDouble());
                    m_speciesMap[specCounter] = spec;
                    ++specCounter;

                    if (vStoich[counter].GetDouble() < 0)
                    {
                        reactants.push_back(make_pair(-1*vStoich[counter].GetDouble(), spec));
                    }
                    else
                    {
                        products.push_back(make_pair(1*vStoich[counter].GetDouble(), spec));
                    }

                }
                params.insert({"reactants", reactants});
                params.insert({"products", products});
                pos.first->first->setParams( params );
            }

        }
        else if (rxnName.compare("Adsorption") == 0)
        {
            // Read parameters for Adsorption
            Value &specie = doc["Process"]["Adsorption"]["Species"];
            Value &stick = doc["Process"]["Adsorption"]["Sticking"];
            Value &mFrac = doc["Process"]["Adsorption"]["MolFraction"];
            Value &ctot = doc["Process"]["Adsorption"]["C_tot"];
            Value &mass = doc["Process"]["Adsorption"]["Mass"];
            Value &ads = doc["Process"]["Adsorption"];

            // Verify presence of each parameter in input file
            logSuccessfulRead(specie.IsArray(), "Adsorption species");
            logSuccessfulRead(stick.IsArray(), "Adsorption sticking coefficients");
            logSuccessfulRead(mFrac.IsArray(), "Adsorption mass fraction");
            logSuccessfulRead(ctot.IsDouble(), "Site density");
            logSuccessfulRead(mass.IsDouble(), "Mass");

            params.insert( {"R", 8.3145 });
            params.insert( {"k", 1.3806503e-23} );
            params.insert( {"Na", 6.0221417930e+23} );    

            for (SizeType spec = 0; spec < specie.Size(); ++spec)
            {
                string f = "f_";
                f.append(specie[spec].GetString());
                params.insert( {f, mFrac[spec].GetDouble()} );
                params.insert( {"C_tot", ctot.GetDouble()} );
                params.insert( {"s0", stick[spec].GetDouble()} );
                params.insert( {"mass", mass.GetDouble()} );

                auto pos = m_processMap.insert( { FactoryProcess::createProcess("AdsorptionSimpleCubic"), emptySet } );
                pos.first->first->setName("Adsorption");
                pos.first->first->setUncoAccepted( true );
                params.insert({"ActivationEnergy", 12.0});
                pos.first->first->setParams( params );
                pos.first->first->setLattice( pLattice );
                pos.first->first->setRandomGen( pRandomGen );
                pos.first->first->setSpecies( m_speciesMap[spec] );
                pos.first->first->setApothesis(this);
                pos.first->first->getParameter("ActivationEnergy");
                
                for ( Site* s:pLattice->getSites() ){
                        if ( pos.first->first->rules( s ) )
                            pos.first->second.insert( s );
                    }
            }
        }
        else if (rxnName.compare("Desorption") == 0)
        {
        
            // Read parameters for Desorption
            Value &vSpecie = doc["Process"]["Desorption"]["Species"];
            Value &vEd = doc["Process"]["Desorption"]["Ed"];
            Value &vEm = doc["Process"]["Desorption"]["Em"];
            Value &vFreq = doc["Process"]["Desorption"]["Frequency"];
            Value &vNeigh = doc["Process"]["Desorption"]["n_neigh"];

            // Verify presence of each parameter in input file
            logSuccessfulRead(vSpecie.IsArray(), "Desorption species");
            logSuccessfulRead(vEd.IsArray(), "Desorption energy");
            logSuccessfulRead(vFreq.IsArray(), "Desorption frequency");
            logSuccessfulRead(vNeigh.IsArray(), "Number of neighbours");
            
            params.insert( {"Na", 6.0221417930e+23} );
            for (SizeType spec = 0; spec < vSpecie.Size(); ++spec)
            {
                params.insert( {"Species" + spec, vSpecie[spec].GetString()} );
                params.insert( {"E_d", vEd[spec].GetDouble()/6.0221417930e+23 } );
                params.insert( {"E_m", vEm[spec].GetDouble()/6.0221417930e+23 } );
                params.insert( {"Freq", vFreq[spec].GetDouble()/6.0221417930e+23 } );
                params.insert ( {"neighs", vNeigh[spec].GetInt()});

                for (int i = 0; i < vNeigh[spec].GetInt(); ++i)
                {
                    auto pos = m_processMap.insert( { FactoryProcess::createProcess("DesorptionSimpleCubic"), emptySet } );
                    string name = "Desorption" ;
                    name.append(to_string(i)).append("N");
                    pos.first->first->setName(name);
                    pos.first->first->setParams( params );
                    pos.first->first->setLattice( pLattice );
                    pos.first->first->setRandomGen( pRandomGen );
                    pos.first->first->setSpecies( m_speciesMap[spec] );

                    for ( Site* s:pLattice->getSites() ){
                        if ( pos.first->first->rules( s ) )
                            pos.first->second.insert( s );
                    }
                }
            }
        }
    }
    

    
        
        
   
    
    
   
    

//    pLattice->print();
    cout << " Lets see! " << endl;
}

void Apothesis::exec()
{
    Site* tempSite = 0;
    //--------------- Open files for writting ---------------------->
    pIO->openRoughnessFile( "testRough" );
    pIO->writeLogOutput("Running " + to_string( m_dEndTime ) + " sec");

    //Calculate first time the total probability (R) --------------------------//
    m_dRTot = 0.0;
    for (pair<Process*, set< Site* > > p:m_processMap)
        m_dRTot += p.first->getProbability()*(double)p.second.size();

    pIO->writeInOutput( "\n" );
    pIO->writeInOutput( "********************************************************************" );
    string output = "Time (s)"s + '\t' + "Micro roughness (-)" + '\t' + "RMS (-)" + '\t' ;

    for ( auto &p:m_processMap)
        output += p.first->getName() + '\t';

    for ( auto &p:m_processMap)
        output += "Coverage " + p.first->getName() + '\t';
    pIO->writeInOutput( output );

    int iTimeStep = 0;
    pIO->writeLatticeHeights( m_dProcTime, iTimeStep );

    double writeLatHeigsEvery = 1e-3; //in s
    double timeToWrite = 0.0;

    output = std::to_string(m_dProcTime) + '\t' + std::to_string( pProperties->getMicroroughness() ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
    for ( auto &p:m_processMap)
        output += std::to_string( p.first->getNumEventHappened() ) + '\t';

    for ( auto &p:m_processMap)
        output += std::to_string( (double)p.second.size()/(double)pLattice->getSize() ) + '\t';
    pIO->writeInOutput( output );


    pLattice->writeXYZ( "lattice.xyz" );

    while ( m_dProcTime <= m_dEndTime ){
        //1. Get a random numbers
        m_dSum = 0.0;
        m_dRandom = pRandomGen->getDoubleRandom();

        for ( auto &p:m_processMap){
            m_dProcRate = p.first->getProbability()*p.second.size(); //TODO: verify getProbabitliy() != 0
            m_dSum += m_dProcRate/m_dRTot;
    
            //2. Pick a process according to the rates
            if ( m_dRandom <= m_dSum ){
                //Get a random number which is the ID of the site where this process can performed
                m_iSiteNum = pRandomGen->getIntRandom(0, p.second.size() - 1 );

                //3. From this process pick the random site with id and perform it:
                Site* s = *next( p.second.begin(), m_iSiteNum );

                p.first->perform( s );
                tempSite = s;
                //Count the event for this class
                p.first->eventHappened();

                // Check if an affected site must enter to a class or not
                for (Site* affectedSite:p.first->getAffectedSites() ){
                    //Erase the affected site from the processes
                    for (auto &p2:m_processMap){
                        if ( !p2.first->isUncoAccepted() ) {
                            //Added if it obeys the rules of this process
                            if ( p2.first->rules( affectedSite ) ) {
                                if (p2.second.find( affectedSite ) == p2.second.end() )
                                    p2.second.insert( affectedSite );
                            }
                            else
                                p2.second.erase( affectedSite );

                            if ( p2.second.empty() && p2.first->getName() == "AdsorptionFCC1102SSimple") {
                                cout << "No place to adsorb! " << endl;
                                exit(0);
                            }

                        }
                    }
                }


                //4. Re-compute the processes rates and re-compute Rtot (see ppt).
                m_dRTot = 0.0;
                for (pair<Process*, set< Site* > > p3:m_processMap)
                    m_dRTot += p3.first->getProbability()*p3.second.size();

                //5. Compute dt = -ln(ksi)/Rtot
                m_dt = -log( pRandomGen->getDoubleRandom()  )/m_dRTot;
                break;
            }
        }

        //6. advance time: time += dt;
        m_dProcTime += m_dt;
        timeToWrite += m_dt;
        cout.precision(17);
        cout << m_dProcTime << endl;

        pLattice->writeXYZ( "test.xzy" );

        //Write the lattice heights
        iTimeStep++;

//        cout << "******************" << endl;
//        cout <<" Site " << tempSite->getID() << endl;
//        cout << "******************" << endl;
//        pLattice->print();
//        cout << "******************" << endl;

        if ( timeToWrite >= writeLatHeigsEvery ) {
            pIO->writeLatticeHeights( m_dProcTime, iTimeStep );

            output = std::to_string(m_dProcTime) + '\t' + std::to_string( pProperties->getMicroroughness() ) + '\t' + std::to_string( pProperties->getRMS() )  + '\t' ;
            for ( auto &p:m_processMap)
                output += std::to_string( p.first->getNumEventHappened() ) + '\t';

            for ( auto &p:m_processMap)
                output += std::to_string( (double)p.second.size() ) + '\t';
            pIO->writeInOutput( output );

            timeToWrite = 0.0;
        }
    }
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

map<int, species_new*> Apothesis::getAllSpecies()
{
    return m_speciesMap;
}

species_new *Apothesis::getSpecies(int speciesID)
{
    return m_speciesMap[speciesID];
}



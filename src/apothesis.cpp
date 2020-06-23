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
#include "lattice.h"
#include "io.h"
#include "read.h"
#include "errorhandler.h"
#include "parameters.h"
#include "process.h"
#include "species.h"
#include "string.h"
#include "adsorption.h"
#include<numeric>

using namespace MicroProcesses;

typedef rapidjson::Document Document;
typedef rapidjson::Value Value;
typedef rapidjson::SizeType SizeType;

//using namespace Utils;

Apothesis::Apothesis( int argc, char* argv[] ):pLattice( 0 ),pRead( 0 )
  {
    m_iArgc = argc;
    m_vcArgv = argv;

    pParameters = new Utils::Parameters(this);

    /* This must be constructed before the input */
    pLattice = new Lattice( this );

    // Create input instance
    pIO = new IO(this);
    //pIO->init( m_iArgc, m_vcArgv );

    pRead = new Read(this);

    vector<string> pName = pRead->getSpeciesNames();

    for (vector<string>::iterator it = pName.begin(); it != pName.end(); ++it)
    {
      int counter = 0;
      // Read the molecular weights from the pName file
      vector<double> pMW = pRead->getMWs();
      m_species.push_back( new Species(*it, pMW[counter], 0));
    }


    

    // Read the input
    //pIO->readInputFile();

    // Build the lattice. This should always follow the read input

    std::cout<<"Building the lattice"<<std::endl;
    pLattice->build();
    std::cout<<"Finished building the lattice"<<std::endl;
  }

Apothesis::~Apothesis()
  {
  delete pIO;
  delete pRead;
  delete pLattice;

  // Delete the processes created by the factory method
  for (vector<Process*>::iterator it = m_vProcesses.begin();
         it != m_vProcesses.end(); it++) {
       delete *it;
    }
  }

void Apothesis::init()
{
   cout<<"Opening output file"<<endl;
  /// The output file name will come from the user and will have the extenstion .log
  /// This would come as a parameter from the user from the args (also the input).
  /// Now both are hard copied.
  pIO->openOutputFile( "Output");

  // Processes in this case
  vector<string> pProc = pRead->getSpeciesNames();

  Document& doc = pRead->getDoc();

  pIO->writeLogOutput("Initializing processes");

  if (std::find(pProc.begin(), pProc.end(), "Adsorption") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Adsorption");
    
    // Read parameters for Adsorption
    Value& spec = doc["Process"]["Adsorption"]["Species"];
    Value& stick = doc["Process"]["Adsorption"]["Sticking"];
    Value& mFrac = doc["Process"]["Adsorption"]["MassFraction"];
  
    // Verify presence of each parameter in input file
    logSuccessfulRead(spec.IsArray(), "Adsorption species");
    logSuccessfulRead(stick.IsArray(), "Adsorption sticking coefficients");
    logSuccessfulRead(mFrac.IsArray(), "Adsorption mass fraction");
    

    // Initialize vectors
    vector<string> species;
    vector<double> sticking;
    vector<double> massFraction;

    for(SizeType i = 0; i < spec.Size(); i++)
    {
      // Output possible errors
      if (!spec[i].IsString())
        pErrorHandler->error_simple_msg("Species format is not a string");
      if(!stick[i].IsNumber())
        pErrorHandler->error_simple_msg("Sticking coefficient format is not a double");
      if(!mFrac[i].IsNumber())
        pErrorHandler->error_simple_msg("Mass fraction format is not a double");

      // Push values to corresponding vectors
      species.push_back(spec[i].GetString());
      sticking.push_back(stick[i].GetDouble());
      massFraction.push_back(mFrac[i].GetDouble());
    }

    // Normalize the values of the mass fraction
    double sum = std::accumulate(massFraction.begin(), massFraction.end(), 0);
    for (vector<double> :: iterator itr = massFraction.begin(); itr != massFraction.end(); ++itr)
    {
      *itr = *itr/sum;
    }

    // Add process to m_vProcesses
    m_vProcesses.push_back(new Adsorption (species, sticking));

  }
  if (std::find(pProc.begin(), pProc.end(), "Desorption") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Desorption");
  }
  if (std::find(pProc.begin(), pProc.end(), "Diffusion") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Diffusion");
  }
  if (std::find(pProc.begin(), pProc.end(), "Reaction") != pProc.end())
  {
    pIO->writeLogOutput("Initializing Reaction");
  }
  

  cout<<"Making a map" << endl;  

  /// Get the processes read and create them
  map< string, vector<double> > tempMap = pParameters->getProcesses();

  // Contruct the process
  int procsCounter = 0;
  map< string, vector<double> >::iterator mapIt = tempMap.begin();

  cout<<"Creating factory processes" << endl;
  for (; mapIt != tempMap.end(); mapIt++ )
  {
    Process* proc = FactoryProcess::createProcess( mapIt->first );
    if ( proc )
      m_vProcesses.push_back( proc );
    else 
    {
      pErrorHandler->error_simple_msg("Unknown process->" + mapIt->first );
      EXIT;
    }
  }

  if ( m_vProcesses.empty() ){
    pErrorHandler->error_simple_msg("No processes found.");
    EXIT;
    }

  /// First the processes that participate in the simulation
  /// that were read from the file input and the I/O functionality
    m_vProcesses[0]->setInstance( this );
    m_vProcesses[0]->activeSites( pLattice );
    m_vProcesses[0]->setProcessMap( &m_processMap );
}

void Apothesis::exec()
{
  ///Perform the number of KMC steps read from the input.
  int iterations = pParameters->getIterations();

  if ( iterations == 0)
  {
    pErrorHandler->error_simple_msg("Zero iterations found.");
    EXIT;
  }
  else
  {
    cout <<m_vProcesses[0]->getName()<< " process is being started..." << endl;
    pIO->writeLogOutput("Running " + to_string(iterations) + " iterations");
  }
  
    
  for ( int i = 0; i< iterations; i++)
  {
    
    ///Here a process should be selected in random. We have only one for now...

    /// Print to output
    pIO->writeLogOutput( "Time step: " + to_string( i ) );


    /// Get the active sites of the process
    list<Site*> lAdsList = m_vProcesses[ 0]->getActiveList();

    /// Check if there are available sites that it can be performed
    if (lAdsList.size() == 0)
    {
      cout << "No more "<<m_vProcesses[0]->getName()<< " site is available. Exiting..." << endl;
      pErrorHandler->error_simple_msg( "No "+ m_vProcesses[0]->getName() + " site is available. Last time step: " + to_string( i ) );
      EXIT;
    }

    /// Select randomly a site
    m_vProcesses[ 0 ]->selectSite();
    /// Perform the process
    m_vProcesses[ 0]->perform();

    /// The frequency that the various information are written in the file
    /// must befined by the user. Fix it ...
    /// The user should also check if the messages are written on the terminal or not.
    pIO->writeLogOutput( m_vProcesses[ 0]->getName() + " " );
    pIO->writeLatticeHeights();

    if(i==0){
      cout << m_vProcesses[0]->getName() << " process is being performed..." << endl;
    }
    }

  }

  void Apothesis::addProcess(string process)
  {
    m_processes.push_back(process);
  }

  vector<Species *> Apothesis::getSpecies()
  {
    return m_species;
  }

  void Apothesis::logSuccessfulRead(bool read, string parameter)
  {
    if(!pIO->outputOpen())
    {
      pIO->openOutputFile("Output");
    }
    
    read ? pIO->writeLogOutput("Reading "  + parameter) 
    :  pIO->writeLogOutput("Can't find "  + parameter) ; //pErrorHandler-> error_simple_msg("No " + parameter + " found in input file");
  }

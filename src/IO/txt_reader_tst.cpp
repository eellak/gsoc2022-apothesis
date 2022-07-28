#include "txt_reader.h"

TxtReader::TxtReader(string inputPath):m_inputPath(inputPath),
                                       m_sBuildKey("build_lattice"),
                                       m_sReadKey("read_lattice"),
                                       m_sNSpeciesKey("nspecies"),
                                       m_sNProcKey("nprocesses"),
                                       m_sPressureKey("pressure"),
                                       m_sTemperatureKey("temperature"),
                                       m_sTimeKey("time"),
                                       m_ssiteKey("*"),
                                       m_sCommentLine("#")
{}

void TxtReader::parseFile(){

    //list<string> lKeywords{m_sBuildKey, m_sReadKey, m_sNSpeciesKey, m_sNProcKey, m_sPressureKey, m_sTemperatureKey, m_sTimeKey, m_sCommentLine};


    openInputFile(m_inputPath);
    vector<string> fLines=inputFileLines();

    for(int i=0; i < fLines.size(); ++i) {
        string currentLine;
        currentLine=fLines[i];

        // split the line using space separator and store them to a vector of tokens
        vector<string> vsTokens;
        vsTokens = split(currentLine, string(" "));

        if (vsTokens[0].compare(m_sBuildKey) == 0)
        {
            m_fsetLattice(vsTokens);
        }
        if (vsTokens[0].compare(m_sReadKey) == 0)
        {
            //Cml reader constructor
            //xyz reader constructor
        }

        if(vsTokens[0].compare(m_sNSpeciesKey)==0){
            vector<string> specLines;
            int k=i+1;
            for(int j=k; j < k+toInt(vsTokens[1]); j++){
                specLines.push_back(fLines[j]);
            }
            m_fsetSpecies(specLines);
        }

        if(vsTokens[0].compare(m_sNProcKey)==0){
            vector<string> procLines;
            int k=i+1;
            for(int j=k; j < k+toInt(vsTokens[1]); j++){
                procLines.push_back(fLines[j]);
            }
            m_fsetProcesses(procLines);
        }

        if (vsTokens[0].compare(m_sPressureKey) == 0)
        {
            m_fsetPressure(vsTokens[1]);
        }

        if (vsTokens[0].compare(m_sTemperatureKey) == 0)
        {
            m_fsetTemperature(vsTokens[1]);
        }

        if (vsTokens[0].compare(m_sTimeKey) == 0)
        {
            m_fsetTime(vsTokens[1]);
        }

        if (vsTokens[0].compare(m_sDebugKey) == 0)
        {
            m_fsetDebugMode(vsTokens[1]);
        }

    }

}

void TxtReader::openInputFile(string path){
    m_inputFile.open(path, ios::in);

    if (!m_inputFile.is_open())
    {
      cout <<"Cannot open file input.kmc" <<endl;
      EXIT;
    }
}

vector<string> TxtReader::inputFileLines()
{

  string line;
  vector<string> lines;
  while (getline(m_inputFile, line))
  {

    // Remove any tabs, weird spaces etc.
    line = simplified(line);

    // split the line using space separator and store them to a vector of tokens
    vector<string> vsTokens;
    vsTokens = split(line, string(" "));

    //We do not care about empty lines
    if (vsTokens.size() == 0)
      continue;

    //We do not care about comments.
    if (startsWith(vsTokens[0], m_sCommentLine))
      continue;

    //Fill the vector
    lines.push_back(line);

  }
  return lines;
}

void TxtReader::m_fsetLattice(vector<string> vsTokens){

    std::cout << "lattice type : "<<vsTokens[1]<< std::endl;

    if (isNumber(vsTokens[2]))
    {
     std::cout << "lattice x : "<<toInt(vsTokens[2]) << std::endl;

    }
    else
    {
      std::cout << "The x dimension of lattice is not a number." << std::endl;
      EXIT;
    }

    if (isNumber(vsTokens[3]))
    {

      std::cout << "lattice y : "<<toInt(vsTokens[3]) << std::endl;
    }
    else
    {
      std::cout <<"The y dimension of lattice is not a number."<< endl;
      EXIT;
    }

    if (isNumber(vsTokens[4]))
    {
       std::cout << "lattice initial height : "<<toInt(vsTokens[4]) << std::endl;
    }
    else
    {
      std::cout<<"The height must be a  number." <<std::endl;
      EXIT;
    }
}

void TxtReader::m_fsetSpecies(vector<string> lines){

    foreach(string line,lines){
        vector<string> vsTokens;
        vsTokens = split(line, string(" "));
        if(vsTokens.size() < 2){
             std::cout << "Missing input fields for species "<< vsTokens[0] << std::endl;
        }
        else {
            if(isNumber(vsTokens[1])){
                m_mSpecies.insert({vsTokens[0],toDouble(vsTokens[1])});
            }
            else{
               std::cout << "Missing mw for species "<< vsTokens[0] << "is not a number"<< std::endl;
            }
        }

    }

}

void TxtReader::m_fsetProcesses(vector<string> lines){
    int procId=0;
    foreach(string line,lines){
        cout << line << endl;
        m_fidentifyProcess(line,procId);
        procId++;
    }
}

void TxtReader::m_fidentifyProcess(string processKey, int id){
    vector<string> process = split(processKey, string(","));
    vector<string> reactants, products;

    tie(reactants,products)=m_fsplitReactionKey(process[0]);
    vector<string> species=m_fprocSpecies(reactants,products);
    vector<double> energetics= m_fprocEnergetics(process[1]);
    vector<double> stoichiometry;

    string procName;
    if(m_bisAdsorption(reactants)){
        procName="Adsoption"+to_string(id);
    }else if(m_bisDesorption(products)){
        procName="Desorption"+to_string(id);
    }else if(m_bisDiffusion(reactants,products)){
        procName="Diffusion"+to_string(id);
    }else{
        procName="Reaction"+to_string(id);
        stoichiometry=m_fprocStoichiometry(reactants,products);
        foreach (double num, stoichiometry){
            cout << num <<endl;
        }
    }
    m_fsetProcInfo(procName,species,energetics,stoichiometry);
}

void TxtReader::m_fsetProcInfo(string procName, vector<string> species, vector<double> energetics, vector<double> stoichiometry){

    m_mProcSpecies.insert({procName,species});
    m_mProcEnergetics.insert({procName,energetics});
    if(!stoichiometry.empty())
        m_mProcStoichiometry.insert({procName,stoichiometry});
}


vector<string> TxtReader::m_fprocSpecies(vector<string> reactants, vector<string> products){
    vector<string> species;
    for(const auto& [key,value]: m_mSpecies){
         //std::cout << key << '\n';
         foreach(string reactant, reactants){
            if(contains(reactant,key)){
                if (std::find(species.begin(), species.end(), key) == species.end()) {
                  species.push_back(key);
                }
            }
         }
         foreach(string product, products){
            if(contains(product,key)){
                if (std::find(species.begin(), species.end(), key) == species.end()) {
                  species.push_back(key);
                }
            }
         }
    }

    return species;
}

vector<double> TxtReader::m_fprocStoichiometry(vector<string> reactants, vector<string> products){
    vector<double> stoichiometry;

    for(const auto& [key,value]: m_mSpecies){
         //std::cout << key << '\n';
         foreach(string reactant, reactants){
            if(contains(reactant,key)){
                eraseSubstring(reactant, key);
                cout << reactant << endl;
                if (isNumber(reactant))
                    stoichiometry.push_back(toDouble(reactant));
            }
         }
         foreach(string product, products){
             if(contains(product,key)){
                eraseSubstring(product, key);
                cout << product << endl;
                if (isNumber(product))
                    stoichiometry.push_back(toDouble(product));
             }
         }
    }

    return stoichiometry;
}


void TxtReader::eraseSubstring(string & mainStr, const string & toErase)
{
    size_t pos = string::npos;
    // Search for the substring in string in a loop untill nothing is found
    while ((pos  = mainStr.find(toErase) )!= string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

vector<double> TxtReader::m_fprocEnergetics(string energeticsKey){
    vector<double> vdEnergetics;
    vector<string> vsEnergetics = split(energeticsKey, string(" "));
    foreach (string str, vsEnergetics) {
        if (isNumber(str)){
            vdEnergetics.push_back(toDouble(str));
        }
    }
    return vdEnergetics;
}

pair<vector<string>,vector<string>> TxtReader::m_fsplitReactionKey(string reactionKey){
    vector<string> reaction = split(reactionKey, string("->"));
    vector<string> reactants=split(reaction[0],"+");
    vector<string> products=split(reaction[1],"+");
    return make_pair(reactants,products);
}

bool TxtReader::m_bisAdsorption(vector<string> reactants){
    return (reactants.size()>1) ? (simplified(reactants[0])==m_ssiteKey || simplified(reactants[1])==m_ssiteKey): false ;
}

bool TxtReader::m_bisDesorption(vector<string> products){
    return (products.size()>1) ? (simplified(products[0])==m_ssiteKey || simplified(products[1])==m_ssiteKey): false ;
}

bool TxtReader::m_bisDiffusion(vector<string> reactants, vector<string> products){
    return (reactants.size()==1 && products.size()==1 && contains(reactants[0],m_ssiteKey) && contains(products[0],m_ssiteKey)) ;
}

void TxtReader::m_fsetTime(string time){
    if (isNumber(time))
    {
      std::cout << "Time "<<toDouble(time)<<std::endl;
    }
    else
    {
      std::cout << "Could not read number of KMC simulation time from input file. Is it a number?"<<std::endl;
      EXIT;
    }
}

void TxtReader::m_fsetPressure(string pressure){
    if (isNumber(pressure))
    {
      std::cout << "Pressure : "<<toDouble(pressure) << std::endl;
    }
    else
    {
      std::cout << "Could not read pressure from input file. Is it a number?"<<std::endl ;
      EXIT;
    }
}

void TxtReader::m_fsetTemperature(string temperature){
    if (isNumber(temperature))
    {
      std::cout << "Temperature : "<<toDouble(temperature) << std::endl;
    }
    else
    {
      std::cout << "Could not read temperature from input file. Is it a number?"<<std::endl ;
      EXIT;
    }
}

void TxtReader::m_fsetDebugMode(string debugMode){
    cout << debugMode << endl;
    if (contains(debugMode,"on", Insensitive) || contains(debugMode,"off", Insensitive))
    {
      std::cout << "Debug mode "<<debugMode<<std::endl;
    }
    else
    {
      std::cout << "Could not read bebug mode."<<std::endl;
      EXIT;
    }
}

string TxtReader::simplified(string str)
{
  string s;
  bool is_white = false;
  bool was_white = false;
  bool append = false;
  size_t n = str.size();
  size_t nm = n - 1;

  for (size_t i = 0; i < n; i++)
  {
    char c;
    c = str[i];
    switch (c)
    {
    case ' ':
    case '\n':
    case '\t':
    case '\v':
    case '\f':
    case '\r':
      is_white = true;
      c = ' ';
      append = (was_white || i == 0 || i == nm) ? false : true;
      break;
    default:
      is_white = false;
      append = true;
      break;
    };
    if (append)
      s.push_back(c);
    was_white = is_white;
  }
  return s;
}

bool TxtReader::isNumber(string str)
{

  size_t n = str.size();

  bool bPlusFound = false;
  bool bPlusFound2 = false;
  bool bMinusFound = false;
  bool bMinusFound2 = false;
  bool be = false;
  bool bE = false;
  bool bd = false;
  bool bD = false;
  bool bf = false;
  bool bF = false;
  bool bDotFound = false;

  for (size_t i = 0; i < n; i++)
  {
    char c;
    c = str[i];
    switch (c)
    {
    case '+':
      if (i == 0)
      {
        bPlusFound = true;
        break;
      }

      if (i != 0 && (bPlusFound || bMinusFound) && !bD && !bd && !be && !bE && !bf && !bF)
        return false;

      if ((bD || bd || be || bE || bf || bF) && (!bPlusFound2))
        bPlusFound2 = true;
      else if ((bD || bd || be || bE || bf || bF) && bPlusFound2)
        return false;
      else if ((bD || bd || be || bE || bf || bF) && bMinusFound2)
        return false;

      break;
    case '-':
      if (i == 0)
      {
        bMinusFound = true;
        break;
      }

      if (i != 0 && (bPlusFound || bMinusFound) && !bD && !bd && !be && !bE && !bf && !bF)
        return false;

      if ((bD || bd || be || bE || bf || bF) && (!bPlusFound2))
        bPlusFound2 = true;
      else if ((bD || bd || be || bE || bf || bF) && bPlusFound2)
        return false;
      else if ((bD || bd || be || bE || bf || bF) && bMinusFound2)
        return false;

      break;
    case 'e':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        be = true;
      else
        return false;
      break;
    case 'E':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bE = true;
      else
        return false;
      break;
    case 'D':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bd = true;
      else
        return false;
      break;
    case 'd':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bD = true;
      else
        return false;
      break;
    case 'F':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bf = true;
      else
        return false;
      break;
    case 'f':
      if (!bD && !bd && !be && !bE && !bf && !bF)
        bF = true;
      else
        return false;
      break;
    case '.':
      if (!bDotFound && !bD && !bd && !be && !bE && !bf && !bF)
      {
        bDotFound = true;
      }
      else if (!bDotFound && (bD || bd || be || bE || bf || bF))
        return false;
      else if (bDotFound)
        return false;
      break;
    default:
      if (!isdigit(c))
        return false;
    };
  }
  return true;
}

double TxtReader::toDouble(string str)
{
  if (isNumber(str))
    return stod(str);
  else
    return 0;
}

vector<string> TxtReader::split(string str, string delim)
{
  vector<string> r;
  size_t prev = 0, pos = 0;
  do
  {
    pos = str.find(delim, prev);
    if (pos == string::npos)
      pos = str.length();
    string token = str.substr(prev, pos - prev);
    if (!token.empty())
      r.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());

  return r;
}

int TxtReader::toInt(string str)
{
  if (isNumber(str))
    return stof(str);
  else
    return 0;
}

bool TxtReader::contains(string str1, string str2, CASE cas)
{
  switch (cas)
  {
  case Insensitive:
  {
    char c;
    char csmall;
    string smstr, sstr;
    size_t n = str1.size();

    for (size_t i = 0; i < n; i++)
    {
      c = str1[i];
      csmall = (char)tolower(c);
      smstr.push_back(csmall);
    }

    for (size_t i = 0; i < str2.size(); i++)
    {
      c = str2[i];
      csmall = (char)tolower(c);
      sstr.push_back(csmall);
    }

    if (smstr.find(sstr) != string::npos)
      return true;
    return false;
  }
  default:
    if (str1.find(str2) != string::npos)
      return true;
    return false;
  }
}

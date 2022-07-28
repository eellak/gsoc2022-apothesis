#ifndef TXTREADER_H
#define TXTREADER_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>

#define EXIT { printf("Apothesis terminated. \n"); exit( EXIT_FAILURE ); }

using namespace std;

class TxtReader
{
public:
    enum CASE{ Sensitive, Insensitive };

    explicit TxtReader(string inputPath);
    void parseFile();

    /// Opens the input file.
    void openInputFile(string);

    /// Returns input file lines without empty lines and comments.
    vector<string> inputFileLines();

    /// Converts a string to double.
    double toDouble(string);

    /// Converts a string to int.
    int toInt(string);

    /// Check if a string contains another string and remove it.
    void eraseSubstring(string &, const string &);


    /// Check if a string contains another string. TODO: This should be transferred to a generic string class).
    bool contains(string, string, CASE cas = Insensitive);

    /// Splits a string to a vector of strings. TODO: This should be transferred to a generic string class).
    vector<string> split(string , string);

    /// Returns if the given string in number. TODO: This should be transferred to a generic string class).
    bool isNumber(string);

    /// Given a string returns a string with all the delimeters replaced. TODO: This should be transferred to a generic string class).
    string simplified( string );

    /// Checks if a file exists.
    bool exists(const string& s);

    inline bool startsWith( string str, string substr )
      {
      if ( str.find( substr ) != string::npos )
        return true;
      else
        return false;
      }




private:
    ///Path of input.kmc
    string m_inputPath;

    ///Filestream for input.kmc
    ifstream m_inputFile;

    /// Keywords:
    /// Build lattice keyword
    string m_sBuildKey;

    /// Read lattice from file keyword
    string m_sReadKey;

    /// Species keyword
    string m_sNSpeciesKey;

    /// Process keyword
    string m_sNProcKey;

    /// Temperature keyword
    string m_sTemperatureKey;

    /// Pressure keyword
    string m_sPressureKey;

    ///  Num of iterations keyword.
    string m_sTimeKey;

    ///  Debug mode keyword.
    string m_sDebugKey;

    /// Reaction site key
    string m_ssiteKey;

    /// Comment
    string m_sCommentLine;

    /// Input variables:
    /// Simulation temperature
    double m_dTemperature;

    /// Simulation time
    double m_dTime;

    /// Simulation pressure
    double m_dPressure;

    /// Debug mode
    string m_sDebugMode;

    /// Species representation in a map species name key and mw as value
    map<string,double> m_mSpecies;

    /// Process maps
    map<string,vector<string>> m_mProcSpecies;
    map<string,vector<double>> m_mProcEnergetics;
    map<string,vector<double>> m_mProcStoichiometry;

    /// Lattice


    /// Set lattice info
    void m_fsetLattice(vector<string>);

    /// Set species info
    void m_fsetSpecies(vector<string>);

    /// Set processes info
    void m_fsetProcesses(vector<string>);

    /// Set simulation time
    void m_fsetTime(string);

    /// Set pressure
    void m_fsetPressure(string );

    /// Set temperature
    void m_fsetTemperature(string);

    /// Set debug mode
    void m_fsetDebugMode(string);

    /// Get left part of process keyword and identify the type of process
    void m_fidentifyProcess(string,int);

    /// Set process species, energetics & stoichiometry map
    void m_fsetProcInfo(string, vector<string>, vector<double>, vector<double>);

    /// Get right part of process keyword and identify the process energetics
    vector<double> m_fprocEnergetics(string);

    /// Get left part of process keyword and identify the process stoichiometry
    vector<double> m_fprocStoichiometry(vector<string>, vector<string>);

    /// Get left part of process keyword and identify the participating species
    vector<string> m_fprocSpecies(vector<string>, vector<string>);

    /// Returns reactants and products for reaction key
    pair<vector<string>,vector<string>> m_fsplitReactionKey(string);

    ///Adsorption
    bool m_bisAdsorption(vector<string>);

    ///Desorption
    bool m_bisDesorption(vector<string>);

    ///Diffusion
    bool m_bisDiffusion(vector<string>,vector<string>);


};

#endif // TXTREADER_H

#ifndef TXTREADER_H
#define TXTREADER_H

#include <vector>
#include <iostream>
#include <string>

using namespace std;

class TxtReader
{
public:
    explicit TxtReader(string inputPath);

private:
    /// The input file path
    string m_sInputPath;
};

#endif // TXTREADER_H

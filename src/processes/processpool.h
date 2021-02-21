#ifndef PROCESSPOOL_H
#define PROCESSPOOL_H

#include <string>
#include <map>

#include "process_new.h"

using namespace std;

namespace newDesign
{

class ProcessPool
{
public:
    ProcessPool();

    inline Process_new* getProcessByName(string name) { return m_mapProcs[ name ]; }
    inline Process_new* getProcessByID( int id ){ return m_mapProcsIDs[ id ]; }

    void addProcess( string name, Process_new* proc);
    void addProcess( int id, Process_new* proc);
    inline int getProcessNum() { return m_mapProcs.size(); }

private:
    map<string, Process_new*> m_mapProcs;
    map<int, Process_new*> m_mapProcsIDs;
};

}

#endif // PROCESSPOOL_H

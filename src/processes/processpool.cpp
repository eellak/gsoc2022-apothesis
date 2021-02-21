#include "processpool.h"

namespace newDesign
{

ProcessPool::ProcessPool(){}

void ProcessPool::addProcess( string name, Process_new* proc)
{
    pair<string, Process_new* > p;
    p.first = name;
    p.second = proc;

    m_mapProcs.insert( p );
}

void ProcessPool::addProcess( int id, Process_new* proc)
{
    pair<int, Process_new* > p;
    p.first = id;
    p.second = proc;

    m_mapProcsIDs.insert( p );
}




}

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

#include "factory_process.h"
#include "process.h"

using namespace MicroProcesses;

FactoryProcess::FactoryProcess(){ }
FactoryProcess::~FactoryProcess(){ }

void FactoryProcess::registerThis( const string& s, AbstractProcess* proc ){
    getTable()[ s] = proc;
}

MicroProcesses::Process* FactoryProcess::createProcess( const string& s){

    map< string, AbstractProcess* >::iterator it = getTable().find( s);

    if  ( it != getTable().end() )
        return it->second->create();
    else
        return (MicroProcesses::Process*)( 0 );
}

std::map< string, AbstractProcess* >& FactoryProcess::getTable(){
    static std::map<std::string, AbstractProcess*> table;
    return table;
}



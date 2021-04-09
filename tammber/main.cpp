/*
   Copyright (c) 2016, Los Alamos National Security, LLC
   All rights reserved.
   Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

   Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
   1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
   2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
   3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */



#include <iostream>
#include <algorithm>

#include "CustomTypedefs.hpp"
#include "Constants.hpp"
#include "XYZSystem.hpp"
#include "LAMMPSSystem.hpp"
#include "LAMMPSEngine.hpp"
#include "Task.hpp"
#include "Graph.hpp"
#include "TaskManager.hpp"
#include "Worker.hpp"

#include "DDS3.hpp"
#include "Log.hpp"

#include <mpi.h>
#include <vector>
#include <ostream>
#include <streambuf>
#include <sstream>
#include <unistd.h>
#include <limits.h>


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/lexical_cast.hpp>

#include <chrono>
#include <thread>


#include "HCDS.hpp"
#include "LocalStore.hpp"
#include "DDS.hpp"

//This file contains the type definitions. Edit this file to change to LAMMPS classes
#include "CustomTypedefs.hpp"

#include "AbstractSystem.hpp"

#include "PullMMbuilder.hpp"
#include "PullWorkManager.hpp"
#include "RankPlacement.hpp"

#include "Log.hpp"


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/program_options.hpp>

bool opt_parser(int argc, char * argv[]);

int main(int argc, char * argv[]) {
	pid_t pid=getpid();

	// Check for special flags
	if(opt_parser(argc, argv)) return 0;


	MPI_Init(&argc, &argv);

	int rank;
	int parent;
	int nranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	std::cout << "TAMMBER-main\n";

	// Create empty property tree object
	boost::property_tree::ptree tree;

	// Parse the XML into the property tree.
	boost::property_tree::read_xml("./input/ps-config.xml", tree,boost::property_tree::xml_parser::no_comments | boost::property_tree::xml_parser::trim_whitespace);
	//boost::property_tree::read_info("./input/ps-config.info", tree);


	std::string dbRoot=tree.get<std::string>("Configuration.DB.RootDirectory","./");

	#ifdef USE_BOOST_LOG
	initLog(rank,1);
	#endif

	int nNodes=boost::lexical_cast<int>(argv[1]);
	int nSlotsPerNode=boost::lexical_cast<int>(argv[2]);
	int nSlotsPerWorker=boost::lexical_cast<int>(argv[3]);
	int nWorkManagers=boost::lexical_cast<int>(argv[4]);
	bool padNodes=boost::lexical_cast<bool>(argv[5]);
	int mseed=boost::lexical_cast<int>(argv[6]);

	RankPlacer placer(nNodes,nSlotsPerNode,nSlotsPerWorker,nWorkManagers,padNodes,mseed);

	std::map<int,std::tuple<int,int,int,std::string,int> > &layout=placer.layout;
	std::map<int,int> myChildren;
	std::set<int> wgroups;

	int r,t,p,seed;
	std::string s;

	for(auto it=layout.begin(); it!=layout.end(); it++) {
		r=std::get<0>(it->second);
		t=std::get<1>(it->second);
		p=std::get<2>(it->second);
		s=std::get<3>(it->second);
		seed=std::get<4>(it->second);

		if(p==rank) {
			if(myChildren.count(t)==0) {
				myChildren[t]=r;
			}
			else{
				myChildren[t]=std::min(r,myChildren[t]);
			}
		}

		if(s=="Worker") {
			wgroups.insert(seed);
		}

	}

	int myGroup=std::get<1>(layout[rank]);
	int myParent=std::get<2>(layout[rank]);
	std::string myType=std::get<3>(layout[rank]);
	int mySeed=std::get<4>(layout[rank]);


	char hostname[_POSIX_HOST_NAME_MAX];
	gethostname(hostname, _POSIX_HOST_NAME_MAX);

	std::cout<<"LAYOUT: "<<pid<<" "<<myGroup<<" "<<myParent<<" "<<myType<<" "<<hostname<<std::endl;
	std::set<int> dbKeys;
	dbKeys.insert(0);
	dbKeys.insert(1);

	std::unordered_map< unsigned int, std::pair<bool,bool> > dbAttributesMin;
	dbAttributesMin[0]=std::make_pair(false,false);
	dbAttributesMin[1]=std::make_pair(false,false);


	MPI_Comm dbComm;
	MPI_Comm workComm;
	MPI_Comm localComm;

	//form local groups
	MPI_Comm_split(MPI_COMM_WORLD, myGroup, 1, &localComm);

	std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
	std::chrono::minutes runTime=std::chrono::minutes( tree.get<unsigned>("Configuration.RunTime",1000000) );



	// Everyone waits for Splicer
	MPI_Barrier(MPI_COMM_WORLD);

	uint64_t maxDataSize=100000000;
	uint64_t maxCacheSize=1000000000;

	//empty placeholder to enforce node placement
	if(myType=="Nothing") {
		//db-comm
		MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,0,&dbComm);
		//work comm
		MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,0,&workComm);
	}


	if(myType=="Splicer") {

		//db-comm
		MPI_Comm_split(MPI_COMM_WORLD,1,2,&dbComm);
		//work comm
		MPI_Comm_split(MPI_COMM_WORLD,1,0,&workComm);




		//form intercommunicator with child

		/*
		        MPI_Comm worker;
		   assert(myChildren.size()==1);
		   for(auto it=myChildren.begin(); it!=myChildren.end(); it++) {
		   std::cout<<it->second<<" "<<it->first<<std::endl;
		   MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, it->second, it->first, &worker);
		   }
		 */


		std::set<int> children;
		int nrankWork;
		MPI_Comm_size(workComm, &nrankWork);
		for(int i=1; i<nrankWork; i++) {
			std::cout<<"CHILDREN: "<<i<<std::endl;
			children.insert(i);
		}
		HDDS3<STLLocalDataStore> minimaStore(dbComm,1,dbRoot+"./db0/","min",maxDataSize,maxCacheSize,dbKeys,dbAttributesMin,true);
		//ParSpliceSplicer splicer(workComm,&minimaStore,children,wgroups.size(),tree);
		TammberModelBuilder mmbuilder(workComm,&minimaStore,children,tree);
		std::cout<<"TammberModelBuilder is established"<<std::endl;
		std::cout<<"SERVER() CALLED BY.. "<<myType<<std::endl<<std::flush;
		mmbuilder.server();
	}


	if(myType=="WorkManager") {

		//db-comm
		MPI_Comm_split(MPI_COMM_WORLD,1,2,&dbComm);
		//work comm
		MPI_Comm_split(MPI_COMM_WORLD,1,1,&workComm);


		//form intercommunicators with children
		std::set<MPI_Comm> workers;
		for(auto it=myChildren.begin(); it!=myChildren.end(); it++) {
			MPI_Comm interComm;
			std::cout<<it->second<<" "<<it->first<<std::endl;
			MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, it->second, it->first, &interComm);
			workers.insert(interComm);
		}

		STLLocalDataStore qsdStore;
		qsdStore.initialize("","");
		qsdStore.setMaximumSize(1000000000);
		HDDS3<STLLocalDataStore> minimaStore(dbComm,1,dbRoot+"./db0/","min",maxDataSize,maxCacheSize,dbKeys,dbAttributesMin,true);
		WorkManagerType workManager(workComm, 0, &minimaStore, &qsdStore,  tree, workers);
		workManager.server();
		std::cout<<"WorkManager is done"<<std::endl;
	}


	if(myType=="PersistentDB") {

		//db-comm
		MPI_Comm_split(MPI_COMM_WORLD,1,0,&dbComm);
		//work comm
		MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,0,&workComm);



		uint64_t maxDataSize=100000000;
		uint64_t maxCacheSize=0;

		//CompressedPersistentLocalStore
		//PersistentLocalStore
		HDDS3<PersistentLocalStore> minimaStore(dbComm,-1,dbRoot+"./db0/","min",maxDataSize,maxCacheSize,dbKeys,dbAttributesMin,false);
		//HDDS3<STLLocalDataStore> minimaStore(dbComm,-1,dbRoot+"./db0/","min",maxDataSize,maxCacheSize,dbKeys,dbAttributesMin,false);

		std::chrono::seconds syncDelay( tree.get<unsigned>("Configuration.DB.SyncDelay",1) );

		std::chrono::high_resolution_clock::time_point lastSync=std::chrono::high_resolution_clock::now();
		while(true) {
			minimaStore.singleServe();
			//minimaStore.printStatus();

			if(std::chrono::high_resolution_clock::now() - lastSync > syncDelay  ) {
				LOGGER("DB SYNC")
				minimaStore.printStatus();

				minimaStore.sync();
				lastSync=std::chrono::high_resolution_clock::now();
			}

			if(std::chrono::high_resolution_clock::now() - start > runTime  ) {
				std::cout<<"DB is going down"<<std::endl;
				minimaStore.sync();
				minimaStore.cancelCommunications();
				break;
			}

		}
		std::cout<<"PersistentDB is done"<<std::endl;
	}


	if(myType=="InMemoryDB") {
		//db-comm
		MPI_Comm_split(MPI_COMM_WORLD,1,1,&dbComm);
		//work comm
		MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,0,&workComm);

		uint64_t sharedBufferSize=tree.get<uint64_t>("Configuration.DB.SharedCacheSize",1000000000);

		HDDS3<STLLocalDataStore> minimaStore(dbComm,0,dbRoot+"./db0/","min",maxDataSize,sharedBufferSize,dbKeys,dbAttributesMin,true);
		std::chrono::seconds syncDelay( tree.get<unsigned>("Configuration.DB.SyncDelay",1) );

		std::chrono::high_resolution_clock::time_point lastSync=std::chrono::high_resolution_clock::now();
		while(true) {
			minimaStore.singleServe();
			if(std::chrono::high_resolution_clock::now() - lastSync > syncDelay  ) {
				LOGGER("DB SYNC")
				minimaStore.printStatus();
				//minimaStore.sync();
				lastSync=std::chrono::high_resolution_clock::now();
			}
			if(std::chrono::high_resolution_clock::now() - start > runTime  ) {
				LOGGERA("DB is going down")
				minimaStore.sync();
				minimaStore.cancelCommunications();
				break;
			}
		}
		std::cout<<"InMemoryDB is done"<<std::endl;

	}

	if(myType=="Worker") {
		//db-comm
		MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,0,&dbComm);
		//work comm
		MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,0,&workComm);

		int nloc;
		MPI_Comm_size(localComm, &nloc);


		std::cout<<myParent<<" "<<myGroup<<std::endl;

		MPI_Comm interComm;
		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, myParent, myGroup, &interComm);


		worker(localComm,interComm,mySeed);
		std::cout<<"Worker is done"<<std::endl;
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

}

bool opt_parser(int argc, char * argv[]){

	// Quick Return if there are no command-line options
	if(argc < 2) return 0;

	namespace po = boost::program_options;
	po::options_description desc("Allowed command-line options");
	desc.add_options()
	// First parameter describes option name/short name
	// The second (if there are 3) is parameter to option
	// The last is description
	        ("help,h",  "print usage message")
	        ("info,d",  "print description of ParSplice")
	;

	po::variables_map vm;
	po::store(parse_command_line(argc, argv, desc), vm);

	// "--help" Option should just print information
	// about the ParSplice code..
	if (vm.count("help")) {
		std::cout << "\n" << desc << "\n";
		return 1;
	}

	// "--info" Option should just print information
	// about the ParSplice code..
	if (vm.count("info")) {
		std::cout <<"\nThis is ParSplice (https://gitlab.com/exaalt/parsplice).\n"
		          <<"\nParSplice is the Accelerated Molecular Dynamics (AMD) engine for the\n"
		          <<"EXAALT (Exascale Atomistics for Accuracy, Length and Time) Framework.\n"
		          <<"EXAALT includes the ParSplice, LAMMPS, and LATTE codes as components,\n"
		          <<"and is funded by the US DOE Exascale Computing Project (ECP).\n";
		std::cout <<"\n";
		return 1;
	}
	return 0;
}

#include <iostream>
#include <algorithm>

#include "CustomTypedefs.hpp"
#include "Constants.hpp"
#include "XYZSystem.hpp"
#include "LAMMPSSystem.hpp"
#include "LAMMPSEngine.hpp"
#include "Task.hpp"
#include "Log.hpp"
#include "Graph.hpp"
#include "TaskManager.hpp"
#include "Worker.hpp"

#include <mpi.h>
#include <vector>
#include <ostream>
#include <streambuf>
#include <sstream>


#include <chrono>
#include <thread>

#include "HCDS.hpp"
#include "LocalStore.hpp"
#include "DDS.hpp"
#include "AbstractSystem.hpp"
#include "ModelWrapper.hpp"


#include <boost/filesystem.hpp>

int main(int argc, char * argv[]) {
	MPI_Init(&argc, &argv);
	int rank;
	int nranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	TaskMapperType mapper;

	std::cout << "TAMMBER-neb\n";


	MPI_Comm localComm, workerComm;

	if(rank==0) {

		MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &localComm);

		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, 1, 1, &workerComm);

		DriverHandleType handle(workerComm);

		ModelWrapper mmbuilder(config);
		std::cout<<"TammberModel loaded"<<std::endl;

		GenericTask label,neb;

		label.type=mapper.type("TASK_LABEL");
		label.flavor=1;

		neb.type=mapper.type("TASK_NEB");
		neb.flavor=1;

		PersistentLocalStore minimaStore;
		minimaStore.initialize("./db0/", "min");
		minimaStore.createDatabase(0, false, false);

		bool found_redo= boost::filesystem::exists("./RedoNEBS.list");

		uint64_t c1,l1,c2,l2;
		Transition t;
		std::set<Transition> trans;
		trans.clear();
		//
		std::set<Transition> rw_trans;
		rw_trans.clear();

		/*
		bool found_rewrite= boost::filesystem::exists("./RewriteNEBS.xml");
		if(found_rewrite) {
			std::list<NEBPathway> paths;
		  parse_xml_path(paths); // TODO
			for(auto p : paths) {
				mmbuilder.add_pathway(p);
				t.first = p.InitialLabels;
				t.second = p.FinalLabels;
				rw_trans.insert(t);
			}
			LOGGERA("Implemented "<<rw_trans.size()<<" rewrite requests");
		}*/

		if(found_redo) {
			std::ifstream infile("./RedoNEBS.list");
			while(infile>>c1>>l1>>c2>>l2) {
				LOGGERA(c1<<" "<<c2)
				t.first.first=c1;
				t.first.second=l1;
				t.second.first=c2;
				t.second.second=l2;
				if(rw_trans.find(t)==rw_trans.end()) trans.insert(t);
			}
			LOGGERA("Found "<<trans.size()<<" redo requests which are not being rewritten");
		}


		if(found_redo) {
			for(auto tran: trans) {

				RawDataVector data;
				SystemType initial,final;
				Transition rl_tran;

				minimaStore.get(LOCATION_SYSTEM_MIN,tran.first.second,data);
				if(data.size()==0) continue;
				unpack(data,initial,data.size());

				label.clearInputs(); label.clearOutputs();
				insert("State",label.inputData,initial);
				handle.assign(label);
				while(not handle.probe(label)) {};
				extract("Labels",label.returns,rl_tran.first);

				minimaStore.get(LOCATION_SYSTEM_MIN,tran.second.second,data);
				if(data.size()==0) continue;
				unpack(data,final,data.size());

				label.clearInputs(); label.clearOutputs();
				insert("State",label.inputData,final);
				handle.assign(label);
				while(not handle.probe(label)) {};
				extract("Labels",label.returns,rl_tran.second);

				if(rl_tran!=tran) LOGGERA("LABEL MISMATCH"<<rl_tran.print_str());

				NEBPathway pathway;
				pathway.InitialLabels = rl_tran.first;
				pathway.FinalLabels = rl_tran.second;
				pathway.pairmap=false;

				LOGGERA("SUBMITTING PATHWAY FOR NEB "<<pathway.submit_info_str())


				neb.clearInputs(); neb.clearOutputs();
				insert("Initial",neb.inputData,initial);
				insert("Final",neb.inputData,final);
				insert("NEBPathway",neb.arguments,pathway);

				handle.assign(neb);
				while(not handle.probe(neb)) {};
				extract("NEBPathway",neb.returns,pathway);

				mmbuilder.add_pathway(pathway);
			}
			mmbuilder.save();
		} else {
			std::cout<<"RedoNEBS.list not found!"<<std::endl;
		}

		label.type=mapper.type("TASK_DIE");
		label.flavor=1;
		handle.assign(label);
	}

	else{
		MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &localComm);
		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, 0, 1, &workerComm);
		worker(localComm,workerComm,1234);
	}





	MPI_Finalize();
	return 0;

};

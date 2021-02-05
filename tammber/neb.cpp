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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/filesystem.hpp>

#include <chrono>
#include <thread>

#include "HCDS.hpp"
#include "LocalStore.hpp"
#include "DDS.hpp"
#include "AbstractSystem.hpp"
#include "DummyMMbuilder.hpp"


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>



int main(int argc, char * argv[]) {
	MPI_Init(&argc, &argv);
	int rank;
	int nranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	TaskMapperType mapper;

	std::cout << "TAMMBER-neb\n";

	// Create empty property tree object
	boost::property_tree::ptree config;
	// Parse the XML into the property tree.
	boost::property_tree::read_xml("./input/ps-config.xml", config, boost::property_tree::xml_parser::no_comments );
	std::map< std::pair<int,int>, std::map<std::string,std::string> > parameters;
	BOOST_FOREACH(boost::property_tree::ptree::value_type &v, config.get_child("Configuration.TaskParameters")) {
		boost::optional<std::string> otype= v.second.get_optional<std::string>("Task");
		if(otype) {
			std::string stype=*otype;
			boost::trim(stype);
		}
	}

	MPI_Comm localComm, workerComm;

	if(rank==0) {

		MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &localComm);

		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, 1, 1, &workerComm);

		DriverHandleType handle(workerComm);

		DummyModelBuilder mmbuilder(config);
		std::cout<<"TammberModel loaded"<<std::endl;

		GenericTask label,neb;

		label.type=mapper.type("TASK_LABEL");
		label.flavor=1;

		neb.type=mapper.type("TASK_NEB");
		neb.flavor=1;

		PersistentLocalStore minimaStore;
		minimaStore.initialize("./db0/", "min");
		minimaStore.createDatabase(0, false, false);

		std::set<uint64_t> keys = minimaStore.availableKeys(LOCATION_SYSTEM_MIN);

		bool found_redo= boost::filesystem::exists("./RedoNEBS.list");
		Transition t;
		std::set<Transition> trans;

		if(found_redo) {
			std::ifstream infile("./RedoNEBS.list");
			trans.clear();
			while(infile>>t.first.first>>t.first.second>>t.second.first>>t.second.second);
				trans.insert(tran);
		} else {
			std::cout<<"RedoNEBS.list not found!"<<std::endl;
		}

				label.clearInputs(); label.clearOutputs();
				insert("State",label.inputData,initial);
				handle.assign(label);
				while(not handle.probe(label)) {};
				extract("Labels",label.returns,rl_tran.first);

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

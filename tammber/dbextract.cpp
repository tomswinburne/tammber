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

	std::cout << "TAMMBER-dbextract\n";
	
	MPI_Comm localComm, workerComm;

	if(rank==0) {

		MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &localComm);

		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, 1, 1, &workerComm);


		DriverHandleType handle(workerComm);

		GenericTask label,write;

		label.type=mapper.type("TASK_LABEL");
		label.flavor=1;

		write.type=mapper.type("TASK_WRITE_TO_FILE");
		write.flavor=1;

		PersistentLocalStore minimaStore;
		minimaStore.initialize("./db0/", "min");
		minimaStore.createDatabase(0, false, false);

		bool readfromlist= boost::filesystem::exists("./Candidates.list");
		uint64_t flab,slab;

		std::set<uint64_t> keys = minimaStore.availableKeys(LOCATION_SYSTEM_MIN);

		if(!readfromlist) {
			std::cout<<"Candidates.list not found, extracting everything.."<<std::endl;
		} else {
			std::ifstream infile("./Candidates.list");
			keys.clear();
			while(infile>>flab>>slab) keys.insert(slab);
		}

		for(auto key : keys) {
			RawDataVector data;
			minimaStore.get(LOCATION_SYSTEM_MIN,key,data);
			std::cout<<key<<" "<<data.size()<<std::endl;
 			if(data.size()==0) continue;
			SystemType s;
			unpack(data,s,data.size());

			label.clearInputs(); label.clearOutputs();
			write.clearInputs(); write.clearOutputs();

			insert("State",label.inputData,s);
			LabelPair labels;
			handle.assign(label);
			while(not handle.probe(label)) {};
			extract("Labels",label.returns,labels);
			std::string filename="states/state-"+std::to_string(labels.first)+"-"+std::to_string(labels.second)+".dat";
			insert("Filename",write.arguments,filename);
			insert("State",write.inputData,s);
			handle.assign(write);
			while(not handle.probe(write)) {};

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

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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>

int main(int argc, char * argv[]) {
	MPI_Init(&argc, &argv);
	int rank;
	int nranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	TaskMapperType mapper;

	std::cout << "TAMMBER-md\n";

	// Create empty property tree object
	boost::property_tree::ptree config;
	// Parse the XML into the property tree.
	boost::property_tree::read_xml("./input/ps-config.xml", config, boost::property_tree::xml_parser::no_comments );

	std::string ic_string = config.get<std::string>("Configuration.InitialConfigurations");
	std::vector<std::string> initialConfigurations;
	boost::split(initialConfigurations,ic_string,boost::is_any_of(" "));



	MPI_Comm localComm, workerComm;

	if(rank==0) {

		MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &localComm);

		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, 1, 1, &workerComm);

		DriverHandleType handle(workerComm);

		GenericTask task,task_seg;
		task.type=mapper.type("TASK_INIT_MIN");
		task.flavor=1;
		task_seg.type=mapper.type("TASK_SEGMENT");
		task_seg.flavor=1;

		LOGGERA("TAMMBER-md::initializeSystems()")

		for(auto initialConfiguration : initialConfigurations) {
			boost::trim(initialConfiguration);
			if(initialConfiguration.length()==0) continue;

			SystemType system;
			std::string file_name;
			std::set<PointShiftSymmetry> self_symmetries;
			LabelPair labels;
			double energy;
			int clusters;
			std::array<double,3> position;
			TADSegment segment;
			bool ProductionRun = false; // for testing
			double temperature = config.get<double>("Configuration.MarkovModel.MinTemperature",100.0);
			segment.InitialLabels = labels;
			segment.temperature = temperature;

			/*
			This returns
		  * returns: Labels, Clusters, Position, SelfSymmetries
		  * outputData State
			*/
			task.clearInputs(); task.clearOutputs();
			insert("Filename",task.arguments,initialConfiguration);
			LOGGERA("REQUESTING "<<file_name)
			handle.assign(task);
			while(not handle.probe(task)) {};


			extract("Labels",task.returns,labels);
			extract("Energy",task.returns,energy);
			extract("Clusters",task.returns,clusters);
			extract("Position",task.returns,position);
			LOGGERA("LABELS: "<<labels.first<<" "<<labels.second<<" E:"<<energy
				<<"eV, Clusters:"<<clusters
				<<" Position:"<<position[0]<<" "<<position[1]<<" "<<position[2])

			#ifdef ISOMORPHIC
			extract("SelfSymmetries",task.returns,self_symmetries);
			LOGGERA("SelfSymmetries:")
			for(auto ss:self_symmetries) LOGGERA(ss.info_str());
			#endif

			extract("State",task.outputData,system);

			LOGGERA("State has "<<system.getNAtoms()<<" atoms!");

			task_seg.clearInputs(); task_seg.clearOutputs();

			insert("TADSegment",task_seg.arguments,segment);
			insert("ProductionRun",task_seg.arguments,ProductionRun);
			insert("Minimum",task_seg.inputData,system);

			LOGGERA("STARTING TASK_SEGMENT AT MarkovModel.MinTemperature = "<<temperature<<"K");
			LOGGERA("TAMMBER-md: SUBMITTING SEGMENT "<<segment.submit_info_str())
			handle.assign(task_seg);
			while(not handle.probe(task_seg)) {};

			extract("TADSegment",task_seg.returns,segment);
			LOGGERA("TAMMBER-md: OUTPUT SEGMENT "<<segment.info_str())
		}

		task.clear();
		task.type=mapper.type("TASK_DIE");
		task.flavor=1;
		handle.assign(task);

	}	else {
		MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &localComm);
		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, 0, 1, &workerComm);
		worker(localComm,workerComm,1234);
	}





	MPI_Finalize();
	return 0;

};

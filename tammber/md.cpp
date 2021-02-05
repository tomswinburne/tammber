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


	MPI_Comm localComm, workerComm;

	if(rank==0) {

		MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &localComm);

		MPI_Intercomm_create( localComm, 0, MPI_COMM_WORLD, 1, 1, &workerComm);

		DriverHandleType handle(workerComm);

		GenericTask task;
		SystemType system;
		PointShiftSymmetry self_symmetries;
		LabelPair labels;
		double energy;
		int clusters;
		std::array<double,3> position;
		/*
		This returns
	  * returns: Labels, Clusters, Position, SelfSymmetries
	  * outputData State
		*/
		task.type=mapper.type("TASK_INIT_MIN");
		task.flavor=1;

		task.clearInputs(); task.clearOutputs();
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
		double temperature = config.get<double>("Configuration.MarkovModel.MinTemperature",100.0);
		LOGGERA("STARTING TASK_SEGMENT AT MarkovModel.MinTemperature = "<<temperature<<"K");
		task.clear();

		task.type=mapper.type("TASK_SEGMENT");
		task.flavor=1;

		TADSegment segment;
		segment.InitialLabels = labels;
		segment.temperature = temperature;
		bool ProductionRun = false; // for testing
		insert("TADSegment",task.arguments,segment);
		insert("ProductionRun",task.arguments,ProductionRun);
		insert("Minimium",task.inputData,system);

		LOGGERA("TAMMBER-md: SUBMITTING SEGMENT "<<segment.submit_info_str())
		handle.assign(task);
		while(not handle.probe(task)) {};

		extract("TADSegment",task.returns,segment);
		LOGGERA("TAMMBER-md: OUTPUT SEGMENT "<<segment.info_str())

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

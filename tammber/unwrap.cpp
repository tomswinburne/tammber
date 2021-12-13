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

#include <boost/filesystem.hpp>

int main(int argc, char * argv[]) {
	MPI_Init(&argc, &argv);
	int rank;
	int nranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);

	TaskMapperType mapper;

	std::cout << "TAMMBER-\n";
	// Create empty property tree object
	boost::property_tree::ptree config;
	// Parse the XML into the property tree.
	boost::property_tree::read_xml("./input/ps-config.xml", config, boost::property_tree::xml_parser::no_comments );
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

		GenericTask label,unwrap;

		label.type=mapper.type("TASK_LABEL");
		label.flavor=1;

		unwrap.type=mapper.type("TASK_UNWRAP");
		unwrap.flavor=1;

		PersistentLocalStore minimaStore;
		minimaStore.initialize("./db0/", "min");
		minimaStore.createDatabase(0, false, false);

		bool found_list = boost::filesystem::exists("./Unwrap.list");

		uint64_t canon,label1,label2;
		Transition t;
		std::set<Transition> trans;
		LabelPair initial_labels,final_labels;
		trans.clear();

		if(found_list) {
			std::ifstream infile("./Unwrap.list");
			while(infile>>canon>>label1>>label2) {
				t.first.first=canon;
				t.first.second=label1;
				t.second.first=canon;
				t.second.second=label2;
				trans.insert(t);
			}
			LOGGERA("Found "<<trans.size()<<" unwrapping tasks")

			std::ofstream out;
			out.open("./Unwrap_Results.list", std::ios::out);

			for(auto tran: trans) {

				RawDataVector data;
				SystemType initial,final;

				std::array<double,NDIM*NDIM> matrix;
				std::array<double,NDIM> shift;
				bool valid;

				minimaStore.get(LOCATION_SYSTEM_MIN,tran.first.second,data);
				if(data.size()==0) {
					LOGGERA("COULDN'T FIND "<<tran.first.second)
					continue;
				}
				unpack(data,initial,data.size());
				label.clearInputs(); label.clearOutputs();
				insert("State",label.inputData,initial);
				handle.assign(label);
				while(not handle.probe(label)) {};
				extract("Labels",label.returns,initial_labels);
				if(initial_labels!=tran.first) {
					LOGGERA("COULDN'T MATCH "<<tran.first.first<<" "<<tran.first.second)
					continue;
				}

				minimaStore.get(LOCATION_SYSTEM_MIN,tran.second.second,data);
				if(data.size()==0) {
					LOGGERA("COULDN'T FIND "<<tran.second.second)
					continue;
				}
				unpack(data,final,data.size());
				label.clearInputs(); label.clearOutputs();
				insert("State",label.inputData,final);
				handle.assign(label);
				while(not handle.probe(label)) {};
				extract("Labels",label.returns,final_labels);
				if(final_labels!=tran.second) {
					LOGGERA("COULDN'T MATCH "<<tran.second.first<<" "<<tran.second.second)
					continue;
				}


				LOGGERA("SUBMITTING PATHWAY FOR UNWRAP ")
				unwrap.clearInputs(); unwrap.clearOutputs();
				insert("InitialState",unwrap.inputData,initial);
				insert("FinalState",unwrap.inputData,final);
				insert("InitialLabels",unwrap.inputData,initial_labels);
				insert("FinalLabels",unwrap.inputData,final_labels);
				handle.assign(unwrap);
				while(not handle.probe(unwrap)) {};

				valid = false;
				extract("Matrix",unwrap.returns,matrix);
				extract("Shift",unwrap.returns,shift);
				extract("Valid",unwrap.returns,valid);
				if(valid) {
					out<<initial_labels.first<<" "<<initial_labels.second<<" "<<final_labels.second<<" ";
					for(int j=0;j<NDIM*NDIM;j++) out<<matrix[j]<<" ";
					for(int j=0;j<NDIM;j++) out<<shift[j]<<" ";
					out<<std::endl;
				}
			}
			out.close();
		} else {
			std::cout<<"Unwrap.list not found!"<<std::endl;
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

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


#ifndef __Tammber__Builder__
#define __Tammber__Builder__

#include <new>
#include <stdio.h>
#include <limits>
#include <memory>
#include <atomic>
#include <future>
#include <deque>
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <thread>

#include <mpi.h>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/optional.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/unordered_multiset_of.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/timer/timer.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/random/random_device.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/functional/hash/hash.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/list.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <thread>

#include "CustomTypedefs.hpp"
#include "Types.hpp"
#include "TammberTypes.hpp"
#include "Task.hpp"
#include "Pack.hpp"
#include "Constants.hpp"
#include "DDS.hpp"
#include "TaskManager.hpp"

#include "TammberModel.hpp"

#include "PullWorkProducer.hpp"
#include "Log.hpp"


class TammberModelBuilder : public AbstractPullWorkProducer {
public:
TammberModelBuilder(MPI_Comm comm_,AbstractDDS *sharedStore_, std::set<int> children_, boost::property_tree::ptree &config) :
AbstractPullWorkProducer(comm_,sharedStore_,children_,config){
	LOGGER("TammberModelBuilder : public AbstractPullWorkProducer")
	carryOverTime=0;
	jobcount = 0;
	uint64_t sharedBufferSize=config.get<uint64_t>("Configuration.MarkovModel.SharedCacheSize",1000000000);
	sharedStore->setMaximumBufferSize(sharedBufferSize);
	bool rs=false;
	if(boost::filesystem::exists("./TammberModel.chk")) {
		rs=true;
		std::ifstream ifs("TammberModel.chk");
		boost::archive::text_iarchive ia(ifs);
		// read class state from archive
		ia >> *this;
	}

	//initialize from config data
	//batchSize=config.get<unsigned>("Configuration.MarkovModel.PredictionSize",nWorkers_);
	maxTaskMesgSize=config.get<int>("Configuration.MaximumTaskMesgSize",1000000);
	reportDelay=std::chrono::milliseconds( config.get<unsigned>("Configuration.MarkovModel.ReportDelay",10000) );
	checkpointDelay=std::chrono::milliseconds( config.get<unsigned>("Configuration.MarkovModel.CheckpointDelay",100000) );
	initialConfigurationString=config.get<std::string>("Configuration.InitialConfigurations");
	defaultFlavor=config.get<int>("Configuration.TaskParameters.DefaultFlavor", 0);
	neb_policy = config.get<int>("Configuration.MarkovModel.NEBPolicy",0);

  	deleteVertex = config.get<uint64_t>("Configuration.MarkovModel.DeleteVertex",0);
	
	// initialConfigurations
	boost::split(initialConfigurations,initialConfigurationString,boost::is_any_of(" "));

	// initialize the model
	markovModel.initialize(config,rs);
};


template<class Archive>
void save(Archive & ar, const unsigned int version) const {
	// note, version is always the latest when saving
	long int co=std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count()+carryOverTime;
	ar & jobcount;
	ar & co;
	ar & markovModel;
};

template<class Archive>
void load(Archive & ar, const unsigned int version){
	ar & jobcount;
	ar & carryOverTime;
	ar & markovModel;
};
BOOST_SERIALIZATION_SPLIT_MEMBER()




virtual void initialize(){
	LOGGER("PullMMbuilder::initialize()")
	initialized=false;
	initializeSystems();
	initialized=true;
};


void initializeSystems() {
	LOGGER("PullMMbuilder::initializeSystems()")
	std::string testfn;
	for(auto initialConfiguration : initialConfigurations) {
		boost::trim(initialConfiguration);
		LOGGER("PullMMbuilder::initializeSystems() : loading "<<initialConfiguration)
		if(initialConfiguration.length()==0) continue;
		TaskDescriptor task;
		task.type=mapper.type("TASK_INIT_MIN");
		task.flavor=jobcount++;
		task.optional=false;
		task.imposeOrdering=false;
		task.nInstances=1;
		task.id=jobcount++;
		insert("Filename",task.arguments,initialConfiguration);
		extract("Filename",task.arguments,testfn);
		LOGGERA("REQUESTING "<<testfn)
		taskQueue.insert(task);
		LOGGERA("IN TQ: "<<taskQueue.count())
	}

	std::list<GenericTask> tasks;
	int counts=0, prev_counts=0, total_counts=initialConfigurations.size();
	std::chrono::high_resolution_clock::time_point now,FirstLog,LastLog;
	std::chrono::milliseconds LogDelay,MaxDelay;
	LogDelay = std::chrono::milliseconds(1000); // 1s
	MaxDelay = std::chrono::milliseconds(60000); // 60s
	now = std::chrono::high_resolution_clock::now();
	FirstLog = now;
	LastLog = now;

	while(counts<initialConfigurations.size()) {
		processSend();
		processRecv();
		prev_counts = tasks.size();
		ready.extract(mapper.type("TASK_INIT_MIN"),tasks);
		counts += tasks.size() - prev_counts;
		now = std::chrono::high_resolution_clock::now();
		if(now - LastLog > LogDelay || tasks.size()>prev_counts) {
			LOGGER("PullMMbuilder::initializeSystems() :  "<<counts<<"/"<<total_counts)
			LastLog = now;
		}
		if(now - FirstLog > MaxDelay && tasks.size()>0) break;
	}

	LOGGER("PullMMbuilder::initializeSystems() : received "<<tasks.size()<<" tasks")

	for(auto &tt: tasks) {
		LabelPair labels;
		double energy;
		std::array<double,3> position = {0.,0.,0.};
		std::set<PointShiftSymmetry> self_symmetries;
		PointShiftSymmetry null;
		self_symmetries.insert(null);
		int clusters=1;

		extract("Labels",tt.returns,labels);
		extract("Energy",tt.returns,energy);
		extract("Clusters",tt.returns,clusters);
		extract("Position",tt.returns,position);
		LOGGERA("LABELS: "<<labels.first<<" "<<labels.second<<" E:"<<energy
			<<"eV, Clusters:"<<clusters
			<<" Position:"<<position[0]<<" "<<position[1]<<" "<<position[2])


		#ifdef ISOMORPHIC
		extract("SelfSymmetries",tt.returns,self_symmetries);
		LOGGERA("SelfSymmetries:")
		for(auto ss:self_symmetries) LOGGERA(ss.info_str());
		#endif
		markovModel.add_vertex(labels,energy,clusters,position,self_symmetries);
	}

	LOGGERA(markovModel.info_str()<<"\nINIT DONE")

};


virtual void processCompletedTasks(){
	LOGGER("PullMMbuilder::processCompletedTasks()")
	//process the completed results
	//
	LOGGER("PullMMbuilder::processCompletedTasks PROCESSING "<<ready.size()<<" TASKS")
	for(auto &seg: ready.segments) {
		LOGGER("ADDING SEGMENT: "<<seg.info_str())
		markovModel.add_segment(seg);

	}
	ready.segments.clear();

	for(auto &path: ready.pathways) {
		LOGGER("PullMMbuilder::processCompletedTasks ADDING PATHWAY: "<<path.info_str())
		markovModel.add_pathway(path);
	}
	ready.pathways.clear();
	LOGGER("END PullMMbuilder::processCompletedTasks()")
};


virtual void report_impl(){
	LOGGER("PullMMbuilder::report_impl()")
	//report at given intervals
	{
		Timer t;
		LOGGER(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count()+carryOverTime)
		LOGGERA(markovModel.info_str())
		LOGGERA("PREDICTION WITH 100 WORKERS :")
		std::list<TADjob> jobs;
		markovModel.generateTADs(jobs,100,true);
		completedTasks.report();
		LOGGERA("In TaskQueue: "<<taskQueue.count())
	}
};

virtual void full_print(bool model) {
	LOGGER("PullMMbuilder::full_print()")
	// Save XML file
	if(deleteVertex!=0) {
		LOGGERA("PullMMbuilder::full_print DELETING "<<deleteVertex)
		markovModel.deleteVertex(deleteVertex);
		{
			std::ofstream ofs("./TammberModelNew.chk");
			{
				boost::archive::text_oarchive oa(ofs);
				oa << *this;
			}
		}
	}
	// Print analysis to screen
	LOGGERA(markovModel.info_str(true))
	if(model) markovModel.write_model("MarkovModel.xml");
	std::list<TADjob> jobs;
	LOGGERA("PREDICTION WITH 100 WORKERS :")
	markovModel.generateTADs(jobs,100,true);

};

virtual void checkpoint_impl() {
	LOGGER("PullMMbuilder::checkpoint_impl()")
	//checkpoint at given intervals
	{
		std::ofstream ofs("./TammberModel.chk");
		// save data to archive
		{
			boost::archive::text_oarchive oa(ofs);
			// write class instance to archive
			//this is done by boost serialization, through a call to the save method of PullSplicer (see above).
			oa << *this;
			// archive and stream closed when destructors are called
		}
	}
};


virtual TaskDescriptorBundle generateTasks(int consumerID, int nTasks){

	LOGGER("PullMMbuilder::generateTasks : GENERATING "<<nTasks<<" TASKS, SPLICER TQ: "<<taskQueue.count())

	TaskDescriptorBundle tasks;
	TaskDescriptor task;

	//get the "global" tasks first
	taskQueue.transferTo(tasks,nTasks);
	int batchSize=tasks.count();
	LOGGER("PullMMbuilder::generateTasks : batchSize="<<batchSize<<", nTasks="<<nTasks<<", initialized:"<<initialized)
	if (batchSize>=nTasks || not initialized) return tasks;

	LOGGER("PullMMbuilder::generateTasks : initialized && batchSize<=nTasks")

	std::list<NEBjob> nebs;
	if(neb_policy==0||neb_policy==1) {
		markovModel.generateNEBs(nebs,nTasks-batchSize);
		for(auto neb=nebs.begin(); neb!=nebs.end();) {
		task.type=mapper.type("TASK_NEB");

		task.imposeOrdering=false;
		task.optional=false;

		task.clearInputs();
		task.nInstances=1;
		task.producer=0;
		task.id=jobcount++;
		task.flavor=defaultFlavor;

		NEBPathway pathway;
		pathway.InitialLabels = neb->TargetTransition.first;
		pathway.FinalLabels = neb->TargetTransition.second;
		pathway.pairmap=false;
		pathway.InitialSymmetries = neb->InitialSymmetries;
		pathway.FinalSymmetries = neb->FinalSymmetries;

		LOGGER("PullMMbuilder::generateTasks : SUBMITTING PATHWAY FOR NEB "<<pathway.submit_info_str())
		if(neb->ExistingPairs.size()>0) {
			LOGGER("ExistingPairs: \n")
			for (auto epl: neb->ExistingPairs)
				LOGGER(neb->TargetTransition.first.first<<","<<epl.first<<" -> "<<neb->TargetTransition.second.first<<","<<epl.second<<"\n")
		}

		insert("Initial",neb->TargetTransition.first.second,LOCATION_SYSTEM_MIN,true,NECESSITY::REQUIRED,task.inputData);
		insert("Final",neb->TargetTransition.second.second,LOCATION_SYSTEM_MIN,true,NECESSITY::REQUIRED,task.inputData);

		for(auto epl: neb->ExistingPairs) {
			insert("ExistingPairs",epl.first,LOCATION_SYSTEM_MIN,true,NECESSITY::REQUIRED,task.inputData);
			insert("ExistingPairs",epl.second,LOCATION_SYSTEM_MIN,true,NECESSITY::REQUIRED,task.inputData);
			pathway.compared_transitions.push_back(epl);
		}

		insert("NEBPathway",task.arguments,pathway);

		tasks.insert(task);
		neb = nebs.erase(neb);
		batchSize++;
		if (batchSize>=nTasks) return tasks;
	}
	}
	if(neb_policy==1) return tasks;

	std::list<TADjob> tads;
	
	#ifdef VERBOSE
	markovModel.generateTADs(tads,nTasks-batchSize,true);
	#else
	markovModel.generateTADs(tads,nTasks-batchSize,false);
	#endif
	
	task.type=mapper.type("TASK_SEGMENT");
	task.imposeOrdering=true;
	task.optional=true;

	for(auto tad=tads.begin(); tad!=tads.end();) {
		task.clearInputs();
		task.nInstances=tad->nInstances;
		task.producer=0;
		task.id=jobcount++;
		task.flavor=defaultFlavor;

		TADSegment segment;
		segment.InitialLabels = tad->InitialLabels;
		segment.temperature = tad->temperature;
		for(auto bl: tad->BasinLabels) segment.BasinLabels.insert(bl);

		LOGGER("PullMMbuilder::generateTasks : SUBMITTING SEGMENT "<<segment.submit_info_str())
		bool ProductionRun=true; // not debug run
		insert("ProductionRun",task.arguments,ProductionRun);

		insert("TADSegment",task.arguments,segment);
		insert("Minimum",tad->InitialLabels.second,LOCATION_SYSTEM_MIN,true,NECESSITY::REQUIRED,task.inputData);
		std::string qsdstr = "QSD"+std::to_string(int(tad->temperature));
		insert(qsdstr,tad->InitialLabels.second,task.flavor,false,NECESSITY::OPTIONAL,task.inputData);

		LOGGER("PullMMbuilder::generateTasks : ADDING "<<task.nInstances<<" OF TASK "<<mapper.type(task.type)<<" "<<task.type)

		tasks.insert(task);
		tad = tads.erase(tad);
		batchSize += task.nInstances;
		if (batchSize>=nTasks) return tasks;
	}
	LOGGER("PullMMbuilder::generateTasks : SPLICER TQ: "<<taskQueue.count())

	return tasks;

};


protected:

TaskDescriptorBundle taskQueue;

TammberModel markovModel;
std::string initialConfigurationString;
std::vector<std::string> initialConfigurations;
unsigned long carryOverTime;
unsigned long jobcount;
unsigned batchSize;
int defaultFlavor,neb_policy;
Label deleteVertex;
std::map< std::pair<int,int>, std::map<std::string,std::string> > taskParameters;
//std::ofstream outTime;
bool initialized;


};

#endif

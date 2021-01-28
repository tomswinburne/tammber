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



#ifndef MDEngine_h
#define MDEngine_h

#include "AbstractEngine.hpp"
#include "SystemModifier.hpp"
#include "Task.hpp"
#include "Log.hpp"
#include "Graph.hpp"
#include "TransitionFilter.hpp"
#include <map>
#include <string>
#include <boost/optional.hpp>
#include <boost/functional/hash/hash.hpp>

/*
	MDTaskMapper class for MDEngine inherited from AbstractTaskMapper
*/

class MDTaskMapper : public AbstractTaskMapper {
public:
MDTaskMapper() : AbstractTaskMapper() {
	// worried about clash with data....
	AbstractTaskMapper::insert("TASK_MD");
	AbstractTaskMapper::insert("TASK_MIN");
	AbstractTaskMapper::insert("TASK_SEGMENT");
	AbstractTaskMapper::insert("TASK_FORCES");
	AbstractTaskMapper::insert("TASK_CENTRO");
	AbstractTaskMapper::insert("TASK_INIT_FROM_FILE");
	AbstractTaskMapper::insert("TASK_WRITE_TO_FILE");
	AbstractTaskMapper::insert("TASK_LABEL");
	AbstractTaskMapper::insert("TASK_REMAP");
	AbstractTaskMapper::insert("TASK_INIT_VELOCITIES");
	AbstractTaskMapper::insert("TASK_INIT_MIN");
	AbstractTaskMapper::insert("TASK_MODIFY");
	AbstractTaskMapper::insert("TASK_FILTER_TRANSITION");
};
};

/*
	MDEngine class inherited from AbstractEngine
*/

template <class System, class EngineTaskMapper>
class MDEngine : public AbstractEngine<EngineTaskMapper> {

public:
typedef AbstractEngine<EngineTaskMapper> BaseEngine;

MDEngine(boost::property_tree::ptree &config, MPI_Comm localComm_, int seed_) : BaseEngine(config,localComm_,seed_)  {
	localComm=localComm_;
	seed=seed_;

	//move this to the tasks
	std::string labelerType=config.get<std::string>("Configuration.StateLabeler.Type", "");
	boost::trim(labelerType);
	labeler=labelerFactory.at(labelerType)();

	labeler->initialize(config);

	std::string modifierType=config.get<std::string>("Configuration.SystemModifier.Type", "");
	boost::trim(modifierType);
	modifier=modifierFactory.at(modifierType)();
	modifier->initialize(config);

	defaultFlavor=config.get<int>("Configuration.TaskParameters.DefaultFlavor", 0);

	BOOST_FOREACH(boost::property_tree::ptree::value_type &v, config.get_child("Configuration.TaskParameters")) {
		boost::optional<std::string> otype= v.second.get_optional<std::string>("Task");
		if(otype) {
			std::string stype=*otype;
			boost::trim(stype);
			int type=BaseEngine::mapper.type(stype);
			int flavor=v.second.get<int>("Flavor");
			BOOST_FOREACH(boost::property_tree::ptree::value_type &vv, v.second.get_child("")) {
				std::string key=vv.first;
				std::string data=vv.second.data();
				boost::trim(key);
				boost::trim(data);
				taskParameters[std::make_pair(type,flavor)][key]=data;
			}
		}
	}

	// Assign implementations to dispatchor
	// Any overwritten functions must be reassigned here as well
	BaseEngine::impls["TASK_MD"] = MDEngine::md_impl; // insert("TASK_MD");
	BaseEngine::impls["TASK_MIN"] = MDEngine::min_impl; // insert("TASK_MIN");
	BaseEngine::impls["TASK_SEGMENT"] = MDEngine::segment_impl; // insert("TASK_SEGMENT");
	BaseEngine::impls["TASK_FORCES"] = MDEngine::forces_impl; // insert("TASK_FORCES");VE");
	BaseEngine::impls["TASK_CENTRO"] = MDEngine::centro_impl; // insert("TASK_CENTRO");VE");
	BaseEngine::impls["TASK_INIT_FROM_FILE"] = MDEngine::file_init_impl; // insert("TASK_INIT_FROM_FILE");
	BaseEngine::impls["TASK_WRITE_TO_FILE"] = MDEngine::file_write_impl; // insert("TASK_WRITE_TO_FILE");
	BaseEngine::impls["TASK_LABEL"] = MDEngine::label_impl; // insert("TASK_LABEL");
	BaseEngine::impls["TASK_REMAP"] = MDEngine::remap_impl; // insert("TASK_REMAP");
	BaseEngine::impls["TASK_INIT_VELOCITIES"] = MDEngine::init_velocities_impl; // insert("TASK_INIT_VELOCITIES");
	BaseEngine::impls["TASK_INIT_MIN"] = MDEngine::init_min_impl; // insert("TASK_INIT_MIN");
	BaseEngine::impls["TASK_MODIFY"] = MDEngine::modify_impl; // insert("TASK_MODIFY");
	BaseEngine::impls["TASK_FILTER_TRANSITION"] = MDEngine::filter_transition_impl; // insert("TASK_FILTER_TRANSITION");
};

int defaultFlavor;
//std::map< std::pair<int,int>, std::map<std::string,std::string> > taskParameters;
std::unordered_map< std::pair<int,int>, std::unordered_map<std::string,std::string>, boost::hash< std::pair<int,int> > > taskParameters;

protected:
std::shared_ptr<AbstractStateLabeler> labeler;


private:

int seed;
MPI_Comm localComm;
std::shared_ptr<AbstractSystemModifier> modifier;

std::function<void(GenericTask&)> md_impl = [this](GenericTask &task) {
	/* to be overwritten (in impls as well) */
};
std::function<void(GenericTask&)> min_impl = [this](GenericTask &task) {
	/* to be overwritten (in impls as well) */
};
std::function<void(GenericTask&)> init_velocities_impl = [this](GenericTask &task) {
	/* to be overwritten (in impls as well) */
};
std::function<void(GenericTask&)> forces_impl = [this](GenericTask &task) {
	/* to be overwritten (in impls as well) */
};
std::function<void(GenericTask&)> centro_impl = [this](GenericTask &task) {
	/* to be overwritten (in impls as well) */
};

std::function<void(GenericTask&)> file_init_impl = [this](GenericTask &task) {
	/* to be overwritten (in impls as well) */
};
std::function<void(GenericTask&)> file_write_impl = [this](GenericTask &task) {
	/* to be overwritten (in impls as well) */
};

/**
 * This should expect:
 * arguments: TrajectoryIndex
 * arguments:InitialLabel
 * inputData: InitialState, the configurations to modify
 *
 * This should return:
 * returns: TrajectoryIndex
 * returns: InitialLabel
 * returns: FinalLabel
 * outputData: FinalState
 */
std::function<void(GenericTask&)> modify_impl = [this](GenericTask &task){
	Label initialLabel,finalLabel;
	int trajectory;

	extract("InitialLabel",task.arguments,initialLabel);
	extract("TrajectoryIndex",task.arguments,trajectory);
	insert("InitialLabel",task.returns,initialLabel);
	insert("TrajectoryIndex",task.returns,trajectory);

	System s;
	extract("InitialState",task.inputData,s);

	std::unordered_map<std::string,std::string> parameters=extractParameters(task.type,task.flavor,defaultFlavor,taskParameters);
	std::string modifierType=parameters["Type"];
	boost::trim(modifierType);

	// DO SOMETHING
	finalLabel=initialLabel;
	insert("FinalLabel", task.returns, finalLabel);

};

std::function<void(GenericTask&)> filter_transition_impl = [this](GenericTask &task){
	std::unordered_map<std::string,std::string>
		parameters = extractParameters(task.type,task.flavor,defaultFlavor,taskParameters);
	task.clearOutputs();

	System reference,state;
	extract("ReferenceState",task.inputData,reference);
	extract("State",task.inputData,state);

	std::string filterType = parameters["Type"];
	if(task.arguments.find("Type")!=task.arguments.end())
		extract("Type",task.arguments,filterType);
	boost::trim(filterType);

	std::shared_ptr<AbstractTransitionFilter>
		filter = transitionFilterFactory.at(filterType)();

	filter->initialize(parameters);

	bool valid = filter->isValid(reference, state, parameters);
	LOGGER("MDEngine::FilterTransition: "<<valid<<" dX: "<<dX<<" dXmax: "<<dXmax)

	insert("Valid", task.returns, valid);
};

virtual void min_label_remap(GenericTask &task){
	task.clearOutputs();

	Label label;
	GenericTask t;
	t.arguments=task.arguments;
	System s;
	t.flavor=task.flavor;

	//minimize
	extract("State",task.inputData,s);
	t.inputData.clear();
	insert("State",t.inputData,s);
	t.type = BaseEngine::mapper.type("TASK_MIN");
	BaseEngine::process(t);

	// label
	extract("State",t.outputData,s);
	t.inputData.clear();
	insert("State",t.inputData,s);
	t.type=BaseEngine::mapper.type("TASK_LABEL");
	BaseEngine::process(t);
	extract("Label",t.returns, label);


	// remap
	t.inputData.clear();
	t.returns.clear();
	insert("ReferenceState",t.inputData,s);
	insert("State",t.inputData,s);
	t.type=BaseEngine::mapper.type("TASK_REMAP");
	BaseEngine::process(t);

	bool remapped;
	extract("Remapped",t.returns, remapped);
	insert("Remapped",task.returns, remapped);

	if(remapped) {
		extract("State",t.outputData,s);
		extract("Label",t.returns, label);
	}
	insert("Label", task.returns, label);
	insert("State",label, LOCATION_SYSTEM_MIN, true, task.outputData, s);
};

/**
 * read a configuration from a file, minimize, and label it.
 *
 * This expects:
 * arguments: Filename
 *
 * This returns
 * returns: Label
 * outputData State
 */
std::function<void(GenericTask&)> init_min_impl = [this](GenericTask &task) {

	Label label;
	GenericTask t;
	t.arguments=task.arguments;
	System s;
	t.flavor=task.flavor;

	//read the file
	t.type=BaseEngine::mapper.type("TASK_INIT_FROM_FILE");
	BaseEngine::process(t);
	extract("State",t.outputData,s);

	//
	t.inputData.clear();
	insert("State",t.inputData,s);
	min_label_remap(t);

	extract("State",t.outputData,s);
	extract("Label",t.returns, label);
	insert("Label", task.returns, label);
	insert("State",label, LOCATION_SYSTEM_MIN, true, task.outputData, s);
};

/**
 * Implement segment generation in terms of other more basic tasks. Can be overridden in derived classes if the engine can generate segments internally.
 *
 * This expects:
 * inputData: Minimum
 * inputData: QSD (optional)
 *
 * This returns
 * returns: Trajectory
 * outputData: FinalMin
 * outputData:
 */
std::function<void(GenericTask&)> segment_impl = [this](GenericTask &task){


	std::unordered_map<std::string,std::string> parameters=extractParameters(task.type,task.flavor,defaultFlavor,taskParameters);

	//read parameters from the task
	double preCorrelationTime=boost::lexical_cast<double>(parameters["PreCorrelationTime"]);
	double postCorrelationTime=boost::lexical_cast<double>(parameters["PostCorrelationTime"]);
	double minimumSegmentLength=boost::lexical_cast<double>(parameters["MinimumSegmentLength"]);
	double blockTime=boost::lexical_cast<double>(parameters["BlockTime"]);
	int nDephasingTrialsMax=boost::lexical_cast<int>(parameters["MaximumDephasingTrials"]);
	bool reportIntermediates=boost::lexical_cast<bool>(parameters["ReportIntermediates"]);

	int segmentFlavor=(taskParameters.count(std::make_pair(task.type,task.flavor))>0 ? task.flavor : defaultFlavor);

	std::set<int> contributeSegmentsTo;
	std::set<int> contributeStatisticsTo;
	if(parameters.count("ContributeSegmentsTo")>0) {
		std::istringstream is( parameters["ContributeSegmentsTo"]);
		contributeSegmentsTo=std::set<int>{ std::istream_iterator<int>( is ), std::istream_iterator<int>() };
	}
	else{
		contributeSegmentsTo.insert(segmentFlavor);
	}

	/*
	      if(parameters.count("ContributeStatisticsTo")>0) {
	              std::istringstream is( parameters["ContributeStatisticsTo"]);
	              contributeStatisticsTo=std::set<int>{ std::istream_iterator<int>( is ), std::istream_iterator<int>() };
	      }
	      else{
	              contributeStatisticsTo.insert(segmentFlavor);
	      }
	 */


	int maximumSegmentLength;
	if(parameters.count("MaximumSegmentLength")>0 ) {
		maximumSegmentLength=boost::lexical_cast<int>(parameters["MaximumSegmentLength"]);
	}
	else{
		maximumSegmentLength=25*minimumSegmentLength;
	}

	//create tasks
	GenericTask md;
	md.type=BaseEngine::mapper.type("TASK_MD");
	md.flavor=task.flavor;

	GenericTask min;
	min.type=BaseEngine::mapper.type("TASK_MIN");
	min.flavor=task.flavor;

	GenericTask label;
	label.type=BaseEngine::mapper.type("TASK_LABEL");;
	label.flavor=task.flavor;

	GenericTask remap;
	remap.type=BaseEngine::mapper.type("TASK_REMAP");;
	remap.flavor=task.flavor;

	GenericTask initVelocities;
	initVelocities.type=BaseEngine::mapper.type("TASK_INIT_VELOCITIES");;
	initVelocities.flavor=task.flavor;

	GenericTask filterTransitions;
	filterTransitions.type=BaseEngine::mapper.type("TASK_FILTER_TRANSITION");;
	filterTransitions.flavor=task.flavor;



	//extract the systems we were provided
	System minimum;
	System reference;
	System qsd;
	System initial;
	System current;
	System currentMin;

	bool gotMin=extract("Minimum",task.inputData,minimum);
	bool gotQsd=extract("QSD",task.inputData,qsd);
	//set the initial state
	if(gotQsd) {
		initial=qsd;
	}
	else{
		initial=minimum;
	}
	currentMin=minimum;

	Label initialLabel,currentLabel;
	insert("State",label.inputData,minimum);
	BaseEngine::process(label);
	extract("Label",label.returns,initialLabel);

	//set labels
	currentLabel=initialLabel;


	bool dephased=gotQsd;
	//sample thermal velocities
	if(not dephased) {
		insert("State",initVelocities.inputData,initial);
		BaseEngine::process(initVelocities);
		extract("State",initVelocities.outputData,initial);
	}

	current=initial;

	int nDephasingTrials=0;
	int nOverheadBlocks=0;

	while(not dephased) {
		//dephasing loop
		current=initial;



		double elapsedTime=0;
		while( elapsedTime  < preCorrelationTime*0.999999999    ) {
			//run md
			md.inputData.clear();
			insert("State",md.inputData,current);
			BaseEngine::process(md);
			elapsedTime+=blockTime;
			nOverheadBlocks++;
			extract("State",md.outputData,current);

			//minimize
			min.inputData.clear();
			insert("State",min.inputData,current);
			BaseEngine::process(min);
			extract("State",min.outputData,currentMin);
			//currentMin=current;

			//label
			label.inputData.clear();
			insert("State",label.inputData,currentMin);
			BaseEngine::process(label);
			extract("Label",label.returns,currentLabel);

			if( currentLabel == initialLabel ) {
				dephased=true;
			}
			else{
				dephased=false;
				break;
			}
		}

		nDephasingTrials++;
		if(nDephasingTrials>=nDephasingTrialsMax) {
			LOGGERA("DEPHASING FAILED")
			if(reportIntermediates) {
				insert("State",currentLabel,LOCATION_SYSTEM_MIN,true,task.outputData,currentMin);
			}
			break;
		}
	}


	bool segmentIsSpliceable=false;
	double lastTransitionTime=-postCorrelationTime;
	int nBlocks=0;
	double elapsedTime=0;
	bool segmentIsValid=true;

	Trajectory trajectory;
	Visit v;
	v.label=initialLabel;
	v.duration=0;
	trajectory.appendVisit(v);
	v.label=currentLabel;
	v.duration=0;
	trajectory.appendVisit(v);
	trajectory.overhead()=nOverheadBlocks;

	while( !segmentIsSpliceable or elapsedTime  < minimumSegmentLength*0.999999999 ) {
		reference=currentMin;

		//take a block of MD
		md.inputData.clear();
		insert("State",md.inputData,current);
		BaseEngine::process(md);
		extract("State",md.outputData,current);
		nBlocks++;
		elapsedTime+=blockTime;

		//append to the trajectory
		if(reportIntermediates) {
			v.label=currentLabel;
		}
		v.duration=1;
		trajectory.appendVisit(v);

		min.inputData.clear();
		insert("State",min.inputData,current);
		BaseEngine::process(min);
		extract("State",min.outputData,currentMin);

		//hash current state
		Label previousLabel=currentLabel;
		label.inputData.clear();
		insert("State",label.inputData,currentMin);
		BaseEngine::process(label);
		extract("Label",label.returns,currentLabel);



		if( currentLabel!=previousLabel ) {

			filterTransitions.clearInputs();
			filterTransitions.clearOutputs();
			insert("State",filterTransitions.inputData,currentMin);
			insert("ReferenceState",filterTransitions.inputData,reference);
			BaseEngine::process(filterTransitions);
			extract("Valid",filterTransitions.returns,segmentIsValid);

			lastTransitionTime=elapsedTime;
			if(reportIntermediates) {
				insert("State",currentLabel,LOCATION_SYSTEM_MIN,true,task.outputData,currentMin);
			}
		}

		//a segment is spliceable if the last transition occurred at least postCorrelationTime in the past
		segmentIsSpliceable=(elapsedTime-lastTransitionTime >= postCorrelationTime*0.999999999);

		if(trajectory.duration()>=maximumSegmentLength) {
			LOGGER("SEGMENT EXCEEDED MAXIMUM LENGTH. BAILING OUT.")
			break;
		}

		if(not segmentIsValid) {
			break;
		}
	}

	if(not segmentIsValid) {
		task.clearOutputs();
		//pack an empty trajectory so that nobody waits for this segment in vain
		trajectory.clear();
		Visit v;
		v.label=initialLabel;
		v.duration=0;
		trajectory.appendVisit(v);
		trajectory.overhead()=nOverheadBlocks+nBlocks;

		SegmentDatabase db;
		std::map<int,TransitionStatistics> stats;
		TransitionStatistics ts;

		for(auto it=contributeSegmentsTo.begin(); it!=contributeSegmentsTo.end(); it++) {
			db.add(*it,trajectory);
			stats[*it]=ts;
		}

		//pack the trajectory
		insert("Statistics",task.returns,stats);
		insert("Trajectory",task.returns,db);
		task.clearInputs();
		bool invalid=true;
		insert("InvalidTransition",task.returns,invalid);
		return;
	}

	if(not reportIntermediates) {
		v.label=currentLabel;
		v.duration=0;
		trajectory.appendVisit(v);
	}

	//leave a marker that we had to end the segment prematurely
	if(!segmentIsSpliceable) {
		v.label=666;
		v.duration=0;
		trajectory.appendVisit(v);
	}

	System remappedMin,remappedQSD;

	Label remappedLabel;
	//remap minimum and qsd to canonical representation
	remap.inputData.clear();
	insert("ReferenceState",remap.inputData,currentMin);
	insert("State",remap.inputData,currentMin);
	BaseEngine::process(remap);
	bool remapped;
	extract("Remapped",remap.returns,remapped);
	if(remapped) {
		extract("Label",remap.returns,remappedLabel);
		extract("State",remap.outputData,remappedMin);

		remap.inputData.clear();
		insert("ReferenceState",remap.inputData,currentMin);
		insert("State",remap.inputData,current);
		BaseEngine::process(remap);
		extract("State",remap.outputData,remappedQSD);
	}
	else{
		remappedLabel=currentLabel;
		remappedMin=currentMin;
		remappedQSD=current;
	}



	insert("FinalMinimum",remappedLabel,LOCATION_SYSTEM_MIN,true,task.outputData,remappedMin);

	//placeholder to signal a switch to a canonical representative
	if(currentLabel!=remappedLabel) {
		v.label=0;
		v.duration=0;
		trajectory.appendVisit(v);
	}
	v.label=remappedLabel;
	v.duration=0;
	trajectory.appendVisit(v);

	SegmentDatabase db;
	std::map<int,TransitionStatistics> stats;
	TransitionStatistics ts;
	ts.update(trajectory.front().label,trajectory.back().label);

	for(auto it=contributeSegmentsTo.begin(); it!=contributeSegmentsTo.end(); it++) {
		insert("FinalQSD",remappedLabel,*it,false,task.outputData,remappedQSD);
		db.add(*it,trajectory);
		stats[*it]=ts;
	}
	//pack the trajectory
	insert("Statistics",task.returns,stats);
	insert("Trajectory",task.returns,db);
	task.clearInputs();
};

std::function<void(GenericTask&)> label_impl = [this](GenericTask &task){
	bool canonical=false;
	System s;
	extract("State",task.inputData,s);
	Label label=labeler->hash(s,canonical);
	insert("Label",task.returns,label);
};

//joint remapping of a number of states. This assumes that all systems in the task correspond to the same state. The first state is used as a reference
std::function<void(GenericTask&)> remap_impl = [this](GenericTask &task){

	std::unordered_map<std::string,std::string>
		parameters=extractParameters(task.type,task.flavor,defaultFlavor,taskParameters);

	bool canonical=boost::lexical_cast<bool>(parameters["Canonical"]);

	bool remapped=false;

	if(canonical) {
		System reference;
		System s;
		bool gotRef=extract("ReferenceState",task.inputData,reference);
		bool gotState=extract("State",task.inputData,s);
		if (gotRef and gotState) {
			std::map<int,int> canonicalMap;
			Label lb;
			labeler->canonicalMap(reference,canonicalMap,lb);
			s.remap(canonicalMap);
			remapped=true;
			insert("Remapped",task.returns,remapped);
			insert("Label",task.returns,lb);
			auto it=task.inputData.find("State");
			insert("State",lb,it->second.location, it->second.shared,task.outputData,s);
		}
	} else insert("Remapped",task.returns,remapped);
};

};

#endif

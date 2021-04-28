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



#ifndef TADEngine_h
#define TADEngine_h

#include "AbstractEngine.hpp"
#include "SystemModifier.hpp"
#include "TammberTypes.hpp"
#include "Task.hpp"
#include "Log.hpp"
#include "Graph.hpp"
#include "MDEngine.hpp"
#include "TransitionFilter.hpp"
#include <map>
#include <string>
#include <boost/optional.hpp>
#include <boost/functional/hash/hash.hpp>

#include <random>


#include <Eigen/Dense> // For Supercell and Hessian
#include <Eigen/Eigenvalues> // Hessian



/*
	TADTaskMapper class for TADEngine inherited from MDTaskMapper
*/

class TADTaskMapper : public MDTaskMapper {
public:
TADTaskMapper() : MDTaskMapper() {
	MDTaskMapper::AbstractTaskMapper::insert("TASK_NEB");
};

};

/*
	TADEngine class inherited from MDEngine
*/

template <class System, class EngineTaskMapper>
class TADEngine : public MDEngine<System,EngineTaskMapper> {

public:
//typedef AbstractEngine<EngineTaskMapper> BaseEngine; // as per MDEngine
typedef MDEngine<System,EngineTaskMapper> BaseMDEngine; // req. for taskParameters, defaultFlavor, labeler
typedef typename BaseMDEngine::BaseEngine BaseEngine; // req. for mapper, impls (could use this-> ....)

TADEngine(boost::property_tree::ptree &config, MPI_Comm localComm_, int seed_) : BaseMDEngine(config, localComm_, seed_)  {
	BaseEngine::impls["TASK_SEGMENT"] = TADEngine::segment_impl; // overwriting
	BaseEngine::impls["TASK_LABEL"] = TADEngine::label_impl; // overwriting
	BaseEngine::impls["TASK_REMAP"] = TADEngine::remap_impl; // overwriting
	BaseEngine::impls["TASK_INIT_MIN"] = TADEngine::init_min_impl; // overwriting
	BaseEngine::impls["TASK_NEB"] = TADEngine::neb_impl; // new task
};

/**
 * Implement segment generation in terms of other more basic tasks. Can be overridden in derived classes if the engine can generate segments internally.
 *
 * This expects:
 * inputData: Minimum
 * inputData: QSD (optional)
 * argument: TADSegment

 */
std::function<void(GenericTask&)> segment_impl = [this](GenericTask &task) {
	std::unordered_map<std::string,std::string> parameters =
		extractParameters(task.type,task.flavor,BaseMDEngine::defaultFlavor,BaseMDEngine::taskParameters);

	//read parameters from the tasks
	double preCorrelationTime=safe_extractor<double>(parameters,"PreCorrelationTime",1.0);
	double minimumSegmentLength=safe_extractor<double>(parameters,"MinimumSegmentLength",1.0);
	double BlockTime=safe_extractor<double>(parameters,"BlockTime",1.0);
	double annealingTime=safe_extractor<double>(parameters,"AnnealingTime",0.25);
	double annealingTemperature=safe_extractor<double>(parameters,"AnnealingTemperature",0.25);

	int nDephasingTrialsMax=safe_extractor<int>(parameters,"MaximumDephasingTrials",3);
	bool reportIntermediates=safe_extractor<bool>(parameters,"ReportIntermediates",false);
	int segmentFlavor=(BaseMDEngine::taskParameters.count(std::make_pair(task.type,task.flavor))>0 ? task.flavor : BaseMDEngine::defaultFlavor);
	int maximumSegmentLength = safe_extractor<int>(parameters,"MaximumSegmentLength",25*minimumSegmentLength);

	//create tasks
	GenericTask md,min,label,carve,initVelocities,filter,write;
	write.type = BaseEngine::mapper.type("TASK_WRITE_TO_FILE");
	write.flavor = task.flavor;

	md.type=BaseEngine::mapper.type("TASK_MD");
	md.flavor=task.flavor;

	min.type=BaseEngine::mapper.type("TASK_MIN");
	min.flavor=task.flavor;

	label.type=BaseEngine::mapper.type("TASK_LABEL");
	label.flavor=task.flavor;

	carve.type=BaseEngine::mapper.type("TASK_CARVE");
	carve.flavor=task.flavor;

	initVelocities.type=BaseEngine::mapper.type("TASK_INIT_VELOCITIES");
	initVelocities.flavor=task.flavor;

	filter.type=BaseEngine::mapper.type("TASK_FILTER_TRANSITION");
	filter.flavor=task.flavor;

	// in input and output data here
	TADSegment segment;
	extract("TADSegment",task.arguments,segment);
	bool ProductionRun=true;
	extract("ProductionRun",task.arguments,ProductionRun);

	//extract the systems we were provided
	System minimum, reference, qsd, initial, current, currentMin, annealingMin;

	bool gotMin = extract("Minimum",task.inputData,minimum);

	// to be completed
	std::string qsdstr = "QSD"+std::to_string(int(segment.temperature));
	bool gotQsd = extract(qsdstr,task.inputData,qsd);

	if(!gotMin) {
		if(gotQsd) {
			insert("State",min.inputData,qsd);
			BaseEngine::process(min);
			gotMin = extract("State",min.outputData,minimum);
		} else {
			if(BaseEngine::local_rank==0) LOGGER("NO QSD OR MINIMUM! EXITING")
			return ;
		}
	}

	// segment temperature
	double temperature = segment.temperature;
	double inittemperature = 2.0 * temperature;

	//set the initial state
	if(gotQsd) initial = qsd;
	else initial = minimum;

	// minimize and find labels
	LabelPair InitialLabels,CurrentLabels,previousLabels;

	// minimize
	min.clearInputs(); min.clearOutputs();
	insert("State",min.inputData,minimum);
	BaseEngine::process(min);
	extract("State",min.outputData,minimum);

	// label
	label.clearInputs(); label.clearOutputs();
	insert("State",label.inputData,minimum);
	BaseEngine::process(label);
	extract("Labels",label.returns,InitialLabels);
	segment.InitialLabels = InitialLabels;

	// carve
	int clusters=1;
	std::array<double,3> position;
	carve.clearInputs(); carve.clearOutputs();
	insert("State",carve.inputData,minimum);
	BaseEngine::process(carve);

	extract("Clusters",carve.returns,clusters);
	extract("Position",carve.returns,position);

	segment.initialClusters = clusters;
	for(int jj=0;jj<3;jj++) segment.initialPosition[jj] = position[jj];


	// ensure basin is populated
	segment.BasinLabels[InitialLabels]=minimum.getEnergy();
	std::map<LabelPair,double> newBasinLabels;

	// transition check
	bool BasinTransition = false;
	bool NewBasinTransition = false;
	bool Transition = false;
	bool Annealed = false;

	// set labels
	CurrentLabels=InitialLabels;

	// DEPHASING ROUTINE
	// Sample velocities if required, then ensure system stays in superbasin for at least 1ps...
	//  current test: do not put back in but must have repeat visits within nDephasingTrials== new tau_c
	// DO REMAPPING AFTER BASIN TRANSITIONS?? no; all done in MarkovModelBuilder
	int nDephasingTrials = 0;
	int nOverheadBlocks = 0;

	double msd_thresh = 0.0;
	double elapsedTime = 0.0;
	segment.dephased = gotQsd; // so no need to do dephasing if we have QSD...
	segment.overhead = 0;
	// sample thermal velocities
	if(not segment.dephased) {
		initVelocities.clearInputs(); initVelocities.clearOutputs();
		insert("InitTemperature",initVelocities.arguments,inittemperature);
		insert("State",initVelocities.inputData,initial);
		BaseEngine::process(initVelocities);
		extract("State",initVelocities.outputData,initial);
	}

	current = initial;
	reference = minimum;

	while(not segment.dephased) {

		elapsedTime = 0.0;
		while( elapsedTime  < preCorrelationTime*0.999999999 ) {
			// one block of MD
			md.clearInputs(); md.clearOutputs();
			insert("BlockTime",md.arguments,BlockTime);
			insert("Temperature",md.arguments,temperature);
			insert("State",md.inputData,current);
			BaseEngine::process(md);
			extract("State",md.outputData,current);
			elapsedTime += BlockTime;
			segment.overhead++;

			// minimize
			min.clearInputs(); min.clearOutputs();
			insert("State",min.inputData,current);
			BaseEngine::process(min);
			extract("State",min.outputData,currentMin);


			// label
			label.clearInputs(); label.clearOutputs();
			insert("State",label.inputData,currentMin);
			BaseEngine::process(label);
			extract("Labels",label.returns,CurrentLabels);
			if(ProductionRun) {
				if(BaseEngine::local_rank==0) LOGGER("Dephase: REFERENCE CURRENTMIN MSD_2: "<<reference.msd(currentMin,false))
				if(BaseEngine::local_rank==0) LOGGER("Dephase: REFERENCE CURRENTMIN MSD_INF: "<<reference.msd(currentMin,true))
				if(BaseEngine::local_rank==0) LOGGER("Dephase: CURRENT LABEL: "<<CurrentLabels.first<<" , "<<CurrentLabels.second)
			} else {
				if(BaseEngine::local_rank==0) LOGGERA("Dephase: REFERENCE CURRENTMIN MSD_2: "<<reference.msd(currentMin,false))
				if(BaseEngine::local_rank==0) LOGGERA("Dephase: REFERENCE CURRENTMIN MSD_INF: "<<reference.msd(currentMin,true))
				if(BaseEngine::local_rank==0) LOGGERA("Dephase: CURRENT LABEL: "<<CurrentLabels.first<<" , "<<CurrentLabels.second)
			}


			// Current "Basin" implementation: anything below MSD thresh from reference

			// check for transition out of "superbasin"
			BasinTransition = bool(segment.BasinLabels.find(CurrentLabels)!=segment.BasinLabels.end());
			NewBasinTransition = bool(newBasinLabels.find(CurrentLabels)!=newBasinLabels.end());

			// See TASK_FILTER_TRANSITION for how "Valid" is determined
			filter.clearInputs(); filter.clearOutputs();
			insert("ReferenceState",filter.inputData,reference);
			insert("State",filter.inputData,currentMin);
			BaseEngine::process(filter);
			extract("Valid",filter.returns,Transition);

			// If a new state is seen it must be revisited for the basin dephased
			// currently we just log the basin transitions, we do not use them
			// segment.dephased = !Transition or BasinTransition or NewBasinTransition;
			segment.dephased = !Transition;
			if(!ProductionRun) if(BaseEngine::local_rank==0) LOGGERA("Dephased: "<<segment.dephased<<" "<<Transition<<" "<<BasinTransition<<" "<<NewBasinTransition);
			if(BaseEngine::local_rank==0) LOGGER("Dephased: "<<segment.dephased<<" "<<Transition<<" "<<BasinTransition<<" "<<NewBasinTransition);

			if(Transition) { // i.e. a transition. We
				carve.clearInputs(); carve.clearOutputs();
				insert("State",carve.inputData,currentMin);
				BaseEngine::process(carve);
				extract("Clusters",carve.returns,clusters);
				//segment.dephased = bool(clusters==segment.initialClusters);
			}

			// i.e. unknown transition
			if(!BasinTransition and !NewBasinTransition and !Transition) {
				newBasinLabels[CurrentLabels] = currentMin.getEnergy();
				if(reportIntermediates)
					insert("State",CurrentLabels.first,CurrentLabels.second,LOCATION_SYSTEM_MIN,true,task.outputData,currentMin);
				if(!ProductionRun) {
					std::string _fn = "states/state-"+std::to_string(CurrentLabels.first)+"-"+std::to_string(CurrentLabels.second)+".dat";
					write.clearInputs(); write.clearOutputs();
		 			insert("Filename",write.arguments,_fn);
		 			insert("State",write.inputData,currentMin);
		 			BaseEngine::process(write);
				}
			}

			if(!Transition) {
				CurrentLabels=InitialLabels;
				currentMin = initial;
			}
			// if the system is in a new basin state is seen, then it gets promoted to 'true' basin state....
		}
		nDephasingTrials++;
		if(nDephasingTrials>=nDephasingTrialsMax) break;

		if(!Transition) break;
	}
	segment.trials = nDephasingTrials;
	// if we have dephased, newBasinLabels have participated in superbasin dephasing so are valid
	if(segment.dephased) for(auto &nbl: newBasinLabels) segment.BasinLabels.insert(nbl);

	if(segment.dephased) {
		if(!ProductionRun) if(BaseEngine::local_rank==0) LOGGERA("DEPHASED WITH "<<newBasinLabels.size()<<" NEW BASIN STATES (NOT CURRENTLY USED)")
		if(BaseEngine::local_rank==0) LOGGER("DEPHASED WITH "<<newBasinLabels.size()<<" NEW BASIN STATES (NOT CURRENTLY USED)")
	} else {
		if(!ProductionRun) if(BaseEngine::local_rank==0) LOGGERA("NOT DEPHASED; TRIED "<<newBasinLabels.size()<<" NEW FAKE BASIN STATES; EXITING")
		if(BaseEngine::local_rank==0) LOGGER("NOT DEPHASED; TRIED "<<newBasinLabels.size()<<" NEW FAKE BASIN STATES; EXITING")
	}

	if(!segment.dephased) {
		insert("TADSegment",task.returns,segment);
		task.clearInputs();
		return ;
	}

	segment.duration = 0;
	segment.elapsedTime = 0.0;

	bool annealing = false;
	double segmentBlockTime=BlockTime;
	double segmentTemperature = segment.temperature;

	reference=currentMin;
	qsd = current;

	previousLabels = InitialLabels;

	while ( segment.duration<=maximumSegmentLength ) {

		if (annealing) {
			segmentBlockTime = annealingTime * BlockTime;
			segmentTemperature = annealingTemperature;
		}
		else {
			segmentBlockTime = BlockTime;
			segmentTemperature = segment.temperature;
		}

		//take a block of MD
		md.clearInputs(); md.clearOutputs();
		insert("BlockTime",md.arguments,segmentBlockTime);
		insert("Temperature",md.arguments,segmentTemperature);
		insert("State",md.inputData,current);
		BaseEngine::process(md);
		extract("State",md.outputData,current);

		if(!annealing) {
			segment.duration += 1;
			segment.elapsedTime += BlockTime;
		}

		previousLabels=CurrentLabels;

		// minimize the current state
		min.clearInputs(); min.clearOutputs();
		insert("State",min.inputData,current);
		BaseEngine::process(min);
		extract("State",min.outputData,currentMin);
		if(!ProductionRun && BaseEngine::local_rank==0) LOGGERA("REFERENCE CURRENTMIN MSD_INF (1st MIN): "<<reference.msd(currentMin,true));
		if(BaseEngine::local_rank==0) LOGGER("REFERENCE CURRENTMIN MSD_INF (1st MIN): "<<reference.msd(currentMin,true));

		// just to avoid future issues..... count this?
		min.clearInputs(); min.clearOutputs();
		insert("State",min.inputData,currentMin);
		BaseEngine::process(min);
		extract("State",min.outputData,currentMin);

		// hash the current minimized state
		label.clearInputs(); label.clearOutputs();
		insert("State",label.inputData,currentMin);
		BaseEngine::process(label);
		extract("Labels",label.returns,CurrentLabels);

		// logs
		msd_thresh = reference.msd(currentMin,true);
		if(!ProductionRun) {
			if(BaseEngine::local_rank==0) LOGGERA("REFERENCE CURRENTMIN MSD_INF (2nd MIN): "<<msd_thresh)
			if(BaseEngine::local_rank==0) LOGGERA("CURRENT LABEL: "<<CurrentLabels.first<<" , "<<CurrentLabels.second)
		}
		if(BaseEngine::local_rank==0) LOGGER("REFERENCE CURRENTMIN MSD_INF (2nd MIN): "<<msd_thresh)
		if(BaseEngine::local_rank==0) LOGGER("CURRENT LABEL: "<<CurrentLabels.first<<" , "<<CurrentLabels.second)

		if(annealing) {

			// check that annealing did not change the minima
			filter.clearInputs(); filter.clearOutputs();
			insert("ReferenceState",filter.inputData,annealingMin);
			insert("State",filter.inputData,currentMin);
			BaseEngine::process(filter);
			extract("Valid",filter.returns,Transition);
			Annealed = !Transition;

			// check that the minima is not the initial minima
			filter.clearInputs(); filter.clearOutputs();
			insert("ReferenceState",filter.inputData,reference);
			insert("State",filter.inputData,currentMin);
			BaseEngine::process(filter);
			extract("Valid",filter.returns,Transition);
			
			if(Annealed and Transition) { // no transition after annealing == exit
				if(!ProductionRun) if(BaseEngine::local_rank==0) LOGGERA("ANNEALED! VALID TRANSITION MADE!");
				if(BaseEngine::local_rank==0) LOGGER("ANNEALED! VALID TRANSITION MADE!");
				insert("FinalMinimum",CurrentLabels.first,CurrentLabels.second,LOCATION_SYSTEM_MIN,true,task.outputData,currentMin);

				// carve
				carve.clearInputs(); carve.clearOutputs();
				insert("State",carve.inputData,currentMin);
				BaseEngine::process(carve);
				extract("Clusters",carve.returns,clusters);
				extract("Position",carve.returns,position);
				segment.initialE = reference.getEnergy();
				segment.finalE = currentMin.getEnergy();
				segment.OutOfBasin = true;
				segment.finalClusters = clusters;
				for(int jj=0;jj<3;jj++) segment.finalPosition[jj] = position[jj];

				segment.transition.first = InitialLabels;
				segment.transition.second = CurrentLabels;
				if(!ProductionRun && BaseEngine::local_rank==0) LOGGERA("INSERTING "<<CurrentLabels.first<<","<<CurrentLabels.second<<" , E="<<currentMin.getEnergy()<<", NClust = "<<clusters<<", Position: "<<position[0]<<" "<<position[1]<<" "<<position[2])
				if(BaseEngine::local_rank==0) LOGGER("INSERTING "<<CurrentLabels.first<<","<<CurrentLabels.second<<" , E="<<currentMin.getEnergy()<<", NClust = "<<clusters<<", Position: "<<position[0]<<" "<<position[1]<<" "<<position[2])

				annealing = false; // not annealing any more
				break;
			} else { // annealing + !transition == continue
				if(!ProductionRun && BaseEngine::local_rank==0) LOGGERA("NO VALID TRANSITION AFTER ANNEALING! CONTINUING")
				if(BaseEngine::local_rank==0) LOGGER("NO VALID TRANSITION AFTER ANNEALING! CONTINUING")
				annealing = false;
			}
		} else {

			// transition check
			//BasinTransition = bool(segment.BasinLabels.find(CurrentLabels)!=segment.BasinLabels.end());
			filter.clearInputs(); filter.clearOutputs();
			insert("ReferenceState",filter.inputData,reference);
			insert("State",filter.inputData,currentMin);
			BaseEngine::process(filter);
			extract("Valid",filter.returns,Transition);


			//if(Transition and (not BasinTransition) ) { // no basin labels yet
			if(Transition) { // not annealing + not no transition == transition has been made
				if(!ProductionRun && BaseEngine::local_rank==0)  LOGGERA("TRANSITION DETECTED! ANNEALING FOR AnnealingTime BLOCKS AT AnnealingTemperature")
				if(BaseEngine::local_rank==0) LOGGER("TRANSITION DETECTED! ANNEALING FOR AnnealingTime BLOCKS AT AnnealingTemperature")
				segment.transition.first = InitialLabels;
				annealing=true;
				annealingMin = currentMin;
				if(!ProductionRun) {
					std::string _fn = "states/state-"+std::to_string(CurrentLabels.first)+"-"+std::to_string(CurrentLabels.second)+".dat";
					write.clearInputs(); write.clearOutputs();
		 			insert("Filename",write.arguments,_fn);
		 			insert("State",write.inputData,currentMin);
		 			BaseEngine::process(write);
				}
			} else { // no transition, continue MD

				// reset labels as the only change was due to small fluctuations
				if(!Transition and (CurrentLabels!=InitialLabels)) {
					if(!ProductionRun && BaseEngine::local_rank==0) LOGGERA("NO TRANSITION FROM MSD CHECK");
					if(BaseEngine::local_rank==0) LOGGER("NO TRANSITION FROM MSD CHECK");
					CurrentLabels=InitialLabels; //
				}

				segment.transition.first = InitialLabels;
				segment.transition.second = CurrentLabels;
			}
		}
		if(segment.duration>=maximumSegmentLength and !annealing) {
			if(!ProductionRun && BaseEngine::local_rank==0) LOGGERA("SEGMENT EXCEEDED MaximumSegmentLength BLOCKS. EXITING.")
			if(BaseEngine::local_rank==0) LOGGER("SEGMENT EXCEEDED MaximumSegmentLength BLOCKS. EXITING.")
			segment.initialE = reference.getEnergy();
			segment.finalE = currentMin.getEnergy();
			segment.OutOfBasin = false;
			segment.transition.first = InitialLabels;
			segment.transition.second = CurrentLabels;

			break;
		}
	}
	//if(LabelError) CurrentLabels = InitialLabels;
	//segment.transition.second = CurrentLabels;
	insert("TADSegment",task.returns,segment);
	task.clearInputs();
	return ;
};

/**
 * Implement NEB routine between two states
 *
 * This expects:
 * inputData: Initial, Final, Saddle (optional)
 * arguments: pathway: [required: (InitialLabels, FinalLabels) , optional: SaddleLabels]
 *
 * ExistingPairs are pairs of states with the same initial and final canonical  labels
 * This returns
 * outputData: Saddle (if converged) else FoundTransitions
 * returns: pathway
 */
 std::function<void(GenericTask&)> neb_impl = [this](GenericTask &task) {

 	task.clearOutputs(); // just in case

 	std::unordered_map<std::string,std::string> parameters=\
 		extractParameters(task.type,task.flavor,BaseMDEngine::defaultFlavor,BaseMDEngine::taskParameters);

 	//read parameters from the task
 	double dt=safe_extractor<double>(parameters,"NEBTimestep",0.001);
	double WellDepth=safe_extractor<double>(parameters,"WellDepth",0.5);
	int nImages=safe_extractor<int>(parameters,"Images",11);
	int maxiter=safe_extractor<int>(parameters,"MaxIterations",1000);
 	double spring=safe_extractor<double>(parameters,"Spring",1.0);
	double ftol=safe_extractor<double>(parameters,"ForceTolerance",0.01);
 	bool doClimb=safe_extractor<bool>(parameters,"Climbing",false);
 	bool writeFiles=safe_extractor<bool>(parameters,"WriteFiles",false);
	int clust_thresh=safe_extractor<int>(parameters,"NEBClusterThresh",-1);
	bool CalculatePrefactor=safe_extractor<bool>(parameters,"CalculatePrefactor",false);
	double ThresholdBarrier=safe_extractor<double>(parameters,"ThresholdBarrier",1.0);

 	// GenericTasks
 	GenericTask label;
 	label.type=BaseEngine::mapper.type("TASK_LABEL");
 	label.flavor=task.flavor;

 	GenericTask min;
 	min.type=BaseEngine::mapper.type("TASK_MIN");
 	min.flavor=task.flavor;

 	GenericTask carve;
	carve.type=BaseEngine::mapper.type("TASK_CARVE");
	carve.flavor = task.flavor;

 	GenericTask symmetry;
 	symmetry.type=BaseEngine::mapper.type("TASK_SYMMETRY");
 	symmetry.flavor=task.flavor;

 	GenericTask write;
 	write.type=BaseEngine::mapper.type("TASK_WRITE_TO_FILE");
 	write.flavor=task.flavor;
 	std::string filename;

 	//Extract Systems
 	System initial, final, saddle;
 	std::list<System> ExistingPairs;


 	LabelPair InitialLabels, FinalLabels, saddle_labels;
 	//std::list<LabelPair> ExistingPair_labels;

 	NEBPathway pathway; // home of all non-configuration return data;
 	bool success, reset, have_saddle, have_pairs, have_pathway,fsuccess;
	int InitialClusters=1,FinalClusters=1;
 	have_pathway = extract("NEBPathway",task.arguments,pathway);


 	std::list<TransitionSymmetry> transitions;
 	std::map< Label, std::set<PointShiftSymmetry> > self_symmetries;

 	success = extract("Initial",task.inputData,initial);
 	fsuccess = extract("Final",task.inputData,final);
 	have_saddle = extract("Saddle",task.inputData,saddle);

 	if(!success) if(BaseEngine::local_rank==0) LOGGERA("COULDN'T FIND INITIAL STATES! EXITING!")
 	if(!fsuccess) if(BaseEngine::local_rank==0) LOGGERA("COULDN'T FIND FINAL STATES! EXITING!")

 	if(!success || !fsuccess) {
 		pathway.valid=false;
 		pathway.saddleE = MAX_BARRIER;
 		pathway.SaddleLabels = pathway.InitialLabels;
 		insert("NEBPathway",task.returns,pathway);
 		return;
 	}

 	if(have_saddle) {
 		label.clearInputs(); label.clearOutputs();
 		insert("State",label.inputData,final);
 		BaseEngine::process(label);
 		extract("Labels",label.returns,saddle_labels);
 		pathway.SaddleLabels = saddle_labels;
 	}
 	// not implemented yet...
 	pathway.priornu = std::make_pair(PRIOR_NU,PRIOR_NU);


 	// minimize everything (found to be useful; really needed?)
 	min.clearOutputs(); min.clearInputs();
 	insert("State",min.inputData,initial);
 	BaseEngine::process(min);
 	extract("State",min.outputData,initial);

 	label.clearInputs(); label.clearOutputs();
 	insert("State",label.inputData,initial);
 	BaseEngine::process(label);
 	extract("Labels",label.returns,InitialLabels);

	carve.clearInputs(); carve.clearOutputs();
 	insert("State",carve.inputData,initial);
 	BaseEngine::process(carve);
 	extract("Clusters",carve.returns,InitialClusters);

 	if(BaseEngine::local_rank==0) LOGGER("Initial E: "<<initial.getEnergy())
 	pathway.initialE = initial.getEnergy();


 	min.clearOutputs(); min.clearInputs();
 	insert("State",min.inputData,final);
 	BaseEngine::process(min);
 	extract("State",min.outputData,final);

 	label.clearInputs(); label.clearOutputs();
 	insert("State",label.inputData,final);
 	BaseEngine::process(label);
 	extract("Labels",label.returns,FinalLabels);

	carve.clearInputs(); carve.clearOutputs();
 	insert("State",carve.inputData,final);
 	BaseEngine::process(carve);
 	extract("Clusters",carve.returns,FinalClusters);

 	if(BaseEngine::local_rank==0) LOGGER("Final E: "<<final.getEnergy())
 	pathway.finalE = final.getEnergy();

	pathway.mismatch = false;

 	if(have_pathway and (pathway.InitialLabels!=InitialLabels or pathway.FinalLabels != FinalLabels)) {
 		if(pathway.InitialLabels!=InitialLabels) {
 			if(BaseEngine::local_rank==0) LOGGER("INITIAL LABEL MISMATCH: "<<pathway.InitialLabels.first<<","<<pathway.InitialLabels.second<<" != "<<InitialLabels.first<<","<<InitialLabels.second)
 		} else {
 			if(BaseEngine::local_rank==0) LOGGER("INITIAL LABEL MATCH: "<<pathway.InitialLabels.first<<","<<pathway.InitialLabels.second<<" == "<<InitialLabels.first<<","<<InitialLabels.second)
 		}
 		if(pathway.FinalLabels != FinalLabels) {
 			if(BaseEngine::local_rank==0) LOGGER("FINAL LABEL MISMATCH: "<<pathway.FinalLabels.first<<","<<pathway.FinalLabels.second<<" != "<<FinalLabels.first<<","<<FinalLabels.second)
 		} else {
 			if(BaseEngine::local_rank==0) LOGGER("FINAL LABEL MATCH: "<<pathway.FinalLabels.first<<","<<pathway.FinalLabels.second<<" == "<<FinalLabels.first<<","<<FinalLabels.second)
 		}

		InitialLabels = pathway.InitialLabels;
		FinalLabels = pathway.FinalLabels;


 		if((pathway.InitialLabels==FinalLabels) and (pathway.FinalLabels==InitialLabels)) {
 			if(BaseEngine::local_rank==0) LOGGER("STRANGE SWAPPING BUG IN std::list<System>. SWAPPING BACK")
 			System temp_S = initial;
 			initial = final;
 			final = temp_S;

 			InitialLabels = pathway.InitialLabels;
 			FinalLabels = pathway.FinalLabels;
 			pathway.initialE = initial.getEnergy();
 			pathway.finalE = final.getEnergy();
 		} else {
 			pathway.mismatch = true;
 			pathway.MMInitialLabels = InitialLabels;
 			pathway.MMFinalLabels = FinalLabels;
 			if(BaseEngine::local_rank==0) LOGGER("STRANGE MISMATCH BUG IN DB. EXITING AND TRYING AGAIN....")
 			pathway.valid=false;
 			pathway.saddleE = MAX_BARRIER + std::max(pathway.initialE,pathway.finalE);
			insert("NEBPathway",task.returns,pathway);
 			return;
 		}
 	}

 	pathway.valid=true; // innocent until proven guilty...
	if(BaseEngine::local_rank==0) LOGGER("InitialClusters: "<<InitialClusters<<" FinalClusters: "<<FinalClusters)
	if(clust_thresh>0 and std::max(InitialClusters,FinalClusters)>clust_thresh) {
		pathway.valid=false; // innocent until proven guilty...
		pathway.saddleE = pathway.initialE + MAX_BARRIER;
 		insert("NEBPathway",task.returns,pathway);
		if(BaseEngine::local_rank==0) LOGGERA("TOO MANY CLUSTERS; EXITING")
		return;
	}

 	if(BaseEngine::local_rank==0) LOGGER("CANONICAL TRANSITION: "<<InitialLabels.first<<" -> "<<FinalLabels.first)

 	// add initial and final states to symmetry
 	symmetry.clearInputs();

 	#ifdef ISOMORPHIC

 	for(auto ss: pathway.InitialSymmetries) self_symmetries[pathway.InitialLabels.first].insert(ss);
 	for(auto ss: pathway.FinalSymmetries) self_symmetries[pathway.FinalLabels.first].insert(ss);

 	insert("Targets",InitialLabels.first,InitialLabels.second,0,false,symmetry.inputData,initial);
 	insert("Targets",FinalLabels.first,FinalLabels.second,0,false,symmetry.inputData,final);

 	insert("SelfSymmetries",symmetry.arguments,self_symmetries);

 	// minimize all ExistingPairs then compare transitions
 	extract("ExistingPairs",task.inputData,ExistingPairs);
 	have_pairs = (ExistingPairs.size()>0) and (ExistingPairs.size()%2==0);

 	if (have_pairs) {
 		if(BaseEngine::local_rank==0) LOGGER("Found "<<ExistingPairs.size()/2<<" pairs")
 		auto sys_ep = ExistingPairs.begin();
 		LabelPair labels;

 		while(sys_ep!=ExistingPairs.end()) {
 			min.clearOutputs(); min.clearInputs();
 			insert("State",min.inputData,*sys_ep);
 			BaseEngine::process(min);
 			extract("State",min.outputData,*sys_ep);

 			// better to check?
 			label.clearOutputs(); label.clearInputs();
 			insert("State",label.inputData,*sys_ep);
 			BaseEngine::process(label);
 			extract("Labels",label.returns,labels);

 			insert("Candidates",labels.first,labels.second,0,false,symmetry.inputData,*sys_ep);
 			if(BaseEngine::local_rank==0) LOGGER("Adding ExistingPairs: ("<<labels.first<<","<<labels.second<<")")

 			sys_ep = ExistingPairs.erase(sys_ep); // delete from master; added to symmetry
 			//labels = std::next(labels);
 		}
 	}


 	BaseEngine::process(symmetry);

 	// extract symmetry returns
 	extract("SelfTransitions",symmetry.returns,pathway.self_transitions);
 	extract("StateIsomorphisms",symmetry.returns,pathway.equivalent_states);
 	extract("TransitionIsomorphisms",symmetry.returns,pathway.equivalent_transitions);

 	self_symmetries.clear();
 	extract("SelfSymmetries",symmetry.returns,self_symmetries);
 	for(auto ss: self_symmetries) {
 		if(ss.first==pathway.InitialLabels.first)
 			for(auto sss: ss.second) pathway.InitialSymmetries.insert(sss);
 		if(ss.first==pathway.FinalLabels.first)
 			for(auto sss: ss.second) pathway.FinalSymmetries.insert(sss);
 	}

 	if (pathway.equivalent_transitions.size() > 0) {
 		pathway.duplicate = true;
 		if(BaseEngine::local_rank==0) LOGGER("Matches found, exiting NEB")
 		task.clearInputs();
 		insert("NEBPathway",task.returns,pathway);
 		if(writeFiles) {
 			if(BaseEngine::local_rank==0) LOGGER("Writing Initial and Final states")
 			// write initial and final states....
 			write.clearInputs(); write.clearOutputs();
 			filename="states/state-"+std::to_string(pathway.InitialLabels.first)+"-"+std::to_string(pathway.InitialLabels.second)+".dat";
 			insert("Filename",write.arguments,filename);
 			insert("State",write.inputData,initial);
 			BaseEngine::process(write);

 			write.clearInputs(); write.clearOutputs();
 			filename="states/state-"+std::to_string(pathway.FinalLabels.first)+"-"+std::to_string(pathway.FinalLabels.second)+".dat";
 			insert("Filename",write.arguments,filename);
 			insert("State",write.inputData,final);
 			BaseEngine::process(write);
 		}
 		return;
 	} else {
 		if(BaseEngine::local_rank==0) LOGGER("No matches found, Starting NEB")
 		pathway.duplicate = false;
 	}
 	#endif

 	// NEB routine proper

 	// wrap final to initial
 	initial.minimumImage(final); // final -> initial + minimumImageVector(final-initial)
 	pathway.dXmax = initial.msd(final,true);
 	pathway.dX = initial.msd(final,false);
 	if(BaseEngine::local_rank==0) LOGGER("dXThresh: "<<pathway.dXmax<<" "<<pathway.dX)

 	// if gap too large return with MAX_BARRIER
 	if (pathway.dXmax>10.0) {
 		pathway.valid=false;
 		if(BaseEngine::local_rank==0) LOGGER("max|dX| > 10A!!")
 		pathway.SaddleLabels = pathway.InitialLabels;
 		pathway.saddleE = pathway.initialE + 2.0 * MAX_BARRIER;
 		pathway.Ftol = 100.0;
 		task.clearInputs();
 		insert("NEBPathway",task.returns,pathway);
 		if(writeFiles) {
 			if(BaseEngine::local_rank==0) LOGGER("Writing Initial and Final states")
 			// write initial and final states....
 			write.clearInputs(); write.clearOutputs();
 			filename="states/state-"+std::to_string(pathway.InitialLabels.first)+"-"+std::to_string(pathway.InitialLabels.second)+".dat";
 			insert("Filename",write.arguments,filename);
 			insert("State",write.inputData,initial);
 			BaseEngine::process(write);

 			write.clearInputs(); write.clearOutputs();
 			filename="states/state-"+std::to_string(pathway.FinalLabels.first)+"-"+std::to_string(pathway.FinalLabels.second)+".dat";
 			insert("Filename",write.arguments,filename);
 			insert("State",write.inputData,final);
 			BaseEngine::process(write);
 		}
 		return;
 	}
 	// if gap too small return with MIN_BARRIER ??

 	// initialize f and set v=0
 	const int nAtoms = initial.getNAtoms();
 	initial.fflag = 1; final.fflag = 1;
 	initial.resize(nAtoms); final.resize(nAtoms);
 	for(int i=0;i<NDIM*nAtoms;i++) initial.v[i] = 0.0;
 	for(int i=0;i<NDIM*nAtoms;i++) final.v[i] = 0.0;
 	ExecutionTimer timer;

 	std::vector<System> neb_systems;
 	std::vector<double> energies;
 	energies.push_back(initial.getEnergy());
 	if(have_saddle) {
 		energies.push_back(saddle.getEnergy());
 		neb_systems.push_back(saddle);
 	}
 	energies.push_back(final.getEnergy());
 	if(!have_saddle) interpolate(neb_systems,initial,final,0,nImages-2,energies);
 	else {
 		interpolate(neb_systems,initial,final,1,nImages/2-1,energies);
 		interpolate(neb_systems,initial,final,0,nImages-2-nImages/2,energies);
 	}

 	neb_forces(neb_systems,initial,final,energies); // no projection here

 	int i=0;
 	if(BaseEngine::local_rank==0) for( auto E: energies) LOGGER(i++<<" "<<E-initial.getEnergy())


 	// Initializations for FIRE
 	bool rcNN = true; // recompute nn list etc
 	bool Climbing = false;
 	bool Extended = false;
 	bool ReDiscretize = false;
 	int ClimbingImage = nImages / 2;
 	std::list<int> InterMinImages;
 	int MinImage = 0;
 	double c_v,c_f,v_n,f_n,x_at,v_at,f_at,vdotf,f_ratio=1.0;

 	// FIRE parameters- TMAX=20 quite aggressive but seems OK
 	const double AL_ST = 0.1, AL_SHR = 0.99, DT_SHR = 0.5, DT_GR = 1.3, TMAX = 100.;
 	const int DELAYSTEP = 5;
 	const double dtmax = TMAX * dt;

 	//std::vector<double> dt_a(nImages,dt), alpha(nImages,AL_ST);
 	//std::vector<int> last_negative(nImages,0);
 	int lnt=0;
 	double dt_at=dt,alt=AL_ST;
 	double max_x_disp=0.0,max_x_c_disp=0.0;
 	double max_f_at_sq = 10.0 * ftol * ftol;
 	timer.start("NEB_ROUTINE");
 	if(BaseEngine::local_rank==0) LOGGERA("NEB STEP: MinImage ClimbingImage InterMinImages? E[ClimbingImage]-E[0] E[MinImage]-E[0] sqrt(max_f_at_sq) sqrt(max_x_c_disp) dt")

 	pathway.FoundTransitions.clear();
 	for(int iter = 0; iter < maxiter; iter++) {
 		timer.start("NEB_ITERATION");

 		if(max_x_c_disp>.1*.1 || max_f_at_sq>1.0*1.0) { // skin depth is 2A, so 1A conservative....
 			if(BaseEngine::local_rank==0) LOGGER("Relisting neighbors due to drift: |dX|_inf="<<sqrt(max_x_c_disp)<<", |F|_inf="<<sqrt(max_f_at_sq))
 			rcNN=true;
 			max_x_c_disp=0.0;
			if(max_f_at_sq>50.0*50.0) {
				if(BaseEngine::local_rank==0) LOGGER("Fmax is very large! Quitting this NEB routine and flagging for restart!")
				if(BaseEngine::local_rank==0) LOGGERA("eV/A -> EXIT NEB")
				pathway.valid=false;
 				pathway.SaddleLabels = pathway.InitialLabels;
		 		pathway.saddleE = pathway.initialE + 2.0 * MAX_BARRIER;
		 		pathway.Ftol = 100.0;
 				task.clearInputs();
		 		insert("NEBPathway",task.returns,pathway);
				return;
			}
 		} else rcNN=bool(iter%20==0); // always every 20 steps in either case

 		rcNN=true; // safer for now...
 		neb_forces(neb_systems,initial,final,energies,false,true,spring,ClimbingImage,Climbing,rcNN); // no reset this time, but projection


 		path_analysis(neb_systems,initial,final,energies,MinImage,ClimbingImage,InterMinImages,WellDepth);

 		// we have intermediate minima...
 		if ((InterMinImages.size()>0) && ((max_f_at_sq<4.0*ftol*ftol)||(iter==maxiter-1))) {
 			ReDiscretize = true; // only happens once..
 			if(BaseEngine::local_rank==0) LOGGERA("NEB "<<iter<<": "<<MinImage<<" "<<ClimbingImage<<" "<<bool(InterMinImages.size()>0)<<" "<<energies[ClimbingImage]-energies[0]<<" "<<energies[MinImage]-energies[0]<<" "<<sqrt(max_f_at_sq)<<" "<<sqrt(max_x_c_disp))
 			if(BaseEngine::local_rank==0) LOGGERA("InterMinImage(s) detected! Profile: E[0]="<<energies[0])
 			double dE = energies[ClimbingImage]-energies[0];
 			// max res: 20 spaces
			int min_im=0;
			for(int im=1;im<nImages;im++) if(energies[im]<energies[min_im]) min_im=im;
			for(int im=0;im<nImages;im++) {
 				std::string res = "\t";
 				for(int ss=0;ss<int(20.0*(energies[im]-energies[min_im])/dE);ss++) res+=" ";
 				res += "| "+std::to_string(energies[im]-energies[min_im]);
 				for(auto InterMinImage: InterMinImages)	if(im==InterMinImage) res+=" InterMinImage";
 				if(BaseEngine::local_rank==0) LOGGERA(res)
 			}

 			Transition ntrans;
 			ntrans.first = pathway.InitialLabels;
 			for(auto InterMinImage: InterMinImages) {
 				min.clearOutputs();
 				min.clearInputs();
 				insert("State",min.inputData,neb_systems[InterMinImage-1]);
 				BaseEngine::process(min);
 				extract("State",min.outputData,saddle);

 				label.clearOutputs();
 				label.clearInputs();
 				insert("State",label.inputData,saddle);
 				BaseEngine::process(label);
 				extract("Labels",label.returns,ntrans.second);

 				pathway.FoundTransitions.push_back(ntrans);
 				if(BaseEngine::local_rank==0) LOGGER("INSERTING "<<ntrans.second.first<<","<<ntrans.second.second)
 				insert("State",ntrans.second.first,ntrans.second.second,LOCATION_SYSTEM_MIN,true,task.outputData,saddle);
 				if(writeFiles) {
 					if(BaseEngine::local_rank==0) LOGGER("Writing Intermediate states")
 					// write initial and final states....
 					write.clearInputs(); write.clearOutputs();
 					filename="states/state-"+std::to_string(ntrans.second.first)+"-"+std::to_string(ntrans.second.second)+".dat";
 					insert("Filename",write.arguments,filename);
 					insert("State",write.inputData,neb_systems[InterMinImage-1]);
 					BaseEngine::process(write);
 				}
 				ntrans.first = ntrans.second;
 			}
 			ntrans.second = pathway.FinalLabels;
 			pathway.FoundTransitions.push_back(ntrans);

 			pathway.valid=false;
 			pathway.saddleE = energies[ClimbingImage];
 			pathway.energies = energies;
 			pathway.Ftol = 100.0;//sqrt(max_f_at_sq);
 			insert("NEBPathway",task.returns,pathway);
 			if(BaseEngine::local_rank==0) LOGGER("Exiting NEB!")

 			if(writeFiles) {
 				if(BaseEngine::local_rank==0)
					LOGGER("Writing Initial and Final states")
 				// write initial and final states....
 				write.clearInputs(); write.clearOutputs();
 				filename="states/state-"+std::to_string(pathway.InitialLabels.first)+"-"+std::to_string(pathway.InitialLabels.second)+".dat";
 				insert("Filename",write.arguments,filename);
 				insert("State",write.inputData,initial);
 				BaseEngine::process(write);

 				write.clearInputs(); write.clearOutputs();
 				filename="states/state-"+std::to_string(pathway.FinalLabels.first)+"-"+std::to_string(pathway.FinalLabels.second)+".dat";
 				insert("Filename",write.arguments,filename);
 				insert("State",write.inputData,final);
 				BaseEngine::process(write);
 			}
 			return;
 		}

 		if(iter%10==0 && BaseEngine::local_rank==0)
			LOGGERA("NEB "<<iter<<": "<<MinImage<<" "<<ClimbingImage<<" "<<bool(InterMinImages.size()>0)<<" "<<energies[ClimbingImage]-energies[0]<<" "<<energies[MinImage]-energies[0]<<" "<<sqrt(max_f_at_sq)<<" "<<sqrt(max_x_c_disp))

 		max_f_at_sq = 0.0; vdotf=0.0; f_n=0.0; v_n=0.0; f_ratio=1.0;
		//for(int l=0;l<nImages-2;l++)
		for(int i=0;i<nAtoms*NDIM;i++) {
 			f_at = neb_systems[std::max(ClimbingImage-1,0)].getForce(i/NDIM,i%NDIM);
 			if(f_at*f_at>max_f_at_sq) max_f_at_sq = f_at*f_at;
 		}


 		// conditioning for minimization
 		if(max_f_at_sq > 2.0*2.0) f_ratio = 2.0/sqrt(max_f_at_sq);

 		for(int l=0;l<nImages-2;l++) for(int i=0;i<nAtoms*NDIM;i++) {
 			x_at = neb_systems[l].getPosition(i/NDIM,i%NDIM);
 			v_at = neb_systems[l].getVelocity(i/NDIM,i%NDIM);
 			f_at = neb_systems[l].getForce(i/NDIM,i%NDIM) * f_ratio;
 			if(f_at*f_at>max_f_at_sq) max_f_at_sq = f_at*f_at;
 			f_n += f_at * f_at;
 			v_n += v_at * v_at;
 			vdotf += v_at * f_at;
 		}
 		f_n = sqrt(f_n);
 		v_n = sqrt(v_n);
  		if ( vdotf > 0.0) {
 	 		// if (v dot f) > 0:
 	 		// v = (1-alpha) v + alpha |v| Fhat
 	 		// |v| = length of v, Fhat = unit f
 	 		// no verlet here, simple to implement though
 	 		c_v = -alt;
 	 		if (f_n == 0.0) c_f = 0.0;
 	 		else c_f = alt * v_n/f_n;
 			max_x_disp=0.0;
 			for(int l=0;l<nImages-2;l++) for(int i=0;i<nAtoms*NDIM;i++) {
 				x_at = neb_systems[l].getPosition(i/NDIM,i%NDIM);
 				v_at = neb_systems[l].getVelocity(i/NDIM,i%NDIM);
 				f_at = neb_systems[l].getForce(i/NDIM,i%NDIM);
 				neb_systems[l].setVelocity(i/NDIM,i%NDIM,v_at + c_v*v_at + c_f*f_at);
 				neb_systems[l].setPosition(i/NDIM,i%NDIM,x_at + v_at * dt_at);
 				max_x_disp = std::max(max_x_disp,v_at * dt_at * v_at * dt_at);
 			}
 			max_x_c_disp += max_x_disp;
 			// if more than DELAYSTEP since v dot f was negative:
 	 		// increase timestep and decrease alpha
  			if(iter - lnt > DELAYSTEP) {
  				dt_at = std::min(dtmax,dt_at*DT_GR);
  				alt *= AL_SHR;
 			}
 		} else {
  			// else (v dot f) <= 0:
  			// decrease timestep, reset alpha, set v = dt * f
  			// FIRE says v=0 but this gets stuck
  			lnt = iter;
  			dt_at *= DT_SHR;
  			alt = AL_ST;
 			for(int l=0;l<nImages-2;l++) for(int i=0;i<nAtoms*NDIM;i++)
 				neb_systems[l].v[i] = dt_at * neb_systems[l].f[i];
 		}

 		if (max_f_at_sq < 4.0 * ftol * ftol) {
  		 if (!Climbing and doClimb) {
  			 if(BaseEngine::local_rank==0) LOGGER("CLIMBING"<<std::endl)
  			 Climbing = true;
  		 } else if(max_f_at_sq<ftol*ftol) break;
  	 }
 	 timer.stop("NEB_ITERATION");
 	}
 	timer.stop("NEB_ROUTINE");

 	i=0;
 	if(BaseEngine::local_rank==0) {
		for( auto E: energies) LOGGERA("E["<<i++<<"]-E[0]="<<E-initial.getEnergy())
	 	LOGGERA(energies[ClimbingImage]-energies[0]<<" "<<neb_systems.size()<<" "<<ClimbingImage)
	}

 	// Add saddle to system
	if(ClimbingImage>=neb_systems.size()) ClimbingImage = neb_systems.size()-1;
	if(ClimbingImage<0) ClimbingImage = 0;
 	saddle = neb_systems[ClimbingImage];
 	label.clearInputs(); label.clearOutputs();
 	insert("State",label.inputData,saddle);
 	BaseEngine::process(label);
 	extract("Labels",label.returns,pathway.SaddleLabels);
 	insert("State",pathway.SaddleLabels.first,pathway.SaddleLabels.second,LOCATION_SYSTEM_SADDLE,true,task.outputData,saddle);
 	pathway.saddleE = energies[ClimbingImage];
 	pathway.energies = energies;
 	pathway.Ftol = sqrt(std::fabs(max_f_at_sq));
 	insert("NEBPathway",task.returns,pathway);
 	task.clearInputs();
 	timer.report();


 	if(writeFiles) {
 		if(BaseEngine::local_rank==0) LOGGER("Writing Initial, Final and Saddle states")
 		// write initial and final states....
 		write.clearInputs(); write.clearOutputs();
 		filename="states/state-"+std::to_string(pathway.InitialLabels.first)+"-"+std::to_string(pathway.InitialLabels.second)+".dat";
 		insert("Filename",write.arguments,filename);
 		insert("State",write.inputData,initial);
 		BaseEngine::process(write);

 		write.clearInputs(); write.clearOutputs();
 		filename="states/state-"+std::to_string(pathway.SaddleLabels.first)+"-"+std::to_string(pathway.SaddleLabels.second)+".dat";
 		insert("Filename",write.arguments,filename);
 		insert("State",write.inputData,saddle);
 		BaseEngine::process(write);

 		write.clearInputs(); write.clearOutputs();
 		filename="states/state-"+std::to_string(pathway.FinalLabels.first)+"-"+std::to_string(pathway.FinalLabels.second)+".dat";
 		insert("Filename",write.arguments,filename);
 		insert("State",write.inputData,final);
 		BaseEngine::process(write);
 	}


 	return;
 };

/**
 * read a configuration from a file, minimize, and label it.
 *
 * This expects:
 *

 * arguments: Filename
 *
 * This returns
 * returns: Labels, Clusters, Position, SelfSymmetries
 * outputData State
 */
std::function<void(GenericTask&)> init_min_impl = [this](GenericTask &task) {

	task.clearOutputs();

	LabelPair labels;
	GenericTask t;
	t.clearOutputs();
	t.arguments=task.arguments;
	System s;
	t.flavor=task.flavor;

	//read the file
	t.type=BaseEngine::mapper.type("TASK_INIT_FROM_FILE");

	BaseEngine::process(t);
	extract("State",t.outputData,s);

	// minimize
	t.clear();
	t.type=BaseEngine::mapper.type("TASK_MIN");
	insert("State",t.inputData,s);
	BaseEngine::process(t);
	extract("State",t.outputData,s);

	double energy = s.getEnergy();
	//extract("Energy",t.returns,energy);
	insert("Energy",task.returns,energy);

	// label it
	t.clear();
	t.type=BaseEngine::mapper.type("TASK_LABEL");
	insert("State",t.inputData,s);
	BaseEngine::process(t);
	extract("State",t.outputData,s);
	extract("Labels",t.returns,labels);
	insert("Labels",task.returns,labels);

	// carve it
	t.clear();
	t.type=BaseEngine::mapper.type("TASK_CARVE");
	insert("State",t.inputData,s);
	BaseEngine::process(t);

	// dump it
	int clusters=1;
	std::array<double,3> position={0.,0.,0.};
	extract("Clusters",t.returns,clusters);
	extract("Position",t.returns,position);

	insert("Clusters",task.returns,clusters);
	insert("Position",task.returns,position);

	#ifdef ISOMORPHIC

	std::unordered_map<std::string,std::string> parameters =
 		extractParameters(BaseEngine::mapper.type("TASK_NEB"),task.flavor,
			BaseMDEngine::defaultFlavor,BaseMDEngine::taskParameters);


	std::map<Label,std::set<PointShiftSymmetry>> self_symmetries;
	t.clear();
	t.type=BaseEngine::mapper.type("TASK_SYMMETRY");
 	insert("Targets",labels.first, labels.second,0,false,t.inputData,s);
	BaseEngine::process(t);
	extract("SelfSymmetries",t.returns,self_symmetries);
	if(self_symmetries.find(labels.first) != self_symmetries.end()) {
		insert("SelfSymmetries",task.returns,self_symmetries[labels.first]);
	}
	#endif

	insert("State",labels.first, labels.second, LOCATION_SYSTEM_MIN, true, task.outputData, s);

};

std::function<void(GenericTask&)> label_impl = [this](GenericTask &task) {
	bool canonical=false;
	System s;
	extract("State",task.inputData,s);
	LabelPair labels;
	labels.first = BaseMDEngine::labeler->hash(s,true);
	labels.second = BaseMDEngine::labeler->hash(s,false);
	if(BaseEngine::local_rank==0)
		LOGGER("TADEngine::label_impl: "<<labels.first<<" "<<labels.second)
	insert("State",labels.first,labels.second,0,false,task.outputData,s);
	insert("Labels",task.returns,labels);
	insert("Label",task.returns,labels.second);
	insert("CanonicalLabel",task.returns,labels.first);
};

//joint remapping of a number of states. This assumes that all systems in the task correspond to the same state. The first state is used as a reference
std::function<void(GenericTask&)> remap_impl = [this](GenericTask &task){
	System csys,sys;
	std::map<int,int> canonicalMap;
	LabelPair labs;
	Label clab;
	extract("State",task.inputData,sys);
	if (extract("State",task.inputData,csys)) {
		labs.first = BaseMDEngine::labeler->hash(sys,false);
		BaseMDEngine::labeler->canonicalMap(csys,canonicalMap,labs.second);
		csys.remap(canonicalMap);
		clab = BaseMDEngine::labeler->hash(csys,false);
		insert("Labels",task.returns,labs);
		insert("State",labs.first,labs.second,LOCATION_SYSTEM_MIN,false,task.outputData,sys);
		insert("CanonState",clab,labs.second,LOCATION_SYSTEM_MIN,false,task.outputData,csys);
	}
};

/*   HELPER FUNCTIONS */

/*
	If final derived engine has singleForceEnergyCall(..) this is overwritten
	Bit of a hack but stable
*/
virtual void singleForceEnergyCall(System &s, bool noforce=false,bool prepost=true) {
	GenericTask forces;
	forces.type = BaseEngine::mapper.type("TASK_FORCES");
	bool reset = false;
	insert("Reset",forces.arguments,reset);
	insert("State",forces.inputData,s);
	BaseEngine::process(forces);
	extract("State",forces.outputData,s);
};

void neb_forces(std::vector<System> &nsv, System &initial, System &final, std::vector<double> &ev, bool reset=true, bool project=false, double spring=1.0, int ClimbingImage=1, bool Climbing=false, bool prepost=true) {

  /* if(fc) {
		GenericTask forces;
		forces.type = BaseEngine::mapper.type("TASK_FORCES");
		insert("Reset",forces.arguments,reset);
		insert("PrePost",forces.arguments,prepost);
		for(auto &ns: nsv) insert("State",forces.inputData,ns);
		BaseEngine::process(forces);

		std::list<System> nsl;
		extract("State",forces.outputData,nsl);
		std::swap_ranges(nsl.begin(), nsl.end(), nsv.begin());
	} else */

	for(auto &ns: nsv) singleForceEnergyCall(ns,false,prepost); // now robust (see above)

	int ei=1;
	for(auto sit=nsv.begin();sit!=nsv.end();sit++,ei++) ev[ei] = sit->getEnergy();

	// do the NEB projection. Following Henkelman et al JCP 2000
	if (project) {
		Cell bc = initial.getCell();
		int nImages=ev.size(), nAtoms=initial.getNAtoms();
		System *sysp;
		std::vector<double> ft(NDIM*nAtoms,0.0), bt(NDIM*nAtoms,0.0);
		double ft_n,bt_n,t_n,mm,ediff[2],tmix[2],tdotf;
		for(int l=0;l<nImages-2;l++) {
			// FT = (X_i+1 - X_i) , |X_i+1 - X_i|
			sysp = ( l==nImages-3 ? &final : &nsv[l+1] );
			bc.minimumImageVector(nsv[l].x,sysp->x,ft,ft_n,mm); ft_n = sqrt(ft_n);

			// BT = (X_i - X_i-1) , |X_i - X_i-1|
			sysp = ( l==0 ? &initial : &nsv[l-1] );
			bc.minimumImageVector(sysp->x,nsv[l].x,bt,bt_n,mm); bt_n = sqrt(bt_n);

			// ft will become true tangent vector
			t_n = ft_n;

			// look at curvature
			ediff[0] = ev[l]-ev[l-1];
			ediff[1] = ev[l+1]-ev[l];

			if(ediff[0]<0.0 and ediff[1]<0.0) { // i.e. both same sign => not at saddle but -ve slope
				t_n = bt_n;
				ft = bt;
			} else { // at 'saddle' so mix
				tmix[0] = (ediff[1]+ediff[0]>0.0 ? std::min(fabs(ediff[0]),fabs(ediff[1])) : std::max(fabs(ediff[0]),fabs(ediff[1])));
				tmix[1] = (ediff[1]+ediff[0]>0.0 ? std::max(fabs(ediff[0]),fabs(ediff[1])) : std::min(fabs(ediff[0]),fabs(ediff[1])));
				t_n = 0.0;
				for(int i=0; i<nAtoms*NDIM; i++) {
					ft[i] = ft[i]*tmix[0] + bt[i]*tmix[1];
					t_n += ft[i] * ft[i];
				}
				t_n = sqrt(t_n);
			}
			for(int i=0; i<nAtoms*NDIM; i++) ft[i] /= t_n; // normalize for projection
			// F_i . T_i
			tdotf = 0.0; for(int i=0; i<nAtoms*NDIM; i++) tdotf += nsv[l].f[i] * ft[i];
			// Climbing: F_i  = F_i - 2T_i  x  [ (F_i.T_i) ]
			// Regular:  F_i  = F_i -  T_i  x  [ (F_i.T_i) - K (|R_{i+1}-R_i| - |R_i-R_{i-1}|) ]
			if(l==ClimbingImage and Climbing) for(int i=0; i<nAtoms*NDIM; i++) nsv[l].f[i] -= 2.0 * ft[i] * tdotf;
			else for(int i=0; i<nAtoms*NDIM; i++) nsv[l].f[i] -= ft[i] * ( tdotf - spring * ( ft_n - bt_n ) );
		}
	}
};

/* inserts nImages Systems from index start-1 of neb_systems. ie start=0 implies index */
virtual void interpolate(std::vector<System> &nsv,System &initial,System &final, unsigned int start,unsigned int nImages,std::vector<double> &energies) {

	double r=0.0, dr=1.0/(double)(nImages+1);
	int n = initial.getNAtoms();
	std::vector<double> dx;
	Cell bc = initial.getCell(); // assumes constant over path

	if(start>nsv.size()) start = nsv.size();

	auto eit = std::next(energies.begin(),start+1); // always shifted by one more
	double refe = energies[start], drefe=energies[start+1]-energies[start]; // always possible
	auto it = std::next(nsv.begin(),start);

	System *a = (start==0 ? &initial : &nsv[start-1]);
	System *b = (start==nsv.size() ? &final : &nsv[start]);
	//a->minimumImage(*b);
	bc.minimumImageVector(a->x,b->x,dx);
	for(int j=0;j<NDIM*n;j++) b->x[j] = a->x[j] + dx[j];

	for(int imi=0;imi<nImages;imi++,it++,eit++) {
		r = dr * (imi+1);
		it = nsv.insert(it,*a);
		eit = energies.insert(eit, refe + r * drefe);
		for(int j=0;j<NDIM*n;j++) it->x[j] += dx[j] * r;
	}

};

/* finds minImage and either one or (if minImage not at an end) two climbing images */
virtual void path_analysis(std::vector<System> &systems,System &initial,System &final,std::vector<double> &energies,int &MinImage,\
	 int &ClimbingImage, std::list<int> &InterMinImages, double e_thresh) {

	// Finds minimia and climbing image of an energy path
	int lastimg = energies.size()-1;

	double tempE = energies[0];

	// MinImage that is not 0
	tempE = energies[0];
	MinImage = 0;
	for (int l=1;l<=lastimg;l++) if (energies[l] < tempE) {
		tempE = energies[l];
		MinImage = l;
	}

	// Climbing Image
	tempE = energies[0];
	for (int l=1;l<=lastimg;l++) if(energies[l] > tempE) {
		tempE = energies[l];
		ClimbingImage = l;
	}

	// InterMinImages
	InterMinImages.clear();
	// intermediate minima at least e_thresh deep
	double relE;
	bool msdc;
	for (int l=1;l<lastimg;l++) {
		relE = energies[l]+e_thresh;
		if((relE<energies[l+1]) && (relE<energies[l-1])) {
			if(InterMinImages.size()==0) msdc = bool(initial.msd(systems[l-1],false)>MSD_THRESH);
			else msdc = bool(systems[*(std::prev(InterMinImages.end()))-1].msd(systems[l-1],false)>MSD_THRESH);
			if(msdc) InterMinImages.push_back(l);
		}
	}
};


};
#endif

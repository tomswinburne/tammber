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
	AbstractTaskMapper::insert("TASK_CARVE");
	AbstractTaskMapper::insert("TASK_SYMMETRY");
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

	MPI_Comm_rank(localComm_,&local_rank);
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
	BaseEngine::impls["TASK_FORCES"] = MDEngine::forces_impl; // insert("TASK_FORCES");
	BaseEngine::impls["TASK_CENTRO"] = MDEngine::centro_impl; // insert("TASK_CENTRO");
	BaseEngine::impls["TASK_CARVE"] = MDEngine::carve_impl; // insert("TASK_CARVE");
	BaseEngine::impls["TASK_SYMMETRY"] = MDEngine::symmetry_impl; // insert("TASK_SYMMETRY");
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
int local_rank;

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

std::function<void(GenericTask&)> symmetry_impl = [this](GenericTask &task) {
	/*
		clabels match by definition
		- find SelfSymmetries of targetInitial: SSI
		- find SelfSymmetries of targetFinal: SSF
		- find transformation targetInitial -> candidateInitial: TI
		- find transformation targetFinal -> candidateFinal: TF
	  - if TI*SSI == TF*SSF for some members of SSI,SSF, we have a match
		- if clabels match across transition, find transformation initial->final
	*/
	std::unordered_map<std::string,std::string> parameters=\
		extractParameters(task.type,task.flavor,defaultFlavor,taskParameters);
	//read parameters from the task
	bool self_check=safe_extractor<bool>(parameters,"SelfCheck",true);
	bool thresh_check=safe_extractor<bool>(parameters,"ThreshCheck",true);
	bool vf2_check=safe_extractor<bool>(parameters,"UseVF2",false);

	std::set<PointShiftSymmetry> self_symmetry, ops, f_ops;
	std::map<Label,std::set<PointShiftSymmetry>> self_symmetries;//,local_self_symmetries;
	std::map<int,int> c_map;
	Label c_lab,lab;
	bool cover;
	Transition transition_label;
	TransitionSymmetry trans;
	System carved;
	std::array<double,NDIM> shift;

	std::list<System> targets,candidates;

	extract("SelfSymmetries",task.arguments,self_symmetries);
	extract("Candidates",task.inputData,candidates);
	extract("Targets",task.inputData,targets);

	Cell bc = targets.begin()->getCell();

	// put everything in the basis of the first state? should be done beforehand??

	// first off, check for self symmetries (only required for first pair)

	if(self_check) {
		for(auto &target: targets) {
			c_lab = labeler->hash(target,true);
			if(self_symmetries.find(c_lab)==self_symmetries.end()) {
				LOGGER("Self isomorphism search for "<<c_lab)

				if(thresh_check) {
					if(vf2_check) {
						LOGGER("Using carved configuration + VF2 UseVF2")
						GenericTask t; t.clear();
						t.type=BaseEngine::mapper.type("TASK_CARVE");
						insert("State",t.inputData,target);
						BaseEngine::process(t);
						extract("ThreshState",t.outputData,carved);
						cover = SelfSymmetries(carved,self_symmetry);
					} else {
						LOGGER("Using carved configuration + brute id map")
						cover = SelfSymmetriesCarve(target,self_symmetry);
					}
				} else cover = SelfSymmetries(target,self_symmetry);
				LOGGER("Found "<<self_symmetry.size()<<" symmetries. Covering set: "<<cover)
				self_symmetries[c_lab] = self_symmetry;
				for(auto &ss : self_symmetry) LOGGER(ss.info_str())
			}
		}
	} else {
		// add identity just in case...
		PointShiftSymmetry null;
		self_symmetry.insert(null);
		if(self_symmetries.find(c_lab)==self_symmetries.end()) {
			self_symmetries[c_lab] = self_symmetry;
		} else {
			self_symmetries[c_lab].insert(null);
		}
		LOGGER("OMITTING self isomorphism searches, replacing with NULL")
	}

	insert("SelfSymmetries",task.returns,self_symmetries);

	if(targets.size()<2) return;

	auto itts = targets.begin();
	auto ittl = task.inputData.equal_range("Targets").first;

	LabelPair t_i_l = std::make_pair(ittl->second.canonical_label,ittl->second.label);
	LabelPair t_f_l = std::make_pair(std::next(ittl)->second.canonical_label,std::next(ittl)->second.label);

	LOGGER("Looking for isomorphisms to I->F == ("<<t_i_l.first<<","<<t_i_l.second<<") -> ("<<t_f_l.first<<","<<t_f_l.second<<")")


	// first, check for self transition
	double self_shift_mag=0.0,shift_mag=0.0;
	transition_label.first = t_i_l;
	transition_label.second = t_f_l;
	if (t_i_l.first == t_f_l.first) { // i.e. same canonical label
		LOGGER("Checking for self transform...")
		labeler->isomorphicMap(*itts,*std::next(itts),c_map); // 2 -> C -> 1
		ops = find_transforms(*itts,*std::next(itts),c_map);
		if(ops.size()>0) {
			auto _op = ops.begin();
			LOGGER("Found!\n"<<(*_op).info_str())
			trans = std::make_pair(transition_label,*_op);
			for(int ix=0;ix<NDIM;ix++) self_shift_mag += (_op->shift[ix]) * (_op->shift[ix]);
			self_shift_mag = sqrt(self_shift_mag);
			insert("SelfTransitions",task.returns,trans);
		}
	}



	if(candidates.size()==0) return;

	// now go through candidate list
	auto itcs = candidates.begin();
	auto itcl = task.inputData.equal_range("Candidates").first;
	int ccount=2;

	while (itcl!=task.inputData.equal_range("Candidates").second) {

		LabelPair c_i_l = std::make_pair(itcl->second.canonical_label,itcl->second.label); // candidate I
		LabelPair c_f_l = std::make_pair(std::next(itcl)->second.canonical_label,std::next(itcl)->second.label); // candidate F

		LOGGER("Target I->F == ("<<t_i_l.first<<","<<t_i_l.second<<") -> ("<<t_f_l.first<<","<<t_f_l.second<<")")
		LOGGER("Candidate I"<<ccount<<"->F"<<ccount<<" == ("<<c_i_l.first<<","<<c_i_l.second<<") -> ("<<c_f_l.first<<","<<c_f_l.second<<")")

		if ((c_i_l.first==t_i_l.first) && (c_f_l.first == t_f_l.first)) { // this is just a check, order should be done beforehand...


			double cand_shift_mag=0.0;
			if ((c_i_l.first==c_f_l.first) && (t_i_l.first==t_f_l.first)) {
				// find candidate self transform also
				LOGGER("Checking for candidate self transform...")
				labeler->isomorphicMap(*itcs,*std::next(itcs),c_map); // 2 -> C -> 1
				ops = find_transforms(*itcs,*std::next(itcs),c_map);
				if(ops.size()>0) {
					auto _op = ops.begin();
					LOGGER("Found!\n"<<_op->info_str())
					trans = std::make_pair(transition_label,*_op);
					for(int ix=0;ix<NDIM;ix++) cand_shift_mag += (_op->shift[ix]) * (_op->shift[ix]);
					cand_shift_mag = sqrt(cand_shift_mag);
					//insert("SelfTransitions",task.returns,trans);
				}
			}

			labeler->isomorphicMap(*itts,*itcs,c_map); // fills c_map : 2 -> C -> 1 for I
			ops = find_transforms(*itts,*itcs,c_map); // finds operations
			transition_label.first = t_i_l;
			transition_label.second = c_i_l;
			for(auto op:ops) { // should only be one....
				trans = std::make_pair(transition_label,op);
				insert("StateIsomorphisms",task.returns,trans);
			}

			labeler->isomorphicMap(*std::next(itts),*std::next(itcs),c_map); // 2 -> C -> 1 for F
			f_ops = find_transforms(*std::next(itts),*std::next(itcs),c_map);
			transition_label.first = t_f_l;
			transition_label.second = c_f_l;
			for(auto op:f_ops) { // should only be one....
				trans = std::make_pair(transition_label,op);
				insert("StateIsomorphisms",task.returns,trans);
			}

			for(auto op:ops) for(auto f_op:f_ops) {
				LOGGER("Map I->I"<<ccount<<":"<<op.info_str()) // MI,SI
			 	LOGGER("Map F->F"<<ccount<<":"<<f_op.info_str()) // MF,SF

				/*
				I2 = MI.I + SI = MI.iMIc.Ic + SI - MI.iMIc.Sc for MIc,Sc in SymI
				F2 = MF.F + SF = MF.iMFc.Ic + SF - MF.iMFc.Sc for MFc,Sc in SymF

				=> MI.iMIc = MF.iMFc, SI - MI.iMIc.Sc = SF - MF.iMFc.Sc
				*/

				// =?
				bool compatible = bool(op.operation==f_op.operation); // assume the same at this point...
				if(compatible) {
					for(int ii=0;ii<NDIM;ii++) shift[ii] = op.shift[ii]-f_op.shift[ii];
					bc.wrap(shift); // this should be the best wrap...
					shift_mag=0.0;
					for(int ii=0;ii<NDIM;ii++) shift_mag += shift[ii]*shift[ii];
					if(sqrt(shift_mag)>0.2) compatible=false; // quite small
				}

				// try to find I->I2 op from F
				PointShiftSymmetry i_sym_i,i_sym_f,cop_i,cop_f;

				/* CHECK TODO
				if(!compatible) for(auto sym_f:self_symmetries[t_f_l.first]) {
					i_sym_f = sym_f.inverse();
					cop_f = i_sym_f.compound(f_op);
					std::cout<<cop_f.info_str()<<"\n";
					std::cout<<"OR";
					cop_f = f_op.compound(i_sym_f);
					std::cout<<cop_f.info_str()<<"\n";
					std::cout<<"--------------------\n";
				}
				*/

				if (!compatible) for(auto sym_f:self_symmetries[t_f_l.first]) for(auto sym_i:self_symmetries[t_i_l.first])  {
					i_sym_i = sym_i.inverse();
					i_sym_f = sym_f.inverse();
					cop_i = i_sym_i.compound(op);
					cop_f = i_sym_f.compound(f_op);
					compatible = bool(cop_i.operation==cop_f.operation);
					/*
					if(compatible) {
						std::cout<<"cop_i==cop_f....";
						for(int ii=0;ii<NDIM;ii++) shift[ii] = cop_i.shift[ii]-cop_f.shift[ii];
						bc.wrap(shift);
						shift_mag=0.0;
						for(int ii=0;ii<NDIM;ii++) shift_mag += shift[ii]*shift[ii];
						if(sqrt(shift_mag)>0.2) {
							compatible=false;
							std::cout<<"but "<<sqrt(shift_mag)<<">0.2A!\n";
						} else std::cout<<"and it's good!\n";
					}
					*/
					if(compatible) break;

				}

				if (compatible) {
					LOGGER("Mapping compatible after accounting for self symmetries!\n")
					transition_label.first = c_i_l;
					transition_label.second = c_f_l;
					TransitionSymmetry trans = std::make_pair(transition_label,op);
					if(c_i_l.second==t_i_l.second) trans.second = f_op;
					insert("TransitionIsomorphisms",task.returns,trans);
					return;
				}
			}
			#ifdef VERBOSE
			if(ops.size()==0 || f_ops.size()==0) LOGGER("No mapping found\n")
			#endif
		} else {
			LOGGER("NOT ISOMORPHIC PAIR ORDERING FALSE")
		}
		itcl = std::next(itcl,2);
		itcs = std::next(itcs,2);
		ccount++;
	}
	return;
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

	double dX=0.0,dXmax=0.0;
	bool valid=false;

	/*
	Label reference_label,state_label;
	GenericTask t;
	t.arguments=task.arguments;
	t.type=BaseEngine::mapper.type("TASK_LABEL");
	t.flavor=task.flavor;

	t.clearInputs(); t.clearOutputs();
	insert("State",t.inputData,reference);
	BaseEngine::process(t);
	extract("Label",t.returns, reference_label);

	t.clearInputs(); t.clearOutputs();
	insert("State",t.inputData,state);
	BaseEngine::process(t);
	extract("Label",t.returns, state_label);

	if(reference_label!=state_label) {
	*/
		std::string filterType = parameters["Type"];
		if(task.arguments.find("Type")!=task.arguments.end())
			extract("Type",task.arguments,filterType);
		boost::trim(filterType);
		std::shared_ptr<AbstractTransitionFilter>
			filter = transitionFilterFactory.at(filterType)();
		filter->initialize(parameters);
		valid = filter->isValid(reference, state, parameters);
		dXmax = boost::lexical_cast<double>(parameters["dXmax"]);
		dX = boost::lexical_cast<double>(parameters["dX"]);
	//}

	LOGGER("MDEngine::FilterTransition: "<<valid<<" dX: "<<dX<<" dXmax: "<<dXmax)
	insert("Valid", task.returns, valid);
	insert("dX", task.returns, dX);
	insert("dXmax", task.returns, dXmax);
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


std::function<void(GenericTask&)> carve_impl = [this](GenericTask &task) {
	task.clearOutputs(); // just in case
	std::unordered_map<std::string,std::string> parameters=\
		extractParameters(task.type,task.flavor,defaultFlavor,taskParameters);
	//read parameters from the task
	int nn=safe_extractor<int>(parameters,"CentroNeighbors",8);
	double thresh=safe_extractor<double>(parameters,"Threshold",5.0);
	double scale=safe_extractor<double>(parameters,"RelativeCutoff",1.0);
	std::list<System> sysl;

	GenericTask centro;
	LOGGER("CSNN: "<<nn);

	centro.type = BaseEngine::mapper.type("TASK_CENTRO");
	insert("CentroNeighbors",centro.arguments,nn);
	extract("State",task.inputData,sysl);
	task.clearInputs();

	LOGGER("TADEngine::carve_impl : FOUND "<<sysl.size()<<" STATES")
	for(auto &s : sysl) {
		Cell bc = s.getCell();
		System thr_s;
		int clusters=1,thr_N=0;
		std::array<double,3> position={0.,0.,0.},fatom={0.,0.,0.},temp={0.,0.,0.};
		double temp_double=0.0;
		std::vector<double> csl;

		if(nn==0) {
			LOGGER("TADEngine::carve_impl : CentroNeighbors==0; assigning center of mass position")
			insert("Clusters",task.returns,clusters);
			insert("Position",task.returns,position);
			thr_N = s.getNAtoms();
			for(int i=0; i<thr_N;i++) for(int j=0;j<3;j++)
				position[j] += s.getPosition(i,j)/double(thr_N);
			thr_N = 0;
		} else {
			centro.clearInputs(); centro.clearOutputs();
			insert("State",centro.inputData,s);
			insert("CentroNeighbors",centro.arguments,nn);
			BaseEngine::process(centro);
			extract("CentroSymmetry",centro.outputData,csl);
			thr_s.setCell(s.getCell());
			for(int i=0; i<csl.size();i++) if(csl[i]>thresh) thr_N++;
			thr_s.setNAtoms(thr_N);
			if(thr_N==0) LOGGERA("TADEngine::carve_impl : No atoms above threshold!")
		}


		if(thr_N==0) {
			insert("ThreshState",task.outputData,s);
		} else {
			thr_N=0;
			LOGGER("TADEngine::carve_impl : Thresholding:")
			// Build new system from above-threshold atoms, scaled by RelativeCutoff
			for(int i=0; i<csl.size();i++) if(csl[i]>thresh) {
				LOGGER(s.getUniqueID(i)<<" "<<s.getPosition(i,0)<<" "<<s.getPosition(i,1)<<" "<<s.getPosition(i,2)<<" "<<csl[i]<<" "<<thresh)
				// new ID
				thr_s.setUniqueID(thr_N,s.getUniqueID(i));
				// snap to first position
				if(thr_N==0) for(int j=0;j<3;j++) fatom[j] = s.getPosition(i,j);
				for(int j=0;j<3;j++) temp[j] = s.getPosition(i,j)-fatom[j];
				bc.wrapc(temp);
				// new_x_j = pbc(x_j-x_i) + x_i
				for(int j=0;j<3;j++) thr_s.setPosition(thr_N,j,(fatom[j]+temp[j])/scale);
				thr_s.setSpecies(thr_N,s.getSpecies(i));
				thr_N++;
			}
			// Find clusters
			std::vector<int> cluster_occ;
			clusters = labeler->connectedComponents(thr_s,cluster_occ);
			// rescale positions
			for(int i=0;i<thr_N;i++) for(int j=0;j<3;j++) {
				temp_double = thr_s.getPosition(i,j) * scale;
				thr_s.setPosition(i,j,temp_double);
			}

			insert("ThreshState",task.outputData,thr_s);

			// find size and center of mass of each cluster
			std::vector<int> rocc(clusters,0); // size
			std::vector<std::array<double,3>> rocp(clusters,{0.,0.,0.}); // position
			for(int i=0;i<cluster_occ.size();i++) {
				for(int j=0;j<3;j++) rocp[cluster_occ[i]][j] += thr_s.getPosition(i,j);
				rocc[cluster_occ[i]]++;
			}

			// Return position of largest cluster
			int max_cl=0;
			for(int i=0;i<clusters;i++) {
				if(rocc[i]>0) for(int j=0;j<3;j++) rocp[i][j] /= float(rocc[i]);
				LOGGER("Cluster "<<i+1<<" : "<<rocc[i]<<" atoms, position :"<<rocp[i][0]<<" "<<rocp[i][1]<<" "<<rocp[i][2])
				if(rocc[i]>rocc[max_cl]) max_cl = i;
			}
			for(int j=0;j<3;j++) position[j] = rocp[max_cl][j];

		}
		LOGGER("TADEngine::carve_impl : NClusters = "<<clusters);
		LOGGER("TADEngine::carve_impl : Position = ["<<position[0]<<" "<<position[1]<<" "<<position[2]<<"]")
		insert("Position",task.returns,position);
		insert("Clusters",task.returns,clusters);
	}
};


/* HELPER FUNCTIONS */

// two is not passed by reference due to remap. finds mapping T,d such that T(y)+d = x  for  y,x == two,one
std::set<PointShiftSymmetry> find_transforms(System &one, System two, std::map<int,int> &map) {
	LOGGER("MDEngine::find_transforms");

	Cell bc = one.getCell();

	int natoms = one.getNAtoms();
	PointShiftSymmetry op;
	std::set<PointShiftSymmetry> ops;
	std::array<double,NDIM*NDIM> T;
	double min_d = 0.0, temp_d=0.0, a_mag=0.0;
	std::array<double,NDIM> temp,shift,cc={0.,0.,0.};
	for(int j=0;j<NDIM;j++) for(int k=0;k<NDIM;k++) cc[j] += bc.rsp[j][k]/2.0;

	two.remap(map);

	LOGGER("cell center: "<<cc[0]<<" "<<cc[1]<<" "<<cc[2])
	LOGGER("cell origin: "<<bc.origin[0]<<" "<<bc.origin[1]<<" "<<bc.origin[2])
	// Rotate around, then shift: X-c = G.(X-c) + d_i

	for(int operation=0; operation<48; operation++) {
		op.transform_matrix(T,operation); // fill transform matrix
		for(int j=0;j<NDIM;j++) {
			shift[j] = 0.0;
			for(int k=0;k<NDIM;k++)
				shift[j] += T[3*j+k]*(two.getPosition(0,k));//-two.getPosition(ca_two,k));//-cc[k]);
			shift[j] -= (one.getPosition(0,j));//-one.getPosition(ca_one,j));//-cc[j]); // shift = T(x1[0]) - x0[0]
		}

		bool complete=true;
		for(unsigned i=0; i<natoms; i++) {
			for(int j=0;j<NDIM;j++) {
				temp[j] = 0.;
				for(int k=0;k<NDIM;k++)
					temp[j] += T[NDIM*j+k] * (two.getPosition(i,k));//-two.getPosition(ca_two,k));//-cc[k]); // T(x1[i])
				temp[j] -= one.getPosition(i,j);//-one.getPosition(ca_one,j);//-cc[j];
				temp[j] -= shift[j]; // T(x1[i]) - x0[i] - shift = T(x1[i]-x1[0]) - (x0[i]-x0[0])
			}
			bc.wrapc(temp);
			for(int j=0;j<NDIM;j++) if(temp[j] * temp[j] > 0.05) {
				complete=false;
				break;
			}
		}
		if(complete) {
			bc.wrap(shift);
			LOGGER("       ["<<T[0]<<" "<<T[1]<<" "<<T[2]<<"]")
			LOGGER("matrix:["<<T[3]<<" "<<T[4]<<" "<<T[5]<<"]")
			LOGGER("       ["<<T[6]<<" "<<T[7]<<" "<<T[8]<<"]")
			LOGGER("shift: ["<<shift[0]<<" "<<shift[1]<<" "<<shift[2]<<"]")
			op.operation=operation;
			op.matrix=T;
			for(int j=0;j<NDIM;j++) op.shift[j] = shift[j];
			op.valid = true;
			ops.insert(op);
		}
	}
	return ops;
};

bool SelfSymmetries(System &one, std::set<PointShiftSymmetry> &syms,bool return_maps=false) {
	LOGGER("MDEngine::SelfSymmetries");
	std::list<std::map<int,int>> amaps;
	syms.clear();
	labeler->isomorphicSelfMaps(one,amaps);
	LOGGER("FOUND "<<amaps.size()<<" SELF MAPS");
	int count = 0;
	for (auto &map: amaps) {
		auto ops = find_transforms(one, one, map);
		for(auto op: ops) {
			if(return_maps) op.map = map;
			syms.insert(op);
		}
		count++;
		if(count>48) {
			LOGGER("TRIED 96 SELF MAPS, FOUND "<<syms.size()<<" OPERATIONS. EXITING")
			break;
		}
	}
	return bool(amaps.size()<=syms.size());
};


// VF2-free implementation

bool find_transforms_direct(System &one, System two, PointShiftSymmetry &op, std::array<double,NDIM> position, int operation) {
	LOGGER("MDEngine::find_transforms_direct");

	Cell bc = one.getCell();
	int natoms = one.getNAtoms();
	unsigned min_i2;

	std::array<double,NDIM*NDIM> T;
	double min_d = 0.0, temp_d=0.0, a_mag=0.0;
	std::array<double,NDIM> temp,cc={0.,0.,0.},cc2,temp2;
	for(int j=0;j<NDIM;j++) for(int k=0;k<NDIM;k++) cc[j] += bc.rsp[j][k]/2.0;

	std::vector<int> count(natoms,0);

	LOGGER("cell center: "<<cc[0]<<" "<<cc[1]<<" "<<cc[2])
	LOGGER("cell origin: "<<bc.origin[0]<<" "<<bc.origin[1]<<" "<<bc.origin[2])

	// Rotate around position with transform matrix
	// one : centered_cluster + com
	// two : T.centered_cluster (no com shift)
	op.transform_matrix(T,operation);
	for(int j=0;j<NDIM;j++) cc[j]=0.0;
	for(int j=0;j<NDIM;j++) cc2[j]=0.0;
	for(unsigned i=0; i<natoms; i++) {
		for(int j=0;j<NDIM;j++) {
			temp[j] = 0.;
			cc[j] += ( one.getPosition(i,j) - position[j] ) / natoms;
			for(int k=0;k<NDIM;k++)
				temp[j] += T[NDIM*j+k] * (one.getPosition(i,k)-position[k]);
			two.setPosition(i,j,temp[j]);
			cc2[j] += two.getPosition(i,j) / natoms;
		}
	}
	LOGGER("one com - position: "<<cc[0]<<" "<<cc[1]<<" "<<cc[2])
	LOGGER("two com - position: "<<cc2[0]<<" "<<cc2[1]<<" "<<cc2[2])

	temp_d=0.; for(int j=0;j<NDIM;j++) temp_d += (cc[j]-cc2[j])*(cc[j]-cc2[j]);
	LOGGER("com difference: "<<sqrt(temp_d))
	if(sqrt(temp_d)>1.0) return false; // quit if not identical here


	// O(N^2) closest remapping- if a symmetry op the com should be unchanged...
	for(unsigned i1=0; i1<natoms; i1++) {
		for(int j=0;j<NDIM;j++) temp[j] = one.getPosition(i1,j)-position[j];
		min_d = 100000000.0;
		min_i2 = i1;
		for(unsigned i2=0; i2<natoms; i2++) {
			if(one.getSpecies(i1)!=two.getSpecies(i2)) continue;
			temp_d = 0.0;
			for(int j=0;j<NDIM;j++)
				temp_d += (two.getPosition(i2,j)-temp[j])*(two.getPosition(i2,j)-temp[j]);
			if(temp_d<min_d) {
				min_i2 = i2;
				min_d = temp_d;
			}
		}
		if(min_d>0.05) break;
		else count[min_i2] += 1;
	}
	for(unsigned i1=0; i1<natoms; i1++) if(count[i1]!=1) return false;


	/*
		Match found!
		r_two = two + com
		one-com = T.(r_two - com)
		one + shift = T.r_two
		shift = T.r_two - one = T.com - com
	*/
	for(int j=0;j<NDIM;j++) {
		temp[j] = 0.;
		for(int k=0;k<NDIM;k++) temp[j] += T[NDIM*j+k] * position[k];
		temp[j] -= position[j];
	}
	bc.wrap(temp);
	LOGGER("       ["<<T[0]<<" "<<T[1]<<" "<<T[2]<<"]")
	LOGGER("matrix:["<<T[3]<<" "<<T[4]<<" "<<T[5]<<"]")
	LOGGER("       ["<<T[6]<<" "<<T[7]<<" "<<T[8]<<"]")
	LOGGER("shift: ["<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<"]")
	op.operation=operation;
	op.matrix=T;
	for(int j=0;j<NDIM;j++) op.shift[j] = temp[j];
	op.valid = true;
	return true;
};




bool SelfSymmetriesCarve(System &one, std::set<PointShiftSymmetry> &syms) {
	LOGGER("MDEngine::SelfSymmetriesCarve");
	GenericTask t;
	System carved;
	int clusters=1;
	std::array<double,NDIM> position;
	bool found_transform;
	PointShiftSymmetry op;

	t.clear();
	t.type=BaseEngine::mapper.type("TASK_CARVE");
	insert("State",t.inputData,one);
	BaseEngine::process(t);

	extract("ThreshState",t.outputData,carved);
	extract("Clusters",t.returns,clusters);
	extract("Position",t.returns,position);

	if(clusters>1) return false;

	// Rotate around, then shift: X-c = G.(X-c) + d_i
	for(int operation=0; operation<48; operation++)
		if(find_transforms_direct(carved,carved,op,position,operation))
			syms.insert(op);
	LOGGER("TRIED 48 SELF MAPS, FOUND "<<syms.size()<<" OPERATIONS. EXITING")
	return true;
};

};

#endif

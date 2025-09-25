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

#include "TammberModel.hpp"

SymmLabelPair::SymmLabelPair(){
	first=0;
	second=0;
};

SymmLabelPair::SymmLabelPair(LabelPair other) {
	first=other.first;
	second=other.second;
};

SymmLabelPair SymmLabelPair::rev() {
	SymmLabelPair res;
	res.second = first;
	res.first = second;
	return res;
};

LabelPair SymmLabelPair::ns() {
	LabelPair res;
	res.first = first;
	res.second = second;
	return res;
};


SymmLabelPair CanonTrans(Transition trlab) {
 SymmLabelPair res(std::make_pair(trlab.first.first,trlab.second.first));
 return res;
};

SymmLabelPair NonCanonTrans(Transition trlab) {
	SymmLabelPair res(std::make_pair(trlab.first.second,trlab.second.second));
	return res;
};

// for submission
NEBjob::NEBjob(){
	LOGGER("NEBjob::NEBjob")
	/* nothing */
};
NEBjob::NEBjob(Transition tt) {
	LOGGER("NEBjob::NEBjob")
	TargetTransition = tt;
	ExistingPairs.clear();
	InitialSymmetries.clear();
	FinalSymmetries.clear();
};
NEBjob::NEBjob(Transition tt,std::list<LabelPair> epl) {
	LOGGER("NEBjob::NEBjob")
	TargetTransition = tt;
	ExistingPairs.clear();
	ExistingPairs = epl;
	InitialSymmetries.clear();
	FinalSymmetries.clear();
};

TADjob::TADjob(){
	LOGGER("TADjob::TADjob")
	nInstances=0;
};
TADjob::TADjob(LabelPair il,double t,int n) {
	LOGGER("TADjob::TADjob")
	InitialLabels = il;
	temperature = t;
	nInstances = n;
};
TADjob::TADjob(LabelPair il,double t,int n,std::map<LabelPair,double> bl) {
	LOGGER("TADjob::TADjob")
	InitialLabels = il;
	temperature = t;
	nInstances = n;
	BasinLabels = bl;
};

// all transition/saddle information here
Connection::Connection(){
	LOGGER("Connection::Connection")
	saddleE = 10.0 * MAX_BARRIER;
	dX = 2.0 * MSD_THRESH;
	dXmax = 2.0 * MSD_THRESH;
	Ftol = 10.0;
	fpt = std::make_pair(0.0,0.0);
	fpT = std::make_pair(0.0,0.0);
	nu = std::make_pair(PRIOR_NU,PRIOR_NU);
	priornu = std::make_pair(PRIOR_NU,PRIOR_NU);
};
Connection::Connection(NEBPathway &path) {
	LOGGER("Connection::Connection")
	Ftol = 10.0;
	add_path(path);
};

void Connection::add_path(NEBPathway &path) {
	LOGGER("Connection::add_path")
	if(path.Ftol < Ftol) {
		saddleE = path.saddleE;
		dX = path.dX;
		dXmax = path.dXmax;
		Ftol = path.Ftol;
		energies = path.energies;
		SaddleLabels = path.SaddleLabels;
		priornu = path.priornu;
	}
	for(auto &nst:path.self_transitions) {
		bool match=false; // why is self_transitions a std::list ??
		for(auto &st:self_transitions) if(st.first==nst.first) match=true;
		if(!match) self_transitions.push_back(nst);
	}
};



void Connection::add_jump(double temperature, double _fpt, bool forwards) {
	LOGGER("Connection::add_jump")
	if(forwards) {
		if(fpt.first==0.0) {
			fpt.first = _fpt;
			fpT.first = temperature;
		}
		counts.first += 1.0; // another count
	} else {
		if(fpt.second==0.0) {
			fpt.second = _fpt;
			fpT.second = temperature;
		}
		counts.second += 1.0; // another count
	}
	LOGGER("END Connection::add_jump")
};

void Connection::assimilate(Connection& other) {
	LOGGER("Connection::assimilate")
	if(saddleE>other.saddleE) {
		saddleE=other.saddleE;
		SaddleLabels=other.SaddleLabels;
		dX=other.dX;
		dXmax=other.dXmax;
		Ftol=other.Ftol;
		energies=other.energies;
		priornu = other.priornu;
		nu = other.nu;
	}

	counts.first += other.counts.first;
	counts.second += other.counts.second;
	/*for(auto count: other.counts) {
		if(counts.find(count.first)==counts.end())
			counts.insert(std::make_pair(count.first,std::make_pair(0.0,0.0)));
		counts[count.first].first += count.second.first;
		counts[count.first].second += count.second.second;
	}*/

	if(other.fpt.first>0.0 && other.fpT.first<=fpT.first && other.fpt.first<fpt.first) {
		fpt.first = other.fpt.first;
		fpT.first = other.fpT.first;
	}

	if(other.fpt.second>0.0 && other.fpT.second<=fpT.second && other.fpt.second<fpt.second) {
		fpt.second = other.fpt.second;
		fpT.second = other.fpT.second;
	}

};

// Canonical labels implicit?
StateEdge::StateEdge(){
	LOGGER("StateEdge::StateEdge()")
	/* do nothing */
};

StateEdge::StateEdge(SymmLabelPair labels_) {
	LOGGER("StateEdge::StateEdge")
	labels = labels_;
};

StateEdge::StateEdge(SymmLabelPair labels_, SymmLabelPair refc_) {
	LOGGER("StateEdge::StateEdge")
	labels = labels_;
	PointShiftSymmetry op; op.valid=true;
	self_edge_map.insert(std::make_pair(refc_,std::make_pair(refc_,op)));
};

std::string StateEdge::info_str(bool seg) {
	LOGGER("StateEdge::info_str")
	std::string res="\n Edge: "+std::to_string(labels.first);
	if(labels.first==labels.second) res+=" (Self Transition)\n";
	else res+=" <-> "+std::to_string(labels.second)+"\n";
	res += "\n   Connections:\n";
	for(auto &conn: connections) {
		if (conn.second.dX < MSD_THRESH) continue;
		res += "  "+std::to_string(conn.first.first)+" -> "+std::to_string(conn.first.second)+": ";
		double initialE = conn.second.energies[0];
		res += std::to_string(conn.second.saddleE-initialE)+",";
		double finalE = *std::prev(conn.second.energies.end());
		res += std::to_string(conn.second.saddleE-finalE)+" ";
		res += "dX Ftol : "+std::to_string(conn.second.dX)+" "+std::to_string(conn.second.Ftol)+"\n";
		if(labels.first==labels.second) res += "  counts: "+std::to_string(conn.second.counts.first+conn.second.counts.second)+"\n";
		else res += "  counts: "+std::to_string(conn.second.counts.first)+" "+std::to_string(conn.second.counts.second)+"\n";

		if(seg) {
			double baseE = initialE;
			for(auto eee: conn.second.energies) if(eee<baseE) baseE=eee;
			res += "     E[]: "+std::to_string(baseE)+" + \n";
			for(auto eee: conn.second.energies) {
				res += "    ";
				for(int ss=0;ss<int(30.0*(eee-baseE)/(conn.second.saddleE-baseE));ss++) res+=" ";
				res += "| "+std::to_string(eee-baseE)+"\n";
			}
			if(conn.second.self_transitions.size()>0) {
				res += "\n Self Symmetries:\n";
				for(auto st: conn.second.self_transitions) res += st.second.info_str();
				res+="\n";
			}
		}
	}
	if(seg && self_edge_map.size()>0)
		res += "\n Seen Equivalent Transitions:"+std::to_string(self_edge_map.size())+"\n";
	res+="\n   Completed NEBS: ";
	for(auto &pn: completedNEBS) res += std::to_string(pn.first)+","+std::to_string(pn.second)+" ";
	res+="\n   Requested NEBS: ";
	for(auto &pn: requestedNEBS) res += std::to_string(pn.first)+","+std::to_string(pn.second)+" ";
	res+="\n   Pending NEBS: ";
	for(auto &pn: pendingNEBS) res += std::to_string(pn.first)+","+std::to_string(pn.second)+" ";
	res+="\n   Failed NEBS: "; // TO BE RENAMED
	for(auto &pn: requestedPMS) res += std::to_string(pn.first)+","+std::to_string(pn.second)+" ";
	res += "\n";
	return res;
};

// This is similar to TADStateStatistics- one for each canonical label
StateVertex::StateVertex(){
	LOGGER("StateVertex::StateVertex()")
	/* nothing */
};
StateVertex::StateVertex(LabelPair ref,double energy_,int id_){
	LOGGER("StateVertex::StateVertex")
	reference_label = ref;
	energy=energy_;
	target_state_time = 1.0;
	duration = 0;
	overhead = 0;
	trials = 0;
	id=id_;
	position = {0.,0.,0.};
	clusters=0;
	elapsedTime.clear();
	self_symmetries.clear();
};

void StateVertex::add_segment(TADSegment &seg) {
	LOGGER("StateVertex::add_segment")
	duration += seg.duration;
	int tint = int(seg.temperature+0.5);
	if( elapsedTime.find(tint)==elapsedTime.end() )
		elapsedTime.insert( std::make_pair(tint,0.0) );
	elapsedTime[tint] += seg.elapsedTime;
	if(clusters==0) {
		clusters=seg.initialClusters;
		position = seg.initialPosition;
	}
};

// recalculate effective state time
double StateVertex::state_time(double targetT) {
	LOGGER("StateVertex::state_time")
	double iT=0.0,time=1.0;
	for (auto rit=elapsedTime.rbegin(); rit!=elapsedTime.rend(); ++rit) {
		if(rit!=elapsedTime.rbegin())
			time *= exp( (log(std::max(time,1.0))+LOG_NU_MIN) * (1.0-rit->first*iT) ); // to new temperature
		time += rit->second; // add time
		iT = 1.0/rit->first; // new iT
	}
	time *= exp( (log(std::max(time,1.0))+LOG_NU_MIN) * (1.0-targetT*iT) );
	if(std::isnan(time)) time=1.0;
	return time;
};

std::string StateVertex::info_str(double targetT,double eshift=0.0,double ku=0.0) {
	LOGGER("StateVertex::info_str")
	double time = state_time(targetT);
	std::string res=" "+std::to_string(reference_label.first)+" ("+std::to_string(reference_label.second)+")";
	double dE = energy-eshift;
	res += " Clusters: "+std::to_string(clusters)+", Position: [";
	res += std::to_string(position[0])+",";
	res += std::to_string(position[1])+",";
	res += std::to_string(position[2])+"],";
	res += " E-Emin:"+std::to_string(dE)+"eV, t:"+std::to_string(time)+"ps, ";
	res += "ku:"+std::to_string(ku)+"THz blocks:"+std::to_string(duration)+", overhead:"+std::to_string(overhead)+"\n";
	return res;
};

void StateVertex::update(Label lab, double energy_,bool force) {
	LOGGER("StateVertex::update")
	if(energy>energy_ or force) energy=energy_;
	//reference_label.second = lab;
};

void StateVertex::update(Label lab, double energy_,int clust, std::array<double,3> pos,bool force) {
	LOGGER("StateVertex::update : UPDATING STATE")
	update(lab,energy_,force);
	if((reference_label.second==lab && clust>0) or force) {
			clusters = clust;
			position = pos;
	}
};


void StateVertex::update(Label lab, double energy_,int clust, std::array<double,3> pos, std::set<PointShiftSymmetry> ss, bool force) {
	LOGGER("StateVertex::update : UPDATING STATE")
	update(lab,energy_,clust,pos,force);
	for(auto s:ss) self_symmetries.insert(s);
};

void Rate::reset(size_t tadTsize) {
	LOGGER("Rate::reset")
	std::pair<double,double> pp = std::make_pair(0.0,-1.0);
	target_k_fp = pp;
	tad_k_fp.clear();
	for(size_t i=0;i<tadTsize;i++) tad_k_fp.push_back(pp);
};

TammberModel::TammberModel(){
	LOGGER("TammberModel::TammberModel()")
	max_id=0;
};

TammberModel::TammberModel(double targetT_, double minT_, double maxT_, int nT_, int pfCT, double mB) {
	LOGGER("TammberModel::TammberModel")
	parametrize(targetT_, minT_, maxT_, nT_, pfCT, mB);
};

void TammberModel::parametrize(double targetT_, double minT_, double maxT_, int nT_, int pfCT, double mB) {
	LOGGER("TammberModel::parametrize")
	max_id=0;
	// TAD parameters. Only inverse temperature counts...
	targetB = 1.0/targetT_/BOLTZ;
	targetT = targetT_;
	tadT.clear();
	if(nT_>1) {
		double dT=(maxT_-minT_)/double(std::max(1,nT_-1));
		for(int i=0; i<nT_; i++) tadT.push_back(minT_+dT*i);
	} else tadT.push_back(minT_);
	// Bayesian parameters / priors
	max_k = 10.0; // Max rate in THz
	min_ku = 1.0e-20; // Min ku in THz - to force wider sampling
	max_ku = 0.1; // Max ku in THz - effective prior over distribution
	prefactorCountThresh = double(pfCT);
	min_barrier = mB;
	HashCost = 250.0;
	NEBCost = 1000.0;
	RhoInitFlavor = 0;
	PredictionSize = 3;
};

void TammberModel::initialize(boost::property_tree::ptree &config,bool restart){
	LOGGER("TammberModel::initialize")
	if(!restart) max_id=0;

	targetT = config.get<double>("Configuration.MarkovModel.TargetTemperature",300.0);

	if(restart) for(auto &v : StateVertices)
		v.second.target_state_time = v.second.state_time(targetT);

	targetB = 1.0/targetT/BOLTZ;
	int nT = config.get<int>("Configuration.MarkovModel.TemperatureSteps",11);
	double minT = config.get<double>("Configuration.MarkovModel.MinTemperature",300.0);
	double maxT = config.get<double>("Configuration.MarkovModel.MaxTemperature",600.0);
	HashCost = config.get<double>("Configuration.MarkovModel.HashCost",250.0);
	NEBCost = config.get<double>("Configuration.MarkovModel.NEBCost",1000.0);
	RhoInitFlavor = config.get<int>("Configuration.MarkovModel.RhoInitFlavor",0);
	AllocScheme = config.get<int>("Configuration.MarkovModel.AllocScheme",0);

	ClusterThresh = config.get<int>("Configuration.MarkovModel.ClusterThresh",0);

	// Need at least a 20% chance of dephasing otherwise we suppress
	DephaseThresh = config.get<double>("Configuration.MarkovModel.DephaseThresh",0.2);

	safe_opt = config.get<bool>("Configuration.MarkovModel.SafeOpt",true);
	PredictionSize = config.get<int>("Configuration.MarkovModel.PredictionSize",10);
	prefactorCountThresh = config.get<double>("Configuration.MarkovModel.PrefactorCountThresh",2.0);

	// Only if specified in input file. Experimental feature
	sim_conn = config.get<bool>("Configuration.MarkovModel.EstimatePendingNEBS",false);
	pNEB_Prior = config.get<double>("Configuration.MarkovModel.PendingNEBSPrior",1.0);

	// include "dead" saddles, produced by erroneous high temperature minimization
	include_shallow_states = config.get<bool>("Configuration.MarkovModel.IncludeShallowStates",false);

	// Only for postprocessing
	LatticeConstant = config.get<std::string>("Configuration.MarkovModel.Lattice","None");
	PrimitiveUnitCell = config.get<std::string>("Configuration.MarkovModel.PrimitiveUnitCell","None");
	UnitCell = config.get<std::string>("Configuration.MarkovModel.UnitCell","None");
	SuperCell = config.get<std::string>("Configuration.MarkovModel.SuperCell","None");

	if(ClusterThresh<=0) ClusterThresh = 100000;

	double dT = (maxT-minT) / double(std::max(1,nT-1));
	tadT.clear();
	for(int i=0;i<nT;i++) tadT.push_back(minT+dT*i);

	// Stitching parameter. Basins should have escape of at least 5ps at minT
	min_barrier = log(5.0) * BOLTZ * minT; //  about 0.04eV at 300K
	max_k = 10.0; // Max rate in THz
	min_ku = 1.0e-20; // Min ku in THz - to force wider sampling when one state dominant
	max_ku = 0.1; // Max ku in THz - effective 'low pass' prior over distribution
	if(restart) for(auto &e: StateEdges) e.second.requestedNEBS.clear();
};

int TammberModel::newIndex() {
	LOGGER("TammberModel::newIndex")
	return max_id++;
};

// to be overloaded
void TammberModel::add_vertex(LabelPair labels, double energy,bool force) {
	LOGGER("TammberModel::add_vertex")
	if(StateVertices.size()==0) InitialLabels = labels;
	auto svp = StateVertices.find(labels.first);
	if(svp==StateVertices.end()) {
		int id = newIndex(); // 0 index, should be unique even if states are deleted
		LOGGER("ADDING VERTEX: "<<labels.first<<","<<labels.second<<" E:"<<energy<<"eV")
		StateVertices.insert(std::make_pair(labels.first,StateVertex(labels,energy,id)));
	} else svp->second.update(labels.second,energy,force);
};

void TammberModel::add_vertex(LabelPair labels, double energy, int clusters,std::array<double,3> pos, bool force) {
	LOGGER("TammberModel::add_vertex")
	if(StateVertices.size()==0) InitialLabels = labels;
	auto svp = StateVertices.find(labels.first);
	if(svp==StateVertices.end()) {
		int id = newIndex(); // 0 index, should be unique even if states are deleted
		LOGGER("ADDING VERTEX: "<<labels.first<<","<<labels.second<<" E:"<<energy<<"eV "<<"Nclusters: "<<clusters<<" pos:["<<pos[0]<<","<<pos[1]<<","<<pos[2]<<"]")
		StateVertices.insert(std::make_pair(labels.first,StateVertex(labels,energy,id)));
	}
	svp = StateVertices.find(labels.first);
	svp->second.update(labels.second,energy,clusters,pos,force);
};

void TammberModel::add_vertex(LabelPair labels, double energy,int clusters,std::array<double,3> pos, std::set<PointShiftSymmetry> ss,bool force) {
	LOGGER("TammberModel::add_vertex")
	if(StateVertices.size()==0) InitialLabels = labels;
	auto svp = StateVertices.find(labels.first);
	if(svp==StateVertices.end()) {
		int id = newIndex(); // 0 index, should be unique even if states are deleted
		LOGGER("ADDING VERTEX: "<<labels.first<<","<<labels.second<<" E:"<<energy<<"eV "<<"Nclusters: "<<clusters<<" pos:["<<pos[0]<<","<<pos[1]<<","<<pos[2]<<"]")
		StateVertices.insert(std::make_pair(labels.first,StateVertex(labels,energy,id)));
	}
	svp = StateVertices.find(labels.first);
	svp->second.update(labels.second,energy,clusters,pos,ss,force);
};

// bare minimum edge
void TammberModel::add_edge(SymmLabelPair el) {
	LOGGER("TammberModel::add_edge")
	if(StateEdges.find(el)==StateEdges.end()){
		StateEdge nse(el);
		StateEdges.insert(std::make_pair(el,nse));
	}
};

// add refconnection as well bare minimum edge
void TammberModel::add_edge(SymmLabelPair el,SymmLabelPair tl) {
	LOGGER("TammberModel::add_edge")
	if(StateEdges.find(el)==StateEdges.end()){
		StateEdge nse(el);
		StateEdges.insert(std::make_pair(el,nse));
	}

	bool forwards = bool(el.first==StateEdges.find(el)->first.first);

	auto semp = &(StateEdges.find(el)->second.self_edge_map);
	PointShiftSymmetry op; op.valid=true;// identity operation
	if(semp->find(tl)==semp->end()) {
		if(forwards) semp->insert(std::make_pair(tl,std::make_pair(tl,op)));
		else semp->insert(std::make_pair(tl.rev(),std::make_pair(tl.rev(),op)));
	}

};

void TammberModel::add_edge(Transition tt) {
	LOGGER("TammberModel::add_edge")
	add_edge(CanonTrans(tt),NonCanonTrans(tt));
};

// could add some more here
bool TammberModel::fault_check(TADSegment &seg) {
	LOGGER("TammberModel::fault_check")
	bool res = false;
	//if(seg.InitialLabels.first==0||seg.InitialLabels.second==0) res=true;
	//if(seg.transition.first.first==0 || seg.transition.first.second==0) res=true;
	//if(seg.transition.second.first==0 || seg.transition.second.second==0) res=true;
	if(!seg.dephased) res=true;
	if(seg.trials>0 && StateVertices.find(seg.InitialLabels.first)!=StateVertices.end()) {
		auto vertex = &(StateVertices.find(seg.InitialLabels.first)->second);
		vertex->trials += seg.trials;
		vertex->overhead += seg.overhead;
	};

	if(res) LOGGER("FAULTY/UNDEPHASED SEGMENT "<<seg.info_str())
	return res;
};

// from TAD
void TammberModel::add_segment(TADSegment &seg) {
	LOGGER("TammberModel::add_segment")
	// Is it valid? (dephased, sensible labels etc)
	if(fault_check(seg)) return;

	LOGGER("ADDING SEGMENT "<<seg.info_str())
	// Add vertices if required- currently treating transition labels as independent (see *)
	//add_vertex(seg.InitialLabels,seg.initialE);
	add_vertex(seg.InitialLabels,seg.initialE,seg.initialClusters,seg.initialPosition);


	if(seg.transition.first!=seg.InitialLabels) {
		//if(StateVertices.find(seg.transition.first.first)==StateVertices.end());
		// BASIN NOT YET IMPLEMENTED...
		LOGGER(seg.transition.first.first<<","<<seg.transition.first.second<<" => "<<seg.InitialLabels.first<<" => "<<seg.InitialLabels.second)
		//add_vertex(seg.transition.first,seg.BasinLabels[seg.transition.first]);
	}

	// basin labels are just for monitoring at the moment
	auto vertex = &(StateVertices.find(seg.InitialLabels.first)->second);
	for(auto bl: seg.BasinLabels) vertex->BasinLabels.insert(bl);

	// (*) every state should have its own vertex which I later collate....
	vertex = &(StateVertices.find(seg.transition.first.first)->second);
	vertex->add_segment(seg); // just for state time etc

	if (seg.OutOfBasin) { // we have a transition...
		//add_vertex(seg.transition.second,seg.finalE);
		add_vertex(seg.transition.second,seg.finalE,seg.finalClusters,seg.finalPosition);

		double fpt = vertex->state_time(seg.temperature);

		if(transitionMap.find(seg.transition)==transitionMap.end()) {
			transitionMap.insert(std::make_pair(seg.transition,seg.transition));
			transitionMap.insert(std::make_pair(seg.transition.rev(),seg.transition.rev()));
			LOGGER("TammberModel::add_segment transitionMap adding map")
		}

		Transition trans =  transitionMap.at(seg.transition);
		LOGGER("TammberModel::add_segment transitionMap mapping: ")
		LOGGER(seg.transition.first.first<<","<<seg.transition.first.second<<" -> "<<seg.transition.second.first<<","<<seg.transition.second.second)
		LOGGER("\t\t|\t\t|\t\t|\n\t\tV\t\tV\t\tV")
		LOGGER(trans.first.first<<","<<trans.first.second<<" -> "<<trans.second.first<<","<<trans.second.second)

		add_edge(trans); // find / create edge (also takes care of initializing self_edge_map)

		bool forwards=bool(CanonTrans(trans).first==StateEdges.find(CanonTrans(trans))->first.first);

		auto ep = &(StateEdges.find(CanonTrans(trans))->second); // pointer to edge

		SymmLabelPair tl = NonCanonTrans(trans); // not nec. parallel

		// No NEB? add to pendingNEBS unless previously failed
		if(ep->completedNEBS.find(tl)==ep->completedNEBS.end())
			if(ep->requestedPMS.find(tl)==ep->requestedPMS.end()) { // TO BE RENAMED!
				if(forwards) ep->pendingNEBS.insert(tl);
				else ep->pendingNEBS.insert(tl.rev());
				LOGGER("TammberModel::add_segment APPENDING TO pendingNEBS")
		}

		// what is the "canonical" noncanonical transition label?
		SymmLabelPair ctl = ep->self_edge_map.at(tl).first; // parallel

		// No connection yet? add to pendingTADS
		if(ep->connections.find(ctl)==ep->connections.end()) {
			if(pendingTADS.find(trans)==pendingTADS.end()) {
				std::list<TADSegment> ttsl;
				ttsl.push_back(seg);
				pendingTADS.insert(std::make_pair(trans,ttsl));
			} else pendingTADS[trans].push_back(seg);
			LOGGER("TammberModel::add_segment APPENDING TO pendingTADS")
		} else { // otherwise we have a normal transition
			ep->connections.find(ctl)->second.add_jump(seg.temperature,fpt,forwards);
		}
	} // OutOfBasin
};

void TammberModel::add_symmetries(NEBPathway &path) {
	LOGGER("TammberModel::add_symmetries")
	// first off, incorporate the symmetry results
	add_vertex(path.InitialLabels,path.initialE,true);
	add_vertex(path.FinalLabels,path.finalE,true);

	// find verticies
	auto initial_vertex = &(StateVertices.find(path.InitialLabels.first)->second);
	auto final_vertex = &(StateVertices.find(path.FinalLabels.first)->second);

	// add self_symmetries;
	for(auto ss: path.InitialSymmetries) initial_vertex->self_symmetries.insert(ss);
	for(auto ss: path.FinalSymmetries) final_vertex->self_symmetries.insert(ss);

	// add self mappings IF from the reference label (currently just as we have it, will be used...)
	for(auto &es: path.equivalent_states) { // ( (LabelPair,LabelPair), PointShiftSymmetry )
		// add mapping if from a reference_label
		if(es.first.first == initial_vertex->reference_label) // i.e. matching LabelPair
			initial_vertex->state_isomorphisms.insert(std::make_pair(es.first.second.second,es.second));
		else if(es.first.first == final_vertex->reference_label)
				final_vertex->state_isomorphisms.insert(std::make_pair(es.first.second.second,es.second));

		// try reverse...
		if(es.first.second == initial_vertex->reference_label) // i.e. matching LabelPair
			initial_vertex->state_isomorphisms.insert(std::make_pair(es.first.first.second,es.second.inverse()));
		else if(es.first.second == final_vertex->reference_label)
				final_vertex->state_isomorphisms.insert(std::make_pair(es.first.first.second,es.second.inverse()));
	}
	// Also add maps if a self transition
	for(auto st: path.self_transitions) { // (LabelPair,LabelPair), PointShiftSymmetry
		if (st.first.first == initial_vertex->reference_label)
			initial_vertex->state_isomorphisms.insert(std::make_pair(st.first.second.second,st.second));
		if (st.first.first == final_vertex->reference_label)
			final_vertex->state_isomorphisms.insert(std::make_pair(st.first.second.second,st.second));
	}
};

bool TammberModel::add_transitionMaps(NEBPathway &path) {
	LOGGER("TammberModel::add_transitionMaps")
	if (path.FoundTransitions.size()==0) return false;
	LOGGER("TammberModel::add_transitionMaps segment intial -> final :")
	LOGGER(path.InitialLabels.first<<","<<path.InitialLabels.second<<") -> "<<path.FinalLabels.first<<","<<path.FinalLabels.second)
	LOGGER("TammberModel::add_transitionMaps Going through FoundTransitions:")

	for(auto ft:path.FoundTransitions) {
		LOGGER(ft.print_str())
		add_edge(ft); // find or create edge
		auto ep = StateEdges.find(CanonTrans(ft));

		// add to pendingNEBS, checking for direction relative to StateEdge
		auto nctp = NonCanonTrans(ft);
		if(CanonTrans(ft).first!=ep->first.first) nctp = NonCanonTrans(ft).rev(); // now parallel
		if(ep->second.requestedNEBS.find(nctp)==ep->second.requestedNEBS.end()) // not requested
			if(ep->second.completedNEBS.find(nctp)==ep->second.completedNEBS.end()) // not completed
				if(ep->second.pendingNEBS.find(nctp)==ep->second.pendingNEBS.end()) // not pending
					if(ep->second.requestedPMS.find(nctp)==ep->second.requestedPMS.end()) // not failed TO BE RENAMED!!
						ep->second.pendingNEBS.insert(nctp);

		// add to transitionMap
		transitionMap[ft] = ft;
		transitionMap[ft.rev()] = ft.rev();
	}

	// add to transitionMap - check if already present....?
	Transition trans(path.InitialLabels,path.FinalLabels);
	transitionMap[trans] = path.FoundTransitions.front();
	transitionMap[trans.rev()] = path.FoundTransitions.back().rev();

	// job is done at this point....
	return true;
};

bool TammberModel::add_duplicate(NEBPathway &path) {
	LOGGER("TammberModel::add_duplicate")
	if(!path.duplicate) return false;

	Transition trans(path.InitialLabels,path.FinalLabels);

	auto ep = StateEdges.find(CanonTrans(trans));

	bool forwards(CanonTrans(trans).first==ep->first.first);
	if(!forwards) LOGGER("TammberModel::add_duplicate : NEBPathway not parallel to StateEdge!")

	std::pair<SymmLabelPair,PointShiftSymmetry> match;

	bool matched=false;
	for(auto &trp_sy: path.equivalent_transitions) { // ( (LabelPair,LabelPair), PointShiftSymmetry
		SymmLabelPair tl = NonCanonTrans(trp_sy.first);
		if(ep->second.connections.find(tl)!=ep->second.connections.end()) {
			matched=true; // first match
			match = std::make_pair(tl,trp_sy.second);
			// map central transition
			ep->second.self_edge_map[NonCanonTrans(trans)] = match;
		}
	}
	if(matched) {
		for(auto &trp_sy: path.equivalent_transitions) { // ( (LabelPair,LabelPair), PointShiftSymmetry
			SymmLabelPair tl = NonCanonTrans(trp_sy.first);
			if(ep->second.connections.find(tl)!=ep->second.connections.end() and tl!=match.first) {
				// if S1 : T->E => iS1 : E->T. SM : T -> M => iS1.SM : E->M
				PointShiftSymmetry iop = trp_sy.second.inverse();
				ep->second.self_edge_map[tl] =	std::make_pair(match.first,match.second.compound(iop));
				// assimilte and delete
				auto cit = ep->second.connections.find(tl);
				ep->second.connections.find(match.first)->second.assimilate(cit->second);
				ep->second.connections.erase(cit);
			}
		}
		return true;
	} else { // ensure correct placeholder but otherwise try again....?
		path.saddleE = MAX_BARRIER + path.initialE;
		path.Ftol = 1.0;
		path.dX = MSD_THRESH * 2.0;
		path.dXmax = MSD_THRESH * 2.0;
		path.SaddleLabels = path.InitialLabels;
		return false;
	}
	return true;
};

// from NEB
void TammberModel::add_pathway(NEBPathway &path) {
	LOGGER("TammberModel::add_pathway")
	if(!path.valid) {
		std::string err_msg = "";
		if(path.mismatch)
			err_msg = "MISMATCHED / ERRONEOUS LABELS IN PATH! \n"
								"CURRENTLY REMOVING FROM REQUESTED AND PENDING\n"
								"CONSIDER INCREASING MINIMIZATION TOLERANCES\n";
		err_msg += "ADDING TO FAILED NEBS";
		LOGGERA("TammberModel::add_pathway : \nINVALID PATH! "<<path.info_str())

		Transition ftrans(path.InitialLabels,path.FinalLabels); // requested labels
		auto fep = StateEdges.find(CanonTrans(ftrans)); // corresponding edge

		if(fep!=StateEdges.end()) {
			// In any case, remove from requestedNEBS
			if(fep->second.requestedNEBS.find(NonCanonTrans(ftrans))!=fep->second.requestedNEBS.end())
				fep->second.requestedNEBS.erase(fep->second.requestedNEBS.find(NonCanonTrans(ftrans)));

			// In any case, remove from pending NEBS
			if(fep->second.pendingNEBS.find(NonCanonTrans(ftrans))!=fep->second.pendingNEBS.end())
				fep->second.pendingNEBS.erase(fep->second.pendingNEBS.find(NonCanonTrans(ftrans)));

			// Add to failed NEBS TO BE RENAMED!
			fep->second.requestedPMS.insert(NonCanonTrans(ftrans));
		}
		// add these anyway?
		//add_transitionMaps(path);
		return;
	}

	// first off, incorporate the symmetry results
	add_symmetries(path);

	// If NEB found multiple jumps- we request the jumps and add maps
	bool multijump = add_transitionMaps(path);

	Transition trans(path.InitialLabels,path.FinalLabels);
	SymmLabelPair tl = NonCanonTrans(trans);
	add_edge(CanonTrans(trans)); // find / create edge - it should already exist!!
	auto ep = &(StateEdges.find(CanonTrans(trans))->second);

	// should always be forwards, no?
	bool forwards = bool(CanonTrans(trans).first==StateEdges.find(CanonTrans(trans))->first.first);
	if(!forwards) LOGGER("TammberModel::add_pathway : NEBPathway not parallel to StateEdge! Ignoring as likely bug")

	if(path.valid && !multijump) { // we have a normal result
		transitionMap[trans] = trans; // i.e. identity. This is only if multijump==false
		transitionMap[trans.rev()] = trans.rev();

		// If NEB found equivalent transitions- incorporate
		bool duplicate = add_duplicate(path);

		if (!duplicate) {
			// overwrite by default
			if(ep->connections.find(tl)!=ep->connections.end())
				ep->connections.erase(tl);
				Connection conn(path);
				ep->connections.insert(std::make_pair(tl,conn));
				PointShiftSymmetry op; op.valid=true;
				ep->self_edge_map[tl] = std::make_pair(tl,op); // i.e. identity operator
		}
	}

	// either way, all is done
	if(ep->pendingNEBS.find(tl)!=ep->pendingNEBS.end())
		ep->pendingNEBS.erase(ep->pendingNEBS.find(tl));

	if(ep->requestedNEBS.find(tl)!=ep->requestedNEBS.end())
		ep->requestedNEBS.erase(ep->requestedNEBS.find(tl));

	ep->completedNEBS.insert(tl);

	// what is the "canonical" noncanoncial transition label?
	for(auto &pta: pendingTADS) if(pta.second.size()>0) {

		Transition ptrans = transitionMap[pta.first];

		ep = &(StateEdges.find(CanonTrans(ptrans))->second);

		forwards=bool(CanonTrans(ptrans).first==StateEdges.find(CanonTrans(ptrans))->first.first);

		tl = ep->self_edge_map.at(NonCanonTrans(ptrans)).first; // parallel even if ptrans not

		LOGGER("TammberModel::add_pathway ep->connections.find(tl)!=ep->connections.end()")
		if(ep->connections.find(tl)!=ep->connections.end()) {
			LOGGER("TammberModel::add_pathway for(auto &pseg: pta.second)")
			for(auto &pseg: pta.second) { // list of TADSegments
				double fpt = StateVertices.find(ptrans.first.first)->second.state_time(pseg.temperature);
				ep->connections.find(tl)->second.add_jump(pseg.temperature,fpt,forwards);
			}
			LOGGER("TammberModel::add_pathway pta.second.clear()")
			pta.second.clear();
		}
	}
	LOGGER("TammberModel::add_pathway pta=pendingTADS.erase(pta);")
	for(auto pta=pendingTADS.begin();pta!=pendingTADS.end();) {
		if(pta->second.size()==0) pta=pendingTADS.erase(pta);
		else pta++;
	}
};

std::string TammberModel::info_str(bool seg){
	LOGGER("TammberModel::info_str")
	std::string res;
	res = "TammberModel status:\n";
	double st=0.0;
	unsigned blocks=0,overhead=0,nebc=0,jc=0;
	double minE = 99999999.0;
	Label minL;
	for (auto &v: StateVertices) if(minE>v.second.energy) {
		minE = v.second.energy;
		minL = v.first;
	}
	res += "Minimum Energy State : "+std::to_string(minL)+" ("+std::to_string(minE)+")\n\nVerticies:\n";

	for (auto &v: StateVertices) {
		if(!allow_allocation(v.first) and !include_shallow_states) continue;
		UnknownRate ku;
		unknown_rate(v.first,ku);
		res += v.second.info_str(targetT,minE,ku.unknown_rate);
		if (seg) {
			res += "Self Symmetries: "+std::to_string( v.second.self_symmetries.size())+"\n";
			for(auto &ss: v.second.self_symmetries) res+=ss.info_str();
		}
		st += v.second.target_state_time;
		blocks += v.second.duration;
		overhead += v.second.overhead;
	}
	res += "Edges:\n";
	for(auto &e: StateEdges) {
		if(!allow_allocation(e.first.first) and !include_shallow_states) continue;
		if(!allow_allocation(e.first.second) and !include_shallow_states) continue;
		res+=e.second.info_str(seg);
		nebc += e.second.connections.size();
		jc += e.second.self_edge_map.size();
	}
	res += "Total t:"+std::to_string(st)+"ps, blocks:"+std::to_string(blocks)+", overhead:"+std::to_string(overhead);
	res += ", NEBS:"+std::to_string(nebc)+", jumps:"+std::to_string(jc)+"\n";

	if(pendingTADS.size()>0) {
		res+="\n  Pending TADS: \n";
		for(auto &ta: pendingTADS) {
			res += "  "+std::to_string(ta.first.first.first)+","+std::to_string(ta.first.first.second);
			res += " -> "+std::to_string(ta.first.second.first)+","+std::to_string(ta.first.second.second);
			res += " : "+std::to_string(ta.second.size())+"\n";
		}
	} else res+="\n  No Pending TADS\n";
	return res;
};

std::map<Label, std::pair<Label,std::set<Label>>> TammberModel::listStates() {
	LOGGER("TammberModel::listStates")
	std::map<Label,std::pair<Label,std::set<Label>>> res;
	for (auto &v: StateVertices) {
		Label lab = v.second.reference_label.second;
		res.insert({v.second.reference_label.first,{lab,{lab}}});
	}
	for (auto &e: StateEdges) {
		LabelPair clab;
		clab.first = e.second.labels.first;
		clab.second = e.second.labels.second;
		for(auto &sem : e.second.self_edge_map) {
			res[clab.first].second.insert(sem.first.first);
			res[clab.first].second.insert(sem.second.first.first);
			res[clab.second].second.insert(sem.first.second);
			res[clab.second].second.insert(sem.second.first.second);
		}
	}
	return res;
};

void TammberModel::write_model(std::string mmfile) {
	LOGGER("TammberModel::write_model")
	std::ofstream res(mmfile.c_str(),std::ofstream::out);
	std::ofstream state_list("StatesToExtract.list",std::ofstream::out);
	std::ofstream pfn_list("PendingFailedNEBS.list",std::ofstream::out);
	pfn_list<<"# InitialLabels FinalLabels 1=Pending/0=Failed\n";
	double st=0.0;
	unsigned blocks=0,overhead=0,nebc=0,jc=0;
	double minE = 99999999.0;
	std::vector<std::string> p_xml; // for parsing
	for (auto v: StateVertices) minE = std::min(minE,v.second.energy);
	res<<"<MarkovModel>\n";

	if(LatticeConstant.compare("None")!=0)
		res<<"  <Lattice>"<<LatticeConstant<<"</Lattice>\n";

	if(PrimitiveUnitCell.compare("None")!=0)
		res<<"  <PrimitiveUnitCell>"<<PrimitiveUnitCell<<"  </PrimitiveUnitCell>\n";

	if(UnitCell.compare("None")!=0)
		res<<"  <UnitCell>"<<UnitCell<<"  </UnitCell>\n";

	if(SuperCell.compare("None")!=0)
		res<<"  <SuperCell>"<<SuperCell<<"  </SuperCell>\n";
	res<<"<GlobalMinEnergy>"<<minE<<"</GlobalMinEnergy>\n";
	res<<"<ResidenceTimeMean>"<<valid_time<<"</ResidenceTimeMean>\n";
	res<<"<ResidenceTimeStd>"<<sqrt(std::fabs(valid_time_sd))*valid_time<<"</ResidenceTimeStd>\n";

	for (auto &v: StateVertices) {
		if(!allow_allocation(v.first) and !include_shallow_states) continue;
		state_list<<v.second.reference_label.first<<" "<<v.second.reference_label.second<<"\n";
		UnknownRate ku;
		unknown_rate(v.first,ku);
		v.second.target_state_time = v.second.state_time(targetT);
		double emin = (log(v.second.target_state_time)+LOG_NU_MIN)*BOLTZ*targetT;
		double dephase_ratio = (double)(v.second.duration)/std::max(1.0,(double)(v.second.duration+v.second.overhead));
		res<<"  <Vertex>\n";
		res<<"    <CanonLabel>"<<v.second.reference_label.first<<"</CanonLabel>\n";
		res<<"    <ReferenceLabel>"<<v.second.reference_label.second<<"</ReferenceLabel>\n";
		res<<"    <Energy>"<<v.second.energy-minE<<"</Energy>\n";
		res<<"    <UnknownRate>"<<ku.unknown_rate<<"</UnknownRate>\n";
		res<<"    <Time>"<<v.second.target_state_time<<"</Time>\n";
		res<<"    <Temperature>"<<targetT<<"</Temperature>\n";
		res<<"    <TADBarrier>"<<emin<<"</TADBarrier>\n";
		res<<"    <NClusters>"<<v.second.clusters<<"</NClusters>\n";
		res<<"    <DephaseRatio>"<<dephase_ratio<<"</DephaseRatio>\n";
		res<<"    <Allocated>"<<int(allow_allocation(v.first))<<"</Allocated>\n";
		res<<"    <Position>"<<v.second.position[0]<<" "<<v.second.position[1]<<" "<<v.second.position[2]<<"</Position>\n";
		if(v.second.self_symmetries.size()>0){
			res<<"    <SelfSymmetries>\n";
			for(auto &ss: v.second.self_symmetries)
				res<<"      "<<ss.operation<<" "<<ss.shift[0]<<" "<<ss.shift[1]<<" "<<ss.shift[2]<<"\n";
			res<<"    </SelfSymmetries>\n";
		}
		res<<"    <SeenEquivalents>\n";
		res<<"      "<<v.second.reference_label.second<<" 0 0 0 0\n";
		for(auto &ss: v.second.state_isomorphisms) {
			if(ss.first==v.second.reference_label.second) continue;
			res<<"      "<<ss.first<<" "<<ss.second.operation<<" ";
			res<<ss.second.shift[0]<<" "<<ss.second.shift[1]<<" "<<ss.second.shift[2]<<"\n";
		}
		res<<"    </SeenEquivalents>\n";
		res<<"  </Vertex>\n";
	}
	for(auto &e: StateEdges) {
		if(!allow_allocation(e.first.first) and !include_shallow_states) continue;
		if(!allow_allocation(e.first.second) and !include_shallow_states) continue;
		for(auto &tstc : modelEdgeParams(e.first) ) {
			state_list<<e.first.first<<" "<<tstc.first.first<<"\n";
			state_list<<e.first.second<<" "<<tstc.first.second<<"\n";
			res<<"  <Edge>\n";
			res<<"    <CanonLabels>"<<e.first.first<<" "<<e.first.second<<"</CanonLabels>\n";
			res<<"    <ReferenceLabels>"<<tstc.first.first<<" "<<tstc.first.second<<"</ReferenceLabels>\n";
			res<<"    <Barriers>"<<tstc.second.first[0]<<" "<<tstc.second.first[1]<<"</Barriers>\n";
			res<<"    <PreFactors>"<<tstc.second.first[2]<<" "<<tstc.second.first[3]<<"</PreFactors>\n";
			int dynamic_estimate=0;
			if(tstc.second.first[2]<0.0 and tstc.second.first[3]<0.0) dynamic_estimate=1;
			if(!sim_conn) dynamic_estimate=0;
			res<<"    <DynamicEstimate>"<<dynamic_estimate<<"</DynamicEstimate>\n";
			res<<"    <dX>"<<tstc.second.first[4]<<"</dX>\n";
			res<<"    <Ftol>"<<tstc.second.first[5]<<"</Ftol>\n";
			if(tstc.second.second.valid) {
				res<<"    <TransitionSymmetry>"<<tstc.second.second.operation<<" ";
				res<<tstc.second.second.shift[0]<<" "<<tstc.second.second.shift[1]<<" ";
				res<<tstc.second.second.shift[2]<<"</TransitionSymmetry>\n";
			}
			res<<"  </Edge>\n";
		}
		/*
		for(auto &sem: e.second.self_edge_map) {
			if (sem.first==sem.second.first) continue;
			if(sem.second.first==tstc.first) {
				res<<"    <SeenEquivalents>";
				res<<sem.first.first<<" "<<sem.first.second<<" ";
				res<<sem.second.second.operation<<" "<<sem.second.second.shift[0];
				res<<" "<<sem.second.second.shift[1]<<" "<<sem.second.second.shift[2];
				res<<"</SeenEquivalents>\n";
			};
		}
		*/
		for(auto pnl : e.second.pendingNEBS)
			pfn_list<<e.first.first<<" "<<pnl.first<<" "<<e.first.second<<" "<<pnl.second<<" 1\n";
		for(auto pnl : e.second.requestedPMS) // TO BE RENAMED!
			pfn_list<<e.first.first<<" "<<pnl.first<<" "<<e.first.second<<" "<<pnl.second<<" 0\n";
	}
	res<<"</MarkovModel>\n";
	res.close();
	state_list.close();
	pfn_list.close();
};

void TammberModel::generateNEBs(std::list<NEBjob> &jobs, int nMax) {
	LOGGER("TammberModel::generateNEBs")
	jobs.clear();
	std::list<LabelPair> existing;

	for(auto e=StateEdges.begin();e!=StateEdges.end();e++) {
		if (e->second.pendingNEBS.size()>0 && e->second.requestedNEBS.size()==0 && (jobs.size()<nMax)) {
			SymmLabelPair el=e->first; // the canonical transition
			existing.clear();
			for(auto en: e->second.connections)
				existing.push_back(std::make_pair(en.first.first,en.first.second)); // noncanonical labels as inside

			for(auto pn = e->second.pendingNEBS.begin(); pn!=e->second.pendingNEBS.end(); pn++) {
				Transition tt(std::make_pair(el.first,pn->first),std::make_pair(el.second,pn->second));
				NEBjob job(tt,existing);
				if(StateVertices.find(el.first)!=StateVertices.end())
					job.InitialSymmetries = StateVertices.find(el.first)->second.self_symmetries;
				if(StateVertices.find(el.second)!=StateVertices.end())
					job.FinalSymmetries = StateVertices.find(el.second)->second.self_symmetries;
				jobs.push_back(job);
				LOGGER("ADDING NEB "<<el.first<<","<<pn->first<<" -> "<<el.second<<","<<pn->second)
				e->second.requestedNEBS.insert(*pn);
				//pn = e->second.pendingNEBS.erase(pn);
				if(e->second.requestedNEBS.size()>0) break; // only one at a time for the moment....
			}
		}
	}
};

void TammberModel::calculate_rates(SymmLabelPair el, Rate &kf, Rate &kb) {
	LOGGER("TammberModel::calculate_rates")
	double aT,tau_term,pf,fpt,N_term,ttau_term,tk,pndE,maxE;
	kf.reset(tadT.size());
	kb.reset(tadT.size());

	if(!allow_allocation(el.first) || !allow_allocation(el.second))
	 	if(!include_shallow_states) return;

	auto isp = StateVertices.find(el.first);
	auto fsp = StateVertices.find(el.second);
	auto ep = StateEdges.find(el);
	if(ep==StateEdges.end()) return;

	bool forwards = bool(el.first==ep->first.first);
	auto kfp = (forwards ? &kf : &kb);
	auto kbp = (forwards ? &kb : &kf);

	LOGGER("MAKING RATES: "<<ep->second.connections.size()<<" ("<<ep->second.pendingNEBS.size()<<" pending)")

	if(ep->second.connections.size()+ep->second.pendingNEBS.size()==0) return; // need saddle points....

	std::pair<double,double> dE;
	pndE = pNEB_Prior;
	maxE = std::max(isp->second.energy,fsp->second.energy);

	for(auto &conn: ep->second.connections) { // all Connections

		conn.second.nu = conn.second.priornu;

		if(conn.second.dX<MSD_THRESH) {
			LOGGER("CONNECTION V SMALL: dX="<<conn.second.dX<<"A; Setting dE=0.05")
			conn.second.saddleE = maxE+0.05;
		}

		dE.first = std::min(MAX_BARRIER+maxE,conn.second.saddleE)-isp->second.energy;
		dE.second = std::min(MAX_BARRIER+maxE,conn.second.saddleE)-fsp->second.energy;

		pndE = std::min(pndE,std::min(dE.first,dE.second));

		if(conn.second.nu.first<1.0e-9) conn.second.nu.first = PRIOR_NU;
		if(conn.second.nu.second<1.0e-9) conn.second.nu.second = PRIOR_NU;

		if(el.first!=el.second){ // Unsymmetric
			// forwards
			tau_term = 0.0;
			for(auto &et: isp->second.elapsedTime) {
				ttau_term = et.second * exp(-std::min(dE.first,MAX_BARRIER)/BOLTZ/double(et.first));
				if(!std::isnan(ttau_term)) tau_term += ttau_term;
			}

			tau_term = 1.0 - tau_term*conn.second.priornu.first/ prefactorCountThresh;
			N_term = 4.0*conn.second.counts.first/prefactorCountThresh;
			conn.second.nu.first = 0.5 * (tau_term + sqrt(tau_term*tau_term + N_term));
			if(std::isnan(conn.second.nu.first)|| conn.second.nu.first<1.0e-4)
				conn.second.nu.first = PRIOR_NU;

			// backwards
			tau_term = 0.0;
			for(auto &et: fsp->second.elapsedTime){
				ttau_term = et.second * exp(-std::min(dE.second,MAX_BARRIER)/BOLTZ/double(et.first));
				if(!std::isnan(ttau_term)) tau_term += ttau_term;
			}
			tau_term = 1.0 - tau_term * conn.second.priornu.second / prefactorCountThresh;
			N_term = 4.0 * conn.second.counts.second / prefactorCountThresh;
			conn.second.nu.second = 0.5 * (tau_term + sqrt(tau_term*tau_term + N_term));
			if(std::isnan(conn.second.nu.second) || conn.second.nu.second<1.0e-4)
				conn.second.nu.second = PRIOR_NU;

		} else {	// Symmetric
			tau_term = 0.0;
			for(auto &et: isp->second.elapsedTime) {
				ttau_term = et.second * exp(-std::min(dE.first,MAX_BARRIER)/BOLTZ/double(et.first));
				if(!std::isnan(ttau_term)) tau_term += ttau_term;
			}
			tau_term = 1.0 - tau_term*conn.second.priornu.first/ prefactorCountThresh;
			N_term = 4.0*conn.second.counts.first/prefactorCountThresh;
			N_term += 4.0*conn.second.counts.second/prefactorCountThresh;
			conn.second.nu.first = 0.5 * (tau_term + sqrt(tau_term*tau_term + N_term));
			if(std::isnan(conn.second.nu.first)|| conn.second.nu.first<1.0e-9)
				conn.second.nu.first = PRIOR_NU;
			conn.second.nu.second = conn.second.nu.first;
		}

		LOGGER("PREFACTORS: "<<conn.second.nu.first<<" "<<conn.second.nu.second)
		LOGGER("BARRIERS: "<<dE.first<<" "<<dE.second)

		tk = conn.second.nu.first * exp(-std::min(dE.first,MAX_BARRIER)/BOLTZ/targetT);
		if(std::isnan(tk)) {
			LOGGER("Rate calculation returns NaN, most likely due to unconverged NEB. Replacing with 0.5eV")
			tk = 10.0 * exp(-0.5/BOLTZ/targetT);
		}
		kfp->target_k_fp.first += tk;

		fpt = conn.second.fpt.first * exp(std::min(dE.first,MAX_BARRIER)/BOLTZ*(1.0/targetT-1.0/conn.second.fpT.first));
		if(std::isnan(fpt)) {
			LOGGER("FPT returns NaN. Replacing with 1/rate")
			fpt = 1.0 / tk;
		}
		if(kfp->target_k_fp.second<0.0 || kfp->target_k_fp.second>fpt) kfp->target_k_fp.second = fpt;

		tk = conn.second.nu.second * exp(-std::min(dE.second,MAX_BARRIER)/BOLTZ/targetT);
		if(std::isnan(tk)) {
			LOGGER("Rate calculation returns NaN, most likely due to unconverged NEB. Replacing with 0.5eV")
			tk = 10.0 * exp(-0.5/BOLTZ/targetT);
		}

		kbp->target_k_fp.first += tk;
		fpt = conn.second.fpt.second * exp(std::min(dE.second,MAX_BARRIER)/BOLTZ*(1.0/targetT-1.0/conn.second.fpT.second));
		if(std::isnan(fpt)) {
			LOGGER("FPT returns NaN. Replacing with 1/rate")
			fpt = 1.0 / tk;
		}
		if(kbp->target_k_fp.second<0.0 || kbp->target_k_fp.second>fpt) kbp->target_k_fp.second = fpt;

		for(size_t ti=0;ti<tadT.size();ti++) {
			tk = conn.second.nu.first * exp(-std::min(dE.first,MAX_BARRIER)/BOLTZ/tadT[ti]);
			if(std::isnan(tk)) tk = 10.0 * exp(-0.5/BOLTZ/tadT[ti]);
			kfp->tad_k_fp[ti].first += tk;

			fpt = conn.second.fpt.first * exp(std::min(dE.first,MAX_BARRIER)/BOLTZ*(1.0/tadT[ti]-1.0/conn.second.fpT.first));
			if(std::isnan(fpt))	fpt = 1.0 / tk;
			if(kfp->tad_k_fp[ti].second<0.0 || kfp->tad_k_fp[ti].second>fpt) kfp->tad_k_fp[ti].second = fpt;


			tk = conn.second.nu.second * exp(-std::min(dE.second,MAX_BARRIER)/BOLTZ/tadT[ti]);
			if(std::isnan(tk)) tk = 10.0 * exp(-0.5/BOLTZ/tadT[ti]);
			kbp->tad_k_fp[ti].first += tk;
			fpt = conn.second.fpt.second * exp(std::min(dE.second,MAX_BARRIER)/BOLTZ*(1.0/tadT[ti]-1.0/conn.second.fpT.second));
			if(std::isnan(fpt))	fpt = 1.0 / tk;
			if(kbp->tad_k_fp[ti].second<0.0 || kbp->tad_k_fp[ti].second>fpt) kbp->tad_k_fp[ti].second = fpt;
		} // tadT
	} // connections

	// Experimental feature - if a NEB has been requested but not returned,
	// we estimate it's value
	if(sim_conn && ep->second.pendingNEBS.size()>0) { // just one connection
		dE.first = pndE+maxE-isp->second.energy;
		dE.second = pndE+maxE-fsp->second.energy;
		for(size_t ti=0;ti<tadT.size();ti++) {
			kfp->tad_k_fp[ti].first += exp(-dE.first/BOLTZ/tadT[ti]);
			if(kfp->tad_k_fp[ti].second<0.0 || kfp->tad_k_fp[ti].second>exp(dE.first/BOLTZ/tadT[ti]))
				kfp->tad_k_fp[ti].second = exp(dE.first/BOLTZ/tadT[ti]);
			kbp->tad_k_fp[ti].first += exp(-dE.second/BOLTZ/tadT[ti]);
			if(kbp->tad_k_fp[ti].second<0.0 || kbp->tad_k_fp[ti].second>exp(dE.second/BOLTZ/tadT[ti]))
				kbp->tad_k_fp[ti].second = exp(dE.second/BOLTZ/tadT[ti]);
		} // tadT
		kfp->target_k_fp = kfp->tad_k_fp[0];
		kbp->target_k_fp = kbp->tad_k_fp[0];
	} // pendingNEBS

};// canonical labels

// this should be abstracted away as large portion used twice...

// modelEdgeParams- list(SymmLabelPair,([dE_f,dE_b,nu_f,nu_b], PointShiftSymmetry))
std::list<std::pair<SymmLabelPair,std::pair< std::array<double,6>,PointShiftSymmetry >>>
	TammberModel::modelEdgeParams(SymmLabelPair el) {

	LOGGER("TammberModel::modelEdgeParams")

	std::list<std::pair<SymmLabelPair,std::pair< std::array<double,6>,PointShiftSymmetry >>> res;

	auto isp = StateVertices.find(el.first);
	auto fsp = StateVertices.find(el.second);
	auto ep = StateEdges.find(el);

	if(isp==StateVertices.end() || fsp==StateVertices.end() || ep==StateEdges.end()) return res;
	bool forwards = bool(el.first==ep->first.first);
	bool return_pending = sim_conn;
	return_pending *=  (ep->second.pendingNEBS.size()>0);
	return_pending *= (ep->second.connections.size()==0);
	if(ep->second.connections.size()==0 and !sim_conn) return res;

	std::array<double,6> sres;
	LabelPair slab;
	PointShiftSymmetry op;
	std::pair<double,double> dE;
	double maxE = std::max(isp->second.energy,fsp->second.energy) + 0.05;
	for(auto &conn: ep->second.connections) { // all Connections
		double sE = std::min(maxE+MAX_BARRIER,conn.second.saddleE);
		for(auto EE : conn.second.energies) sE = std::max(sE,EE);
		dE.first = sE -  isp->second.energy;
		dE.second = sE - fsp->second.energy;
		sres[int(!forwards)]=dE.first;
		sres[int(forwards)] =dE.second;
		sres[2+int(!forwards)] = conn.second.nu.first < 1.0e-4 ? PRIOR_NU : conn.second.nu.first;
		sres[2+int(forwards)] = conn.second.nu.second < 1.0e-4 ? PRIOR_NU : conn.second.nu.second;
		sres[4] = conn.second.dX;
		sres[5] = conn.second.Ftol;
		slab.first = (forwards ? conn.first.first : conn.first.second);
		slab.second = (forwards ? conn.first.second : conn.first.first);
		if(el.first==el.second && conn.second.self_transitions.size()>0) {
			op.operation = conn.second.self_transitions.begin()->second.operation;
			op.shift = conn.second.self_transitions.begin()->second.shift;
			op.valid=true;
		} else {
			op.valid=false;
			op.operation = 0;
			op.shift = {0.,0.,0.};
		}
		res.push_back(std::make_pair(SymmLabelPair(slab),std::make_pair(sres,op)));
	}

	// return incomplete transitions if no completed NEBS for this edge
	if(return_pending) {
		double maxT = tadT[tadT.size()-1];
		auto pslab = *(ep->second.pendingNEBS.begin());
		slab.first = (forwards ? pslab.first : pslab.second);
		slab.second = (forwards ? pslab.second : pslab.first);
		dE.first = std::max(0.25,log(10.0*isp->second.state_time(maxT))*BOLTZ*maxT) + isp->second.energy;
		dE.second = std::max(0.25,log(10.0*fsp->second.state_time(maxT))*BOLTZ*maxT) + fsp->second.energy;
		sres[int(!forwards)]=std::max(dE.first,dE.second)-isp->second.energy;
		sres[int(forwards)] =std::max(dE.first,dE.second)-fsp->second.energy;
		sres[2+int(!forwards)]=-1.0;
		sres[2+int(forwards)] =-1.0;
		sres[4] = 10.0;
		sres[5] = 100.0;
		op.valid=false;
		op.operation = 0;
		op.shift = {0.,0.,0.};
		res.push_back(std::make_pair(SymmLabelPair(slab),std::make_pair(sres,op)));
	}
	return res;
};

bool TammberModel::allow_allocation(Label lab) {

	if(StateVertices.find(lab)==StateVertices.end()) {
		LOGGER("SUPPRESSING ALLOCATION TO "<<lab<<" : NOT FOUND")
		return false;
	}

	auto v = &(StateVertices.find(lab)->second);

	// false if ClusterThresh==0
	if(int(ClusterThresh)>0 and v->clusters>ClusterThresh) {
		LOGGER("SUPPRESSING ALLOCATION TO "<<lab<<" : TOO MANY CLUSTERS")
		return false;
	}

	/*
		if duration / (duration+overhead) < DephaseThresh : cancel
		rearrange: duration < alpha * overhead
		where alpha = DephaseThresh/(1-DephaseThresh)
	*/
	if (DephaseThresh>0.0) {
		double alpha = DephaseThresh/std::max(TINY,1.0-DephaseThresh);
		if( bool(v->overhead>=5) && bool(v->duration < alpha*v->overhead ) ) {
			LOGGER("SUPPRESSING ALLOCATION TO "<<lab<<" : FAILED DEPHASING THRESHOLD")
			return false;
		}
	}

	return true;
};

void TammberModel::unknown_rate(Label lab, UnknownRate &ku) {
	LOGGER("TammberModel::unknown_rate")

	double emin, TAD_time, TAD_min_rate, TAD_kuvar, TAD_gradient;
	double var_correction,opt_rate, tar_rate, max_TAD_gradient=-10.,Tr;
	std::vector< std::pair<double,double> > TAD_ku_kt, k_fp;
	std::vector<double> TAD_gradients;

	// default values- assume small residence time with high certainty
	ku.observed_rate = 10.0;
	ku.unknown_rate = TINY;
	ku.unknown_variance = 2.0*TINY;
	ku.optimal_temperature = tadT[0];
	ku.optimal_temperature_index = 0;
	ku.optimal_rate = TINY;
	ku.min_rate = TINY;

	if(StateVertices.find(lab)==StateVertices.end()) {
		LOGGER("TammberModel::unknown_rate "<<lab<<" NOT CALCULATING"
		" AS STATE DOESNT HAVE VERTEX! ERROR!")
		return;
	}
	// we know this exists
	auto v = &(StateVertices.find(lab)->second);

	if(!allow_allocation(lab)) {
		LOGGERA("TammberModel::unknown_rate "<<lab<<" NOT CALCULATING"
		" DUE TO DEPHASE OR CLUSTER LIMIT")
		return;
	}
	// low temperature time
	v->target_state_time = v->state_time(targetT);

	// TAD barrier
	emin = log(std::max(v->target_state_time,1.0))+LOG_NU_MIN;

	auto exit_label_rate_map = Rates.find(lab);
	auto self_label_rate_map = SelfRates.find(lab);
	bool have_exit_rates = bool(exit_label_rate_map!=Rates.end());
	bool have_self_rates = bool(self_label_rate_map!=SelfRates.end());

	size_t total_rate_count=0;
	if(have_exit_rates) total_rate_count += exit_label_rate_map->second.size();
	if(have_self_rates) total_rate_count += self_label_rate_map->second.size();


	// If we have no rates at all use analytic solution
	if (total_rate_count == 0) {

		tar_rate = std::min(1.0/std::max(1.0,v->target_state_time),max_ku);
		ku.observed_rate = 0.0;
		ku.unknown_rate = tar_rate;
		ku.unknown_variance = tar_rate * tar_rate;
		ku.optimal_temperature = tadT[tadT.size()-1];
		ku.optimal_temperature_index = tadT.size()-1;
		Tr = targetT/tadT[tadT.size()-1];
		TAD_time = v->target_state_time * exp( emin * (Tr-1.0) );
		opt_rate = std::min(1.0/std::max(1.0,TAD_time),max_ku);
		ku.optimal_rate = opt_rate;
		ku.min_rate = tar_rate;

	} else {
		// target T from just the exit rates
		k_fp.clear();
		if(have_exit_rates) for(auto &label_rate : exit_label_rate_map->second)
			k_fp.push_back(label_rate.second.target_k_fp);
		bayes_ku_kuvar(k_fp, ku.unknown_rate, ku.unknown_variance, ku.min_rate, ku.observed_rate, v->target_state_time);

		// TAD optimization using exit and self rates
		TAD_ku_kt.clear();
		for(int ii=0;ii<tadT.size();ii++) {
			k_fp.clear();
			if(have_exit_rates) for(auto &label_rate : exit_label_rate_map->second)
				k_fp.push_back(label_rate.second.tad_k_fp[ii]);
			if(have_self_rates) for(auto &label_rate : self_label_rate_map->second)
				k_fp.push_back(label_rate.second.tad_k_fp[ii]);

			TAD_ku_kt.push_back(std::make_pair(0.0,0.0));

			TAD_time = v->target_state_time * exp( emin * (targetT/tadT[ii]-1.0) );

			TAD_min_rate=0.0; TAD_kuvar=0.0; // not used...

			bayes_ku_kuvar(k_fp, TAD_ku_kt[ii].first, TAD_kuvar, TAD_min_rate, TAD_ku_kt[ii].second, TAD_time);
		}

		// now find best temperature
		TAD_gradients.clear();
		for(int ii=0;ii<tadT.size();ii++) {

			// target quantity : P(dephase in 1 ps) * dku/dt / (dc/dt)

			// dku/dt = min_k_lt * ku_ht + (tau_lt/tau_ht * T_h/T_l - ku_ht/ku_lt)*ku_var_lt
			TAD_gradient = ku.min_rate * TAD_ku_kt[ii].first;
			var_correction = (tadT[ii]/targetT) * exp( emin * (1.0-targetT/tadT[ii]) ) * ku.unknown_variance;
			var_correction += TAD_ku_kt[ii].first / ku.unknown_rate * ku.unknown_variance;
			if(var_correction>0.0 && !safe_opt) TAD_gradient += var_correction;

			// /= (dc/dt)
			TAD_gradient /= 1. + HashCost * TAD_ku_kt[ii].second + (HashCost + NEBCost) * TAD_ku_kt[ii].first;

			// *= P(dephase in 1 ps)
			TAD_gradient *= exp(- (TAD_ku_kt[ii].first + TAD_ku_kt[ii].second) );

			// error catching (typically float overflow)
			if(std::isnan(TAD_gradient)) {
				TAD_gradient = 1.0/v->target_state_time/v->target_state_time;
				LOGGER("Cost return is NaN! state,time,T = "<<lab<<", "<<v->target_state_time<<", "<<tadT[ii]<<" new kc:"<<TAD_gradient)
			}

			TAD_gradients.push_back(TAD_gradient);

			if(TAD_gradient > max_TAD_gradient) max_TAD_gradient = TAD_gradient;
		}

		// only run at higher temperature if *some* sampling has been accumulated
		int optimal_temperature_index=0;
		if((double)(v->duration+v->overhead) >= 5.0) {
			for(int ii=1;ii<tadT.size();ii++) {
				optimal_temperature_index = ii;
				if(TAD_gradients[ii] >= 0.95*max_TAD_gradient) break; // avoid plateaus
			}
		}
		// use unknown rate estimate accounting for self rates
		ku.optimal_temperature = tadT[optimal_temperature_index];
		ku.optimal_temperature_index = optimal_temperature_index;
		ku.optimal_rate = TAD_ku_kt[optimal_temperature_index].first;
		ku.optimal_gradient = TAD_gradients[optimal_temperature_index];

		// final overide check...
		if(ku.optimal_gradient<0.0) { // default to low temperature
			LOGGER("TammberModel::unknown_rate : -VE Max TAD_gradients! "<<lab)
			ku.optimal_temperature = tadT[0];
			ku.optimal_temperature_index = 0;
			ku.optimal_rate = 1.0/v->target_state_time;
			ku.optimal_gradient = 1.0/v->target_state_time;
		}
	}
	if(ku.unknown_rate<0.0) ku.unknown_rate = TINY;
};

void TammberModel::bayes_ku_kuvar(std::vector<std::pair<double,double>> &k_fp,double &ku, double &kuvar, double &min_k, double &tot_k, double _time){
	LOGGER("TammberModel::bayes_ku_kuvar")
	// all ensemble_average for the moment
	std::vector<double> sum_kobs;
	ku=0.0;
	kuvar=0.0;
	min_k=max_ku;
	tot_k = 0.0;
	double t_ko;
	// push back all rates and fptimes, sum total rate and find minimum observed rate
	for(auto &k: k_fp) {
		min_k = std::min(min_k,k.first);
		tot_k += k.first;
	}
	// sort in ascending fptime
	auto fpindex = pairvecsort(k_fp);

	// create terms for analytic expansion
	sum_kobs.clear();
	for(auto fpi: fpindex) {
		t_ko = 0.0; for(auto &k: k_fp) t_ko += (k.first)*(k.first) / (k.first + k_fp[fpi].first);
		sum_kobs.push_back(tot_k-t_ko);
	}
	int nr = sum_kobs.size();

	if (nr > 0) {
		// analytic
		ku = 0.0;
		kuvar = 0.0;

		std::vector<double> coeffs((unsigned)(nr+1),0.);
		// series of N terms a_i, i E [0,N-1]
		// prod_i(x+a_i) = sum_n(c_n x^n), n E [0,N], where c_N == 1
		// x + a_0 => c_0 = a_0, c_1 = 1
		coeffs[0] = sum_kobs[0];
		coeffs[1] = 1.0;
		// for each new term a_N
		// ( x + a_N ) * prod_i(x+a_i) = sum_n( [a_N*c_n]x^n) + sum_n(c_n x^(n+1) )
		for(unsigned i=1;i<nr;i++) {
			coeffs[i+1] = 1.0; // c_(N+1) = 1
			for(unsigned j=i; j>0;j--) coeffs[j] = coeffs[j-1] + sum_kobs[i] * coeffs[j]; // c_n ->  a_N*c_n + c_(n-1)  (n>0)
			coeffs[0] = coeffs[0] * sum_kobs[i]; // c_0 -> a_N * c_0
		}

		std::vector<double> inc_fac((unsigned)(nr+3),0.);
		inc_fac[0] = 1.;
		for (int i=1;i<nr+3;i++) inc_fac[i] = inc_fac[i-1] * (double)(i); // i!

		double tpow=1., norm = 0.;
		for (int i=0;i<nr+1;i++) {
			norm += inc_fac[i] * coeffs[i] / tpow ; // i! / t^{i+1} * C[i] * t
			ku += inc_fac[i+1] * coeffs[i] / tpow; // (i+1)! / t^{i+2} * C[i] * t^2
			kuvar += inc_fac[i+2] * coeffs[i] / tpow; // (i+2)! / t^{i+3} * C[i] * t^3
			tpow *= _time; // t^(i)
		}

		// <ku> = res[0] / t^2 / (norm/t) = res[0]/norm/t
		ku /= norm * _time;
		kuvar /= norm * _time * _time;
		kuvar -= ku*ku;
	} else {
		ku = 1.0/_time;
		kuvar = 1.0/_time/_time;
	}

	if(std::isnan(ku) || std::isnan(kuvar)) {
		ku = 1.0/_time;
		kuvar = 1.0/_time/_time;
	}
	min_k = std::min(ku,min_k);
};

void TammberModel::generateTADs(std::list<TADjob> &jobs, int nMax, bool screen) {
	LOGGER("TammberModel::generateTADs")
	std::map<Label,std::pair<double,double>> weights; // Label : (weight, temperature)
	predict(weights,screen);
	jobs.clear();
	int tot_count=0;
	double sub_weight=0.0;

	//double ww = 1.0 / (double)(weights.size());
	//for(auto wit=weights.begin();wit!=weights.end();wit++) wit->second.first = ww;

	auto keysort = mapkeysort(weights);

	for(auto &k: keysort) {
		sub_weight += weights[k].first;
		if(tot_count++ >= PredictionSize) break;
	}

	if(screen) LOGGERA("=============PREDICTION===============")


	tot_count=0;
	for(auto &k: keysort) { // sorted in descending counts
		auto v = &(StateVertices.find(k)->second);
		auto labs = v->reference_label;
		double temperature = weights.at(k).second;
		auto bl= v->BasinLabels;
		int count = std::max(1,int(weights.at(k).first*nMax / sub_weight));
		std::string line = std::to_string(labs.first)+" "+std::to_string(labs.second);
		line += " "+std::to_string(weights[k].first/sub_weight)+" "+std::to_string(count);
		line += " "+std::to_string(weights[k].second)+"K "+std::to_string(tot_count);

		if(count>0 && tot_count < nMax) {
			//TADjob job(labs,temperature,count,bl); // not adding this in yet
			if(tot_count+count>nMax) count = nMax-tot_count;
			TADjob job(labs,temperature,count);
			jobs.push_back(job);
			line += " ALLOC.";
		}
		if(screen) LOGGERA(line)
		tot_count += count;
	}
	if(screen) LOGGERA("======================================")
};

void TammberModel::predict(std::map<Label,std::pair<double,double>> &weights,bool screen) {

	// (re)calculate all rates TODO make this more efficient....
	Rates.clear();
	SelfRates.clear();

	unsigned bar_count=0;

	for(auto &el: StateEdges) {
		Rate kf,kb;

		calculate_rates(el.first,kf,kb);

		LOGGER(el.first.first<<" "<<el.first.second<<" : "<<kf.target_k_fp.first<<" "<<kb.target_k_fp.first)

		if(kf.target_k_fp.first<1.0e-30 || kb.target_k_fp.first<1.0e-30) continue;

		if(el.first.first==el.first.second) { // SelfRates

			if(SelfRates.find(el.first.first)==SelfRates.end())
				SelfRates.insert( std::make_pair(el.first.first, *(new std::map<Label,Rate>) ) );

			SelfRates.at(el.first.first).insert(std::make_pair(el.first.second,kf));

		} else {

			if(Rates.find(el.first.first)==Rates.end())
				Rates.insert( std::make_pair(el.first.first, *(new std::map<Label,Rate>) ) );

			Rates.at(el.first.first).insert(std::make_pair(el.first.second,kf));

			if(Rates.find(el.first.second)==Rates.end())
				Rates.insert( std::make_pair(el.first.second, *(new std::map<Label,Rate>) ) );

			Rates.at(el.first.second).insert(std::make_pair(el.first.first,kb));

			bar_count++;
		}
	}

	std::vector<Label> IndexLabel;
	std::map<Label,int> LabelIndex;

	unsigned ms = 0;
	UnknownRates.clear();

	for(auto &v: StateVertices) { // make all rates
		UnknownRate ku;
		unknown_rate(v.first,ku);
		UnknownRates.insert(std::make_pair(v.first,ku));
		// being very careful- only when unknown_rate exists do we do anything...
		if(ku.unknown_rate>0.0 && allow_allocation(v.first)) {
			IndexLabel.push_back(v.first);
			LabelIndex.insert(std::make_pair(v.first,IndexLabel.size()-1));
			ms++;
		} else if(screen) LOGGERA("STATE "<<v.first<<" NOT ALLOCATED");
	}

	bool solved_one = false;
	bool solved_two = false;

	std::vector<double> allocation(ms,0.0), pabs(ms,0.0), ku(ms,0.0), kuvar(ms,0.0), times(ms,0.0);



	if(screen) LOGGERA("bar_count , state_count: "<<bar_count<<" , "<<ms)

	if(bar_count>0 && ms>0) {
		// make matricies
		LowTR.setZero(); LowTR.resize(ms,ms);
		LowTRT.setZero(); LowTRT.resize(ms,ms);
		Eigen::SparseLU<SpMatCM> solver;

		EigVec rhoinit(ms),PiQ(ms),iQI(ms),ones(ms),flat(ms),selfk(ms),rhoboltz(ms);
		for(int si=0; si<ms; si++) {
			rhoinit[si] = 0.0;
			PiQ[si] = 0.0;
			iQI[si] = 0.0;
			ones[si] = 0.0;
			flat[si] = 0.0;
			selfk[si] = 0.0;
		}

		// Boltzmann
		int min_clust=100000;
		double rho_minE = 10.,rho_norm=0.0;
		for(int si=0; si<ms; si++) {
			if(StateVertices.at(IndexLabel[si]).energy<rho_minE) rho_minE = StateVertices.at(IndexLabel[si]).energy;
			if(StateVertices.at(IndexLabel[si]).clusters<min_clust) min_clust = StateVertices.at(IndexLabel[si]).clusters;
		}
		if(screen) LOGGERA("Minimum Energy: "<<rho_minE<<"eV")
		if(screen) LOGGERA("Minimum Cluster Count: "<<min_clust)

		for(int si=0; si<ms; si++) {
			rhoboltz[si] = exp(-targetB*(StateVertices.at(IndexLabel[si]).energy - rho_minE));
			rho_norm += rhoboltz[si];
		}
		for(int si=0; si<ms; si++) rhoboltz[si] /= rho_norm;
		rho_norm=0.0;
		for(int si=0; si<ms; si++) rho_norm += rhoboltz[si];
		if(screen) LOGGERA("Boltzmann Norm: "<<rho_norm)
		int i,j;
		double k;
		std::vector<Eigen::Triplet<double>> LowTR_trip, LowTRT_trip;
		LowTR_trip.reserve(bar_count);
		LowTRT_trip.reserve(bar_count);

		// go through all states
		for(auto &l_i: LabelIndex) {
			i = l_i.second;

			if(SelfRates.find(l_i.first)!=SelfRates.end())
				for(auto &jr: SelfRates.find(l_i.first)->second)
					selfk[i] += jr.second.target_k_fp.first;

			auto ku = UnknownRates.find(l_i.first)->second;

			ku.observed_rate = 0.0;

			for(auto &jr: Rates.find(l_i.first)->second) { // std::map<Label,Rate>
				k = jr.second.target_k_fp.first;
				if(LabelIndex.find(jr.first)==LabelIndex.end()) {
					ku.observed_rate += k;
					continue;
				}
				j = LabelIndex[jr.first];
				if(i!=j) {
					LowTR_trip.push_back(Eigen::Triplet<double>(i,j,k));
					LowTRT_trip.push_back(Eigen::Triplet<double>(j,i,k));
					ku.observed_rate += k;
				} else {
					LOGGER("i<->i in Rates()! : "<<l_i.first<<" "<<jr.first<<" "<<i)
					selfk[i] += k;
				}
			}
			k = ku.unknown_rate+ku.observed_rate;
			LowTR_trip.push_back(Eigen::Triplet<double>(i,i,-k));
			LowTRT_trip.push_back(Eigen::Triplet<double>(i,i,-k));
		}

		LowTR.setFromTriplets(LowTR_trip.begin(), LowTR_trip.end());
		LowTRT.setFromTriplets(LowTRT_trip.begin(), LowTRT_trip.end());


		// populate initial distribution;
		double st_norm, max_allo, min_allo;
		bool found_init=false,solved_one,solved_two;
		for(int si=0; si<ms; si++) {
			ones[si] = 1.0;
			flat[si] = 1.0/ms;
			if(IndexLabel[si]==InitialLabels.first){
				rhoinit[si] = 1.0;
				found_init = true;
			}
			times[si] = StateVertices.at(IndexLabel[si]).target_state_time;
			ku[si] = UnknownRates.at(IndexLabel[si]).unknown_rate;
			kuvar[si] = UnknownRates.at(IndexLabel[si]).unknown_variance;
		}
		// overwrite rho_init if Boltzmann chosen or initialstate not available
		if (!found_init) {
			for(int si=0; si<ms; si++) rhoinit[si] = rhoboltz[si];
			if(screen) LOGGERA("INITIAL STATE NOT FOUND; USING BOLTZ INITIAL DIST")
		}

		//LOGGER("Q:\n"<<LowTR<<"\nQ.T:\n"<<LowTRT<<"\nSelfQ.diagonal:\n"<<selfk)

		// Prepare Matricies
		LowTR.makeCompressed();
		LowTRT.makeCompressed();



		// Solve for PiQ
		solved_one = true;
		if(screen) LOGGERA("STARTING SPARSE LINEAR SOLVE FOR P.iQ = iQT.P")
		solver.compute(LowTRT);
		if(solver.info()==Eigen::Success) {
			if(screen) LOGGERA("FACTORIZATION DONE")
			if(RhoInitFlavor==1) PiQ = solver.solve(rhoboltz);
			else if (RhoInitFlavor==2) PiQ = solver.solve(flat);
			else PiQ = solver.solve(rhoinit);
			if(solver.info()!=Eigen::Success) solved_one=false;
			else if(screen) LOGGERA("SOLVED FOR P.iQ")
		} else solved_one=false;

		// Solve for iQ.1
		solved_two = true;
		if(screen) LOGGERA("STARTING SPARSE LINEAR SOLVE FOR iQ.1")
		solver.compute(LowTR);
		if(solver.info()==Eigen::Success) {
			if(screen) LOGGERA("FACTORIZATION DONE")
			iQI = solver.solve(ones);
			if(solver.info()!=Eigen::Success) solved_two=false;
			else if(screen) LOGGERA("SOLVED FOR iQl")
		} else solved_two=false;



		if (solved_one) {
			st_norm = 0.0;
			for(int si=0; si<ms; si++) {
				allocation[si] = std::max(0.0,iQI[si] * PiQ[si] * UnknownRates.at(IndexLabel[si]).optimal_gradient);
				st_norm += allocation[si];
			}
			max_allo=0.;
			min_allo=0.;
			for(int si=0; si<ms; si++) {
				allocation[si] /= st_norm;
				if(max_allo<allocation[si]) max_allo = allocation[si];
				if(min_allo>allocation[si]) min_allo = allocation[si];
			}
			if(screen) LOGGERA("SOLVED; MAX/MIN ALLOCATION WEIGHT: "<<max_allo<<"/"<<min_allo)

			if(min_allo < 0. || max_allo <= 0.) {
				if(screen) LOGGERA("NEGATIVE/NONPOSITIVE ALLOCATION WEIGHTS!")
				solved_one = false;
				solved_two = false;
			}
		}

		if(solved_two) {
			valid_time = 0.0;
			st_norm = 0.0;
			for(int si=0; si<ms; si++) {
				valid_time += -PiQ[si];
				pabs[si] = std::max(-PiQ[si] * ku[si],0.0);
				st_norm += pabs[si];
			}
			for(int si=0; si<ms; si++) pabs[si] /= st_norm;

			if(RhoInitFlavor==1) {
				valid_time = -rhoboltz.dot(iQI);
			} else if(RhoInitFlavor==2) {
				valid_time = -flat.dot(iQI);
			}	else {
				valid_time = -rhoinit.dot(iQI);
			}
			if(valid_time>0.0) valid_time_sd = PiQ.dot(iQI)*2./valid_time/valid_time-1.;
			else valid_time_sd=0.0;
			if(screen) {
			  LOGGERA("SOLVED; VALIDITY TIME MEAN (RIF="<<RhoInitFlavor<<"): "<<valid_time<<"ps, "<<"VARIANCE/MEAN^2: "<<valid_time_sd)
			  LOGGERA("VALIDITY TIME (BOLTZ): "<<-rhoboltz.dot(iQI)<<"ps")
			  LOGGERA("VALIDITY TIME (DELTA) : "<<-rhoinit.dot(iQI)<<"ps")
			  LOGGERA("VALIDITY TIME (ONES) : "<<-flat.dot(iQI)<<"ps")
			}
		} else {
			for(int si=0; si<ms; si++) valid_time += 1.0 / ku[si];
			if(screen) LOGGERA("COULD NOT SOLVE; SUM OF INVERSE UNKNOWN RATES: "<<valid_time<<"ps")
		}

		LOGGER("ku, ku*t, kuvar*t, pabs, alloc.:")
		for(int si=0; si<ms; si++)
			LOGGER(ku[si]<<", "<<ku[si]*times[si]<<", "<<kuvar[si]*times[si]<<", "<<pabs[si]<<", "<<allocation[si])
		LOGGER("rhoboltz, rhoinit:")
		for(int si=0; si<ms; si++) LOGGER(rhoboltz[si]<<", "<<rhoinit[si])



		if(solved_one && solved_two) {
			if (AllocScheme==1) {
				if(screen) LOGGERA("USING Pabs ALLOCATION")
				for(int si=0;si<ms;si++) {
					weights.insert(std::make_pair(IndexLabel[si],std::make_pair(pabs[si],UnknownRates.at(IndexLabel[si]).optimal_temperature)));
				}
			} else { // default to 0
				if(screen) LOGGERA("USING GRADIENT ALLOCATION")
				for(int si=0;si<ms;si++) {
					weights.insert(std::make_pair(IndexLabel[si],std::make_pair(allocation[si],UnknownRates.at(IndexLabel[si]).optimal_temperature)));
				}
			}

			return;
		}

		if(solved_two) {
			if(screen) LOGGERA("COULD NOT SOLVE GRADIENT: USING Pabs ALLOCATION")
			for(int si=0;si<ms;si++)
				weights.insert(std::make_pair(IndexLabel[si],std::make_pair(pabs[si],tadT[0])));
			return;
		}
	}
	if(screen) LOGGERA("COULD NOT SOLVE OR NETWORK TOO SMALL; USING 1/time ALLOCATION WEIGHTS")
	double vt=0.0,Zsum=0.0,minE=10.0,rho_norm=0.0;
	for(auto &v: StateVertices) minE = std::min(minE,v.second.energy);
	for(auto &v: StateVertices) rho_norm += exp(-targetB*(v.second.energy-minE));
	for(auto &v: StateVertices)
		vt += exp(-targetB*(v.second.energy-minE))/rho_norm/std::max(1.0,v.second.target_state_time);
	valid_time = 1.0/vt;
	valid_time_sd = 1.0;
	if(screen) LOGGERA("APPROXIMATE VALIDITY TIME: "<<1.0/vt<<"ps")

	Label lab;
	double tw=0.0,mt=0.0,t,temp,emin,sw=0.0;

	for(auto &v: StateVertices) {
		t = std::max(1.0,v.second.target_state_time);
		if(allow_allocation(v.first)) {
			tw += 1.0/t;
			mt = std::max(mt,t);
		}
	}

	for(auto &v: StateVertices) {
		temp = log(PRIOR_NU*std::max(1.0,v.second.target_state_time))*targetT;
		temp = std::max(tadT[0],temp/log(2.0*PRIOR_NU));
		for(auto tt:tadT) if(temp<tt) {
			temp=tt;
			break;
		}
		temp = std::min(*std::prev(tadT.end()),temp);
		lab = v.first;
		t = std::max(1.0,v.second.target_state_time);
		sw = allow_allocation(lab) ? 1.0/tw/t : 0.0;
		weights.insert(std::make_pair(lab,std::make_pair(sw,temp)));
	}

	return;
};



void TammberModel::deleteVertex(Label deleteVertex) {
	for(auto iv=StateVertices.begin();iv!=StateVertices.end();) {
		if(iv->first==deleteVertex) iv = StateVertices.erase(iv);
		else iv++;
	}
	for(auto ie=StateEdges.begin();ie!=StateEdges.end();) {
		if (ie->first.first==deleteVertex || ie->first.second==deleteVertex) {
			ie = StateEdges.erase(ie);
		} else ie++;
	}
	for(auto pt=pendingTADS.begin();pt!=pendingTADS.end();) {
		if(pt->first.first.first==deleteVertex || pt->first.second.first==deleteVertex) pt = pendingTADS.erase(pt);
		else pt++;
	}
};

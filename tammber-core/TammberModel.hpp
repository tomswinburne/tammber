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



#ifndef __Tammber__MM__
#define __Tammber__MM__

#include "Types.hpp"
#include "TammberTypes.hpp"
#include "CustomTypedefs.hpp" // for task Mapper
#include "Log.hpp"

#include <set>
#include <limits>
#include <algorithm>
#include <numeric>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat;
typedef Eigen::SparseMatrix<double> SpMatCM;
typedef Eigen::VectorXd EigVec;
typedef Eigen::Triplet<double> EigTriplet; // col,row,double



// symmetric label pair for MarkovModel

// edge and vertex structs for MarkovModel


/*
// NOT YET IMPLEMENTED
struct BasinVertex {
  Label BasinLabel; // canonical. set from start
  // TODO: BasinLabels is for basin building. But if a known_labels state has lower energy, it should become the new reference.
  std::unordered_map<LabelPair,StateVertex> BasinLabels;
};
struct BasinEdge {
  SymmLabelPair BasinLabels;
};
*/

// TODO template template this for any container type
// argsort implementation for weight map. returns keys
// typically V=std::pair , where operator< ==  lhs.first<rhs.first || (!(rhs.first<lhs.first) && lhs.second<rhs.second)
template <typename L, typename V>
std::vector<L> mapkeysort(const std::map<L,V> &m) {
	std::vector<L> lv;
	for(auto &mm: m) lv.push_back(mm.first);
	// sort indexes based on comparing values in v
	sort(lv.begin(), lv.end(), [&m](L i1, L i2) {return m.at(i1) > m.at(i2);});
	return lv;
};

// sort ascending according to second index of pair by default
template <typename F, typename S>
std::vector<size_t> pairvecsort(const std::vector<std::pair<F,S>> &m,bool ssi=true) {
	std::vector<size_t> idx(m.size());
	iota(idx.begin(), idx.end(), 0);
	// sort indexes based on comparing values in v
	if(ssi)	sort(idx.begin(), idx.end(), [&m](size_t i1, size_t i2) {return m[i1].second < m[i2].second;});
	else sort(idx.begin(), idx.end(), [&m](size_t i1, size_t i2) {return m[i1].first < m[i2].first;});
	return idx;
};


template<typename F, typename S>
std::pair<S,F> revp(std::pair<F,S> other) {
	return std::make_pair(other.second,other.first);
};

struct SymmLabelPair {
	Label first, second;
	SymmLabelPair();
	SymmLabelPair(LabelPair other);
	// comparitors for std::set
	bool operator==( const SymmLabelPair& other) const {
		bool res= bool(std::max(first,second)==std::max(other.first,other.second));
		res *= bool(std::min(first,second)==std::min(other.first,other.second));
		return res;
	};
	bool operator!=( const SymmLabelPair& other) const {
		bool res= bool(std::max(first,second)!=std::max(other.first,other.second));
		res += bool(std::min(first,second)!=std::min(other.first,other.second));
		return res;
	};
	bool operator<( const SymmLabelPair& other) const {
		if(std::max(first,second)<std::max(other.first,other.second)) return true;
		if(std::max(first,second)==std::max(other.first,other.second))
			if(std::min(first,second)<std::min(other.first,other.second)) return true;
		return false;
	};
	SymmLabelPair rev();
	LabelPair ns();
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & first;
		ar & second;
	};
};


SymmLabelPair CanonTrans(Transition trlab);
SymmLabelPair NonCanonTrans(Transition trlab);

// for submission
struct NEBjob {
	Transition TargetTransition;
	bool pairmap;
	std::set<PointShiftSymmetry> InitialSymmetries, FinalSymmetries;
	Label Saddle; // optional atm
	std::list<LabelPair> ExistingPairs;
	NEBjob();
	NEBjob(Transition tt);
	NEBjob(Transition tt,std::list<LabelPair> epl);
};

struct TADjob {
	LabelPair InitialLabels;
	double temperature;
	int nInstances;
	std::map<LabelPair,double> BasinLabels; // currently only recorded, not implemented
	TADjob();
	TADjob(LabelPair il,double t,int n=1);
	TADjob(LabelPair il,double t,int n,std::map<LabelPair,double> bl);
};

// all transition/saddle information here
struct Connection {
	double saddleE, dX, dXmax, Ftol;
	std::vector<double> energies;
	LabelPair SaddleLabels;
	std::list< TransitionSymmetry > self_transitions;
	std::pair<double,double> priornu, nu, fpt, fpT;
	std::pair<double,double> counts;

	Connection();
	Connection(NEBPathway &path);

	void add_path(NEBPathway &path);
	void add_jump(double temperature, double fpt, bool forwards);
	void assimilate(Connection& other);
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & saddleE;
		ar & dX;
		ar & dXmax;
		ar & Ftol;
		ar & energies;
		ar & SaddleLabels;
		ar & self_transitions;

		ar & priornu;
		ar & nu;
		ar & fpt;
		ar & fpT;

		ar & counts;
	};
};

// Canonical labels implicit?
struct StateEdge {
  SymmLabelPair labels; // canon label pair. non_symm sets direction of counts

	StateEdge();
  StateEdge(SymmLabelPair labels_);
	StateEdge(SymmLabelPair labels_, SymmLabelPair refc_);

  // All NEB/TAD data, using the pair that was submitted to the NEB
  std::map<SymmLabelPair,Connection> connections;

  // maps any observed connection to a connection in transition_data + symmetry for the mapping.
	// (If the a self edge (isomorphic) this symmetry is a member of the state's point group mod shift)

	// slight duplicate but important for isomorphic remapping
	std::map<SymmLabelPair,std::pair<SymmLabelPair,PointShiftSymmetry>> self_edge_map;

	std::set<SymmLabelPair> requestedNEBS,completedNEBS,pendingNEBS; // requested / completed / pending NEBS
	std::set<SymmLabelPair> requestedPMS,completedPMS,pendingPMS; // requested / completed / pending PairMaps
	std::map<SymmLabelPair,std::set<SymmLabelPair>> PMtestedPairs;
  // TODO: unknown rate and splicing, basically copying old stuff. Using symmetry to compress
  std::string info_str(bool seg=false);

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & labels;
		ar & connections;
    ar & self_edge_map;
  	ar & requestedNEBS;
		ar & completedNEBS;
		ar & pendingNEBS;
		ar & requestedPMS;
		ar & completedPMS;
		ar & pendingPMS;
		ar & PMtestedPairs;
	};

};


// This is similar to TADStateStatistics- one for each canonical label
struct StateVertex {

  LabelPair reference_label; // starting state for simulations. In principle, all observed transitions will be from here
	double energy,target_state_time;
	int id,clusters;
	std::array<double,3> position;
	// at all the different temperatures. good for reconstruction / Bayes / AHTST...
	unsigned duration, overhead, trials;

	std::map<int,double> elapsedTime;

	double state_time(double targetT);

	std::set<PointShiftSymmetry> self_symmetries; // list of symmetries, added from NEB;

	std::map<Label,PointShiftSymmetry> state_isomorphisms; // transformations from reference label to key

  std::map<LabelPair,double> BasinLabels; // currently only recorded, not implemented

  std::set<Label> connections; // canonical labels to find edge (self,connection)

  std::vector<double> unknown_rate; // at simulation temperatures;
	StateVertex();
  StateVertex(LabelPair ref,double energy_,int id_);
	std::string info_str(double targetT, double eshift, double ku);
  // TEST THIS
  void update(Label lab, double energy_,bool force);
	void update(Label lab, double energy_,int clust, std::array<double,3> pos,bool force);
	void update(Label lab, double energy_,int clust, std::array<double,3> pos, std::set<PointShiftSymmetry> ss,bool force);
	void add_segment(TADSegment &seg);

  template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & reference_label;
		ar & BasinLabels;
		ar & self_symmetries;
		ar & state_isomorphisms;
		ar & energy;
		ar & connections;
    ar & target_state_time;
		ar & elapsedTime;
		ar & overhead;
		ar & trials;
    ar & duration;
		ar & unknown_rate;
		ar & id;
		ar & clusters;
		ar & position;

	};
};


// holds the target rate and all TAD rates (as multiArrhenius) for a StateEdge. Building block of MarkovModel.
struct Rate {
	std::pair<double,double> target_k_fp; // some redundancy but fine..
	std::vector< std::pair<double,double> > tad_k_fp;
	void reset(size_t tadTsize);
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & target_k_fp;
		ar & tad_k_fp;
	};
};

// holds the target rate and all TAD rates (as multiArrhenius) for a StateEdge. Building block of MarkovModel.
struct UnknownRate {
	double unknown_rate; // some redundancy but fine..
	double unknown_variance; // some redundancy but fine..
	double min_rate;
	double observed_rate;

	double optimal_temperature;
	int optimal_temperature_index;
	double optimal_rate;
	double optimal_gradient;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & unknown_rate;
		ar & min_rate;
		ar & unknown_variance;
		ar & optimal_temperature_index;
		ar & optimal_temperature;
		ar & optimal_rate;
		ar & observed_rate;
		ar & optimal_gradient;
	};
};

class TammberModel {
public:
	TammberModel();

  TammberModel(double targetT_, double minT_, double maxT_, int nT_, int pfCT=2, double mB=0.05);

  void parametrize(double targetT_, double minT_, double maxT_, int nT_, int pfCT=2, double mB=0.05);

	void initialize(boost::property_tree::ptree &config,bool restart=false);

	int newIndex();
  // to be overloaded
  void add_vertex(LabelPair labels, double energy,bool force=false);
	void add_vertex(LabelPair labels, double energy,int clusters,std::array<double,3> pos,bool force=false);
	void add_vertex(LabelPair labels, double energy,int clusters,std::array<double,3> pos, std::set<PointShiftSymmetry> ss,bool force=false);

  void add_edge(SymmLabelPair el);
  void add_edge(SymmLabelPair el,SymmLabelPair tl);
	void add_edge(Transition tt);

	//Transition transition(Transition &trans);

  // from TAD
  void add_segment(TADSegment &seg);
	bool fault_check(TADSegment &seg);

	// from NEB
	void add_pathway(NEBPathway &path);
	void add_symmetries(NEBPathway &path);
	bool add_transitionMaps(NEBPathway &path);
	bool add_duplicate(NEBPathway &path);


	std::string info_str(bool seg=false);

	void write_model(std::string mmfile);

	std::list<std::pair<SymmLabelPair,std::pair< std::array<double,6>,PointShiftSymmetry >>> modelEdgeParams(SymmLabelPair el);

	std::map<Label,std::pair<Label,std::set<Label>>> listStates();

	void generateNEBs(std::list<NEBjob> &jobs, int nMax);

	void generateTADs(std::list<TADjob> &jobs, int nMax,bool screen);

	void predict(std::map<Label,std::pair<double,double>> &weights,bool screen);

	void calculate_rates(SymmLabelPair el, Rate &kf, Rate &kb);// canonical labels

	bool allow_allocation(Label lab);

	void unknown_rate(Label lab, UnknownRate &ku); // canonical label

	void bayes_ku_kuvar(std::vector<std::pair<double,double>> &k_fp, double &ku, double &kuvar, double &min_k, double &tot_k, double _time); // Bayes/Poisson evaluation

	void deleteVertex(Label deleteVertex);

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & min_barrier;
		ar & prefactorCountThresh;
		ar & max_ku;
		ar & min_ku;
		ar & max_k;
		ar & HashCost;
		ar & NEBCost;
		ar & targetT;
		ar & targetB;
		ar & valid_time;
		ar & valid_time_sd;
		ar & max_id;
		ar & RhoInitFlavor;
		ar & PredictionSize;
		ar & tadT;
		ar & InitialLabels;

		ar & transitionMap;
		ar & basinMap;
		ar & StateVertices;
		ar & StateEdges;

		ar & Rates;
		ar & SelfRates;
		ar & UnknownRates;

		//ar & pendingTADS;

	};

protected:
  // Bayes parameters
  double min_barrier, prefactorCountThresh, max_ku, min_ku, max_k, HashCost, NEBCost;
  // TAD parameters
  double  targetT, targetB, valid_time, valid_time_sd,pNEB_Prior,DephaseThresh;
	int max_id, RhoInitFlavor, PredictionSize, AllocScheme;
	int ClusterThresh;
  std::vector<double> tadT;
	LabelPair InitialLabels;
	bool safe_opt,sim_conn,include_shallow_states;
	std::string LatticeConstant, PrimitiveUnitCell, UnitCell, SuperCell;

	std::map<Transition,Transition> transitionMap; // not symmetric ?

	std::map<Label,Label> basinMap; // noncanonical -> canonical NOT IMPLEMENTED YET

	std::map<Label,StateVertex> StateVertices; // map of verticies

	std::map<SymmLabelPair, StateEdge> StateEdges; // map of edges

	std::map<Transition,std::list<TADSegment>> pendingTADS; // unclassified MD jumps

	std::map<Label, std::map<Label,Rate> > Rates, SelfRates; // map of rates for MM
	std::map<Label, UnknownRate> UnknownRates; // map of rates for MM

	SpMatCM LowTRT, LowTR;

	#ifdef USE_BOOST_LOG
	boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
	#endif

};

#endif

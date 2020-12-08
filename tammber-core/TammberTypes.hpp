#ifndef TammberTypes_hpp
#define TammberTypes_hpp

#include "Types.hpp"

// Tammber specific types
typedef std::pair<Label,Label> LabelPair;

struct Transition {
	LabelPair first,second;
	Transition();
	Transition(LabelPair f,LabelPair s);
	Transition rev();
	std::string print_str();

	// comparitors for std::set, std::map etc
	bool operator==( const Transition& other) const {
		bool f = bool((first.first==other.first.first) && (first.second==other.first.second) );
		bool s = bool((second.first==other.second.first) && (second.second==other.second.second) );
		return f && s;
	};

	bool operator!=( const Transition& other) const {
		bool f = bool((first.first!=other.first.first) || (first.second!=other.first.second) );
		bool s = bool((second.first!=other.second.first) || (second.second!=other.second.second) );
		return f || s;
	};


	bool operator<( const Transition& other) const {
		if(first.first!=other.first.first) return bool(first.first<other.first.first);
		if(second.first!=other.second.first) return bool(second.first<other.second.first);
		if(first.second!=other.first.second) return bool(first.second<other.first.second);
		return bool(second.second<other.second.second);
	};


	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & first;
		ar & second;
	};
};

struct PointShiftSymmetry {
	bool valid; // true when initialized
	int operation;
	std::array<double,NDIM> shift;//,dcom;
	std::array<double,NDIM*NDIM> matrix;
	std::map<int,int> map; // optional

	// comparitors for std::set, std::map etc
	bool operator==( const PointShiftSymmetry& other) const {
		return bool(operation==other.operation);
	};

	bool operator<( const PointShiftSymmetry& other) const {
		return bool(operation<other.operation);
	};

	PointShiftSymmetry();

	/* fill fM matrix which is a cubic symmetry group implementation */
	void transform_matrix(std::array<double,NDIM*NDIM> &fM, int op);

	// could do this analytically
	int find_inverse_op(int op_);
	int find_compound_op(int op1_,int op2_);
	// combine two transforms....
	// I->I2 : A => A.I = I2 => if I == C.I, we need A.C = B => C_I = iA.B
	// F->F2 : B => B.F = F2 => if F == C.F, we need B.C = A => C_F = iB.A = iC_I
	// if F == C.F
	// Need to find B.C = A, C.B = A, C.A = B, B.C=A.....
	// TODO add shift!

	void combine_transform(int one, int two, std::set<int> &onesym, std::set<int> &twosym);


	PointShiftSymmetry inverse();
	// return effect of other*this:
	// this.X == P.X + d == P,d
	// other.Y = Q.Y+e == Q,e
	// other.this.X = Q.P.X + Q.d+e == QP, Qd+e => P -> QP, d->Qd+e
	PointShiftSymmetry compound(PointShiftSymmetry &other);

	std::string info_str() const {
		std::string res = "\nOperation: "+std::to_string(operation)+"\n";
		for(int i=0;i<NDIM;i++) {
			if(i!=1) res += "\t";
			else res += "Matrix:\t";
			for(int j=0;j<NDIM;j++) {
				if(int(matrix[NDIM*i+j])>=0) res+="+";
				res += std::to_string(int(matrix[NDIM*i+j]))+" ";
			}
			res += "\n";
		}
		double smag=0.0;
		res+="shift: [ ";
		for(int i=0;i<NDIM;i++) {
			res+= std::to_string(shift[i])+" ";
			smag += shift[i]*shift[i];
		}
		res+="] |shift| = "+std::to_string(sqrt(smag))+"\n";

		/*
		smag=0.0;
		res+="dcom: [ ";
		for(int i=0;i<NDIM;i++) {
			res+= std::to_string(dcom[i])+" ";
			smag += dcom[i]*dcom[i];
		}
		res+="] |dcom| = "+std::to_string(sqrt(smag))+"\n";
		*/

		return res;
	};

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
		ar & operation;
		ar & shift;
		//ar & dcom;
		ar & matrix;
		ar & map;
		ar & valid;
	};
};

typedef std::pair<Transition,PointShiftSymmetry> TransitionSymmetry;

// Return struct of TASK_SEGMENT
struct TADSegment {
	double initialE,finalE,temperature,elapsedTime;
	std::array<double,3> initialPosition,finalPosition;
	int initialClusters,finalClusters;
	bool dephased,OutOfBasin;
	std::unordered_map<LabelPair,double,boost::hash<LabelPair>> BasinLabels;
	LabelPair InitialLabels;
	Transition transition;
	unsigned int trials,duration,overhead;
	// Shouldn't need multiple....
	//std::set<PointShiftSymmetry> ops; // if transition == isomorphism
	TADSegment();
	std::string info_str();
	std::string submit_info_str();
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & initialE;
		ar & finalE;
		ar & temperature;
		ar & elapsedTime;
		ar & dephased;
		ar & BasinLabels;
		ar & InitialLabels;
		ar & transition;
		ar & trials;
		ar & duration;
		ar & overhead;
		ar & OutOfBasin;
		ar & initialPosition;
		ar & finalPosition;
		ar & initialClusters;
		ar & finalClusters;
		//ar & ops;
	};
};

// Return struct of TASK_NEB
struct NEBPathway {
	LabelPair InitialLabels,FinalLabels,SaddleLabels,MMInitialLabels,MMFinalLabels;
	std::list<Transition> FoundTransitions;
	std::set<PointShiftSymmetry> InitialSymmetries, FinalSymmetries;

	bool duplicate,valid,mismatch,pairmap;
	double initialE,finalE,saddleE;

	std::vector<double> energies;
	std::pair<double,double> priornu;

	std::list< TransitionSymmetry > self_transitions; // when isomorphic
	std::list< TransitionSymmetry > equivalent_transitions; // after spacemap
	std::list< TransitionSymmetry > equivalent_states;
	std::list<LabelPair> compared_transitions;

	double dX,Ftol,dXmax;
	NEBPathway();
	std::string info_str(bool full=false);
	std::string submit_info_str();
	//void reverse();
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & InitialLabels;
		ar & FinalLabels;
		ar & SaddleLabels;
		ar & MMInitialLabels;
		ar & MMFinalLabels;
		ar & InitialSymmetries;
		ar & FinalSymmetries;
		ar & FoundTransitions;

		ar & self_transitions;
		ar & equivalent_transitions;
		ar & equivalent_states;
		ar & compared_transitions;

		ar & priornu;
		ar & dX;
		ar & dXmax;
		ar & Ftol;
		ar & initialE;
		ar & finalE;
		ar & saddleE;

		ar & mismatch;
		ar & duplicate;
		ar & pairmap;
		ar & valid;
		ar & energies;
	};
};

#endif

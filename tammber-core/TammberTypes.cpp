#include "TammberTypes.hpp"

// Tammber specific types
Transition::Transition(){
	/* do nothing */
};
Transition::Transition(LabelPair f,LabelPair s) {
	first = f;
	second = s;
};
Transition Transition::rev(){
	Transition res;
	res.second = first;
	res.first = second;
	return res;
};
std::string Transition::print_str() {
	std::string res="\t";
	res+=std::to_string(first.first)+","+std::to_string(first.second)+" -> ";
	res+=std::to_string(second.first)+","+std::to_string(second.second)+"\n";
	return res;
};

PointShiftSymmetry::PointShiftSymmetry() {
	for(int i=0;i<NDIM*NDIM;i++) matrix[i]=0.0;
	for(int i=0;i<NDIM;i++) {
		matrix[NDIM*i+i] = 1.0;
		shift[i] = 0.0;
	}
	operation=8+40*(NDIM-2);
	valid=false;
};

/* fill fM matrix which is a cubic symmetry group implementation */
void PointShiftSymmetry::transform_matrix(std::array<double,NDIM*NDIM> &fM, int op) {

	// M = Identity
	std::array<double,NDIM*NDIM> M;
	for(int i=0;i<NDIM*NDIM;i++) fM[i] = 0.0;
	for(int i=0;i<NDIM*NDIM;i++) M[i] = 0.0;
	for(int i=0;i<NDIM;i++) M[i*NDIM+i] = 1.0;

	#if NDIM==3
	/*
		48 SYMMETRIES OF THE CUBE = 3!=6 permutations * 2^3=8 reflections
		operation = 6R + P,  R E [0,7],  P E [0,5]
		R = 0 nothing
		R = 1-3 single reflection x,y,z
		R = 4-6 double reflection xy yz zx NOTE ORDER
		R = 7 triple reflection
		P = 0 nothing
		P = 1-2 cyclic permutation: nothing, xyz -> yzx,  xyz -> zxy
		P = 3-5 x<->y then cyclic permutation
	*/
	if(op>=48) op=0;
	unsigned int R = op/6;
	unsigned int P = op%6;
	// Reflections
	// Diagonal elements 0,4,8
	// R = [1-3] : M[(R-1)*4] = -1 														==>    M[4*((R-1+i)%3)] = -1 for i = 0
	// R = [4-6] : M[(R-4)*4] = -1 and M[4*(R-4+1)%3] = -1    ==>    M[4*((R-1+i)%3)] = -1 for i = 0,1
	// R = 7 : M[0] = -1 and M[4] = -1 and M[8] = -1					==>    M[4*((R-1+i)%3)] = -1 for i = 0,1,2
	int ii = (R>0)*(1+(R-1)/3); // R->ii : 0->0, [1-3]->1, [4-6]->2, 7->3
	for(int i=0;i<ii; i++) M[4*((R-1+i)%3)] = -1.;
	// Permutations
	for(int i=0;i<3;i++) {
		// i->ii     P<3: 012->012     else:   012->102   (swap x<->y if P >= 3)
		// P<3:   P/3 =0, ii = i
		// P>=3:  P/3 =1, ii = 1-i%2 + i/2  [  i<2:   ii = 1-i  (01->10)      i==2:   ii = 1+1 = 2 ]
		ii = (1-i%2+i/2)*(P/3)+i*(1-(P/3)); // swap x<->y if P >= 3
 		for(int j=0;j<3;j++) fM[3*ii+j] = M[ 3*((i+P)%3) + j]; // ii swaps xy, (i+P)%3 performs cyclic permutation P%3 times
	}
	#else
	/*
		8 SYMMETRIES OF THE SQUARE = 2!=2 permutations * 2^2=4 reflections
		operation = 2R + P,  R E [0,3],  P E [0,1]
		R = 0,1,2,3 : nothing, single reflection x,y, double reflection
		P = 0,1 nothing, permutation: xy -> yx
	*/
	if(op>=8) op=0;
	unsigned int R = op/2;
	unsigned int P = op%2;
	fM[0][0] = ((R==0)-(R%2))*(P==0);
	fM[1][0] = ((R==0)-(R%2))*(P==1);
	fM[0][1] = ((R==0)-(R>1))*(P==1);
	fM[1][1] = ((R==0)-(R>1))*(P==0);
	#endif
};

// could do this analytically
int PointShiftSymmetry::find_inverse_op(int op_) {
	std::array<double,NDIM*NDIM> M,N,P;
	double cost;
	transform_matrix(M,op_);
	for(int op=0;op<8+40*(NDIM-2);op++) {
		for(int i=0;i<NDIM*NDIM;i++) P[i]=0.0;

		transform_matrix(N,op);
		for(int i=0;i<NDIM;i++) for(int j=0;j<3;j++) for(int k=0;k<NDIM;k++)
			P[i*NDIM+j] += M[i*NDIM+k] * N[k*NDIM+j];
		cost = 0.0;
		for(int i=0;i<NDIM;i++) for(int j=0;j<3;j++)
			cost += (P[i*NDIM+j]-double(i==j))*(P[i*NDIM+j]-double(i==j));
		if(cost<0.001) return op;
	}
	return -1;
};

int PointShiftSymmetry::find_compound_op(int op1_,int op2_) {
	std::array<double,NDIM*NDIM> M,N,P;
	double cost;
	transform_matrix(M,op1_);
	transform_matrix(N,op2_);
	// form compound
	for(int i=0;i<NDIM;i++) for(int j=0;j<3;j++) {
		P[i*NDIM+j] = 0.0;
		for(int k=0;k<NDIM;k++) P[i*NDIM+j] += M[i*NDIM+k] * N[k*NDIM+j];
	}
	// which transform is it?
	for(int op=0;op<8+40*(NDIM-2);op++) {
		transform_matrix(M,op);
		cost=0.0;
		for(int i=0;i<NDIM*NDIM;i++) cost+=(M[i]-P[i])*(M[i]-P[i]);
		if(cost<0.001) return op;
	}
	return -1;
};

// combine two transforms....
// I->I2 : A => A.I = I2 => if I == C.I, we need A.C = B => C_I = iA.B
// F->F2 : B => B.F = F2 => if F == C.F, we need B.C = A => C_F = iB.A = iC_I
// if F == C.F
// Need to find B.C = A, C.B = A, C.A = B, B.C=A.....
// TODO add shift!

void PointShiftSymmetry::combine_transform(int one, int two, std::set<int> &onesym, std::set<int> &twosym) {
	std::array<double,NDIM*NDIM> A,B,C,D;
	transform_matrix(A,one);
	transform_matrix(B,two);
	bool match=false;
	for(int op=0;op<8+40*(NDIM-2);op++) {
		transform_matrix(C,op);
		// D = A.C
		for(int i=0;i<NDIM*NDIM;i++) D[i] = 0.0;
		for(int i=0;i<NDIM;i++)
			for(int j=0;j<NDIM;j++)
				for(int k=0;k<NDIM;k++) D[3*i+j] += A[NDIM*i+k] * C[NDIM*k+j];
		// D == B?
		match=true;
		for(int i=0;i<NDIM*NDIM;i++) if((D[i]-B[i])*(D[i]-B[i])>0.000000001) match = false;
		if(match) onesym.insert(op);

		// D = B.C
		for(int i=0;i<NDIM*NDIM;i++) D[i] = 0.0;
		for(int i=0;i<NDIM;i++)
			for(int j=0;j<NDIM;j++)
				for(int k=0;k<NDIM;k++) D[3*i+j] += B[NDIM*i+k] * C[NDIM*k+j];
		// D == A?
		match=true;
		for(int i=0;i<NDIM*NDIM;i++) if((D[i]-A[i])*(D[i]-A[i])>0.000000001) match = false;
		if(match) twosym.insert(op);
	}
};



PointShiftSymmetry PointShiftSymmetry::inverse() {
	PointShiftSymmetry res;
	res.valid=valid;
	res.operation=find_inverse_op(operation);
	std::array<double,NDIM> ts;
	for(int i=0;i<NDIM;i++) res.shift[i] = 0.0;
	transform_matrix(res.matrix,res.operation);
	for(int i=0;i<NDIM;i++) for(int j=0;j<NDIM;j++) res.shift[i] += -res.matrix[NDIM*i+j]*shift[i];
	if(map.size()>0) for(auto amap:map) res.map[amap.second] = amap.first;
	return res;
};

// return effect of other*this:
// this.X == P.X + d == P,d
// other.Y = Q.Y+e == Q,e
// other.this.X = Q.P.X + Q.d+e == QP, Qd+e => P -> QP, d->Qd+e
PointShiftSymmetry PointShiftSymmetry::compound(PointShiftSymmetry &other) {
	PointShiftSymmetry res;
	res.valid=valid;
	// P -> QP
	res.operation = find_compound_op(other.operation,operation);
	transform_matrix(res.matrix,res.operation);
	// d -> Qd+e
	for(int i=0;i<NDIM;i++) {
		for(int j=0;j<NDIM;j++) res.shift[i] += other.matrix[i*NDIM+j]*shift[j];
		res.shift[i] += other.shift[i];
	}
	if (map.size()>0 && other.map.size()>0) {
		res.map = map;
		for(auto amap:map) res.map[amap.first] = other.map[amap.second];
	}
	return res;
};

TADSegment::TADSegment() {
	dephased=false;
	OutOfBasin=false;
	trials=0;
	duration=0;
	overhead=0;
	elapsedTime=0.0;
	initialE=0.0;
	finalE=0.0;
	InitialLabels = std::make_pair(0,0);
	transition.first = InitialLabels;
	transition.second = InitialLabels;
	initialClusters = 0;
	finalClusters = 0;
	initialPosition = {0.,0.,0.};
	finalPosition = {0.,0.,0.};
};

std::string TADSegment::info_str() {
	std::string str="TADSegment:\n";
	str += "InitialLabels: "+std::to_string(InitialLabels.first)+","+std::to_string(InitialLabels.second)+"\n";
	str += " Clusters: "+std::to_string(initialClusters)+" Position: [";
	str += std::to_string(initialPosition[0])+",";
	str += std::to_string(initialPosition[1])+",";
	str += std::to_string(initialPosition[2])+"],";
	str += "dephased: "+std::to_string(dephased)+"\n";
	str += "OutOfBasin: "+std::to_string(OutOfBasin)+"\n";
	str += "transition: "+std::to_string(transition.first.first)+","+std::to_string(transition.first.second)+" -> "+std::to_string(transition.second.first)+","+std::to_string(transition.second.second)+"\n";
	str += "energy: "+std::to_string(initialE)+" -> "+std::to_string(finalE)+" ("+std::to_string(finalE-initialE)+")\n";
	if(OutOfBasin){
		str += " Final Clusters: "+std::to_string(finalClusters)+" Position: [";
		str += std::to_string(finalPosition[0])+",";
		str += std::to_string(finalPosition[1])+",";
		str += std::to_string(finalPosition[2])+"],";
	}
	str += "Temperature: "+std::to_string(temperature)+"\n";
	str += "elapsedTime, trials, duration, overhead: "+std::to_string(elapsedTime)+" "+std::to_string(trials)+", "+std::to_string(duration)+", "+std::to_string(overhead)+"\n";
	str += "BasinLabels: \n";
	for(auto bb: BasinLabels)
		str +=std::to_string(bb.first.first)+","+std::to_string(bb.first.second)+": "+std::to_string(bb.second)+"eV\n";
	return str;
};

std::string TADSegment::submit_info_str() {
	std::string str="TADSegment:\n";
	str += "InitialLabels: "+std::to_string(InitialLabels.first)+","+std::to_string(InitialLabels.second)+"\n";
	str += "Temperature: "+std::to_string(temperature)+"\n";
	str += "BasinLabel Size:"+std::to_string(BasinLabels.size())+" \n";
	for(auto bb: BasinLabels)
		str +=std::to_string(bb.first.first)+","+std::to_string(bb.first.second)+": "+std::to_string(bb.second)+"eV\n";
	return str;
};


NEBPathway::NEBPathway() {
	LOGGER("NEBPathway::NEBPathway")
	duplicate=false;
	dX = 10.0*MSD_THRESH;
	valid=true;
	mismatch=false;
	Ftol=100.0;

};

std::string NEBPathway::info_str(bool full){
	LOGGER("NEBPathway::info_str")
	std::string res="NEBPathway:\n";
	res += "("+std::to_string(InitialLabels.first)+","+std::to_string(InitialLabels.second)+") -> ";
	res += "("+std::to_string(SaddleLabels.first)+","+std::to_string(SaddleLabels.second)+") -> ";
	res += "("+std::to_string(FinalLabels.first)+","+std::to_string(FinalLabels.second)+")\n";
	res += "Duplicate: "+std::to_string(duplicate)+"\n";
	if (mismatch) {
		res += "MISMATCHED PATH! AFTER MIN: ";
		res += "("+std::to_string(MMInitialLabels.first)+","+std::to_string(MMInitialLabels.second)+") -> ";
		res += "("+std::to_string(MMFinalLabels.first)+","+std::to_string(MMFinalLabels.second)+")\n";
	}
	if(!duplicate) {
		res += "Found Transitions:\n";
		for(auto ft:FoundTransitions) {
			res += "("+std::to_string(ft.first.first)+","+std::to_string(ft.first.second)+")";
			res += " -> ("+std::to_string(ft.second.first)+","+std::to_string(ft.second.second)+")\n";
		}
	}
	if(full) {
		res += "Self Symmetries for "+std::to_string(InitialLabels.second)+":";
		for(auto ss: InitialSymmetries) res += ss.info_str();
		if(FinalLabels.first!=InitialLabels.first) {
			res += "Self Symmetries for "+std::to_string(FinalLabels.second)+":\n";
			for(auto ss: FinalSymmetries) res += ss.info_str();
		} else {
			res += "Self Transition Symmetry:\n";
			for(auto lss: self_transitions) res += lss.second.info_str();
		}
		if(equivalent_transitions.size()>0) {
			res += "Equivalent Transitions:\n";
			for(auto lss: equivalent_transitions) {
				auto lab = lss.first.first;
				res += "("+std::to_string(lab.first)+","+std::to_string(lab.second)+") -> ";
				lab = lss.first.second;
				res += "("+std::to_string(lab.first)+","+std::to_string(lab.second)+")\n";
				res += lss.second.info_str();
			}
		}
	}
	if(!duplicate) {
		res += "NEB results:\n";
		res += "Energy: 0.0 -> "+std::to_string(saddleE-initialE)+" -> "+std::to_string(finalE-initialE)+" (+"+std::to_string(initialE)+")\n";
		res += "Prefactors: "+std::to_string(priornu.first)+","+std::to_string(priornu.second)+"\n";
		res += "dX, dXmax, Ftol: "+std::to_string(dX)+","+std::to_string(dXmax)+","+std::to_string(Ftol)+"\n";
	}
	return res;
};

std::string NEBPathway::submit_info_str(){
	LOGGER("NEBPathway::submit_info_str")
	std::string res="NEBPathway:\n";
	res += "("+std::to_string(InitialLabels.first)+","+std::to_string(InitialLabels.second)+") -> ";
	res += "("+std::to_string(FinalLabels.first)+","+std::to_string(FinalLabels.second)+")\n";
	return res;
};


/*
// this should be avoided
void NEBPathway::reverse() {

	LabelPair temp_lp = InitialLabels;
	InitialLabels = FinalLabels;
	FinalLabels = temp_lp;

	std::set<PointShiftSymmetry> temp_psl = InitialSymmetries;
	InitialSymmetries = FinalSymmetries;
	FinalSymmetries = temp_psl;

	double tempE = initialE;
	initialE = finalE;
	finalE = tempE;

	for(auto it=FoundTransitions.begin();it!=FoundTransitions.end();it++) {
		temp_lp = it->first;
		it->first = it->second;
		it->second = temp_lp;
	}

};*/

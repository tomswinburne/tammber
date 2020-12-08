
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <map>
#include <list>
using Eigen::MatrixXd;
using Eigen::VectorXd;


class TICAEstimator {
public:

double dt;
int lag;
int nDescriptors;
int nReplicas;

std::map<int, std::list<VectorXd> > rData;
VectorXd rsum;
int nsum;
MatrixXd rcross;
MatrixXd rcross0;
MatrixXd c00;
MatrixXd c0t;
MatrixXd c0tc;
MatrixXd transform;
MatrixXd w;
MatrixXd v;

int ncross;
int ncross0;


TICAEstimator(double dt_, int lag_, int nDescriptors_, int nReplicas_){
	dt=dt_;
	lag=lag_;
	nDescriptors=nDescriptors_;
	nReplicas=nReplicas_;


	rsum=VectorXd::Zero(nDescriptors);
	nsum=0;
	rcross=MatrixXd::Zero(nDescriptors, nDescriptors);
	ncross=0;

	c00=MatrixXd::Zero(nDescriptors, nDescriptors);
	c0t=MatrixXd::Zero(nDescriptors, nDescriptors);
	c0tc=MatrixXd::Zero(nDescriptors, nDescriptors);
};

void appendData(int index, VectorXd &data){
	rData[index].push_back(data);

	rsum+=data;
	nsum+=1;

	rcross0+=data*data.transpose();
	ncross0+=1;
	if(rData[index].size()>lag) {
		rcross+=rData[index].front()*rData[index].back().transpose();
		rcross+=rData[index].back()*rData[index].front().transpose();
		ncross+=2;
	}
	while(rData[index].size()>=lag) {
		rData[index].pop_front();
	}
};

void update(){
	double th=1e-12;

	VectorXd mu=rsum/nsum;
	//form c00
	c00=rcross0/ncross0;
	c00-=mu*mu.transpose();

	//diagonalize
	Eigen::SelfAdjointEigenSolver<MatrixXd> solver0(c00);


	VectorXd scale=VectorXd::Zero(nDescriptors);

	//build the whitening matrix
	for(int i=0; i<nDescriptors; i++) {
		if(solver0.eigenvalues()[i]>th) {
			scale(i)=1/sqrt(solver0.eigenvalues()[i]);
		}
		else{
			scale(i)=1;
		}
	}

	transform=scale.asDiagonal()*solver0.eigenvectors();

	//transform the means
	mu=transform*mu;

	//builds the correlation matrix
	c0tc=rcross/ncross;
	c0t=transform*c0tc*transform.transpose() - mu*mu.transpose();

	//diagonalize
	Eigen::SelfAdjointEigenSolver<MatrixXd> solver(c0t);
	w=solver.eigenvalues();
	v=solver.eigenvectors();

};

double distance(double dTime, VectorXd &ci, VectorXd &cj){
	double texp=dTime/dt;

	VectorXd si=transform*ci;
	VectorXd sj=transform*cj;

	double d=0;
	for(int i=0; i<nDescriptors; i++) {
		double c=v.col(i).transpose()*(si-sj);
		d+=pow(pow(w(i),texp)*c,2 );
	}
	return sqrt(d);

};

};

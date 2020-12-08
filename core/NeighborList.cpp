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


#include "NeighborList.hpp"
#include "BoundaryConditions.hpp"
#include <array>
#include <numeric>
#include <cmath>
#include <set>


std::vector<std::vector<int> >  generateShifts(int nDim,int iMin,int iMax, int increment, std::vector<std::vector<int> > shifts){

	std::vector<std::vector<int> > exp;
	if(shifts.size()==0) {
		shifts.push_back(std::vector<int>());
	}

	if(nDim==0) {
		return shifts;
	}
	else{
		for(int i=0; i<shifts.size(); i++) {
			for(int j=iMin; j<=iMax; j+=increment) {
				std::vector<int> tv=shifts[i];
				tv.push_back(j);
				exp.push_back(tv);
			}
		}
		return generateShifts(nDim-1,iMin,iMax,increment,exp);
	}

};

NeighborList::NeighborList(){
	//Generate a vector of all shifts in NDIM dimensions
	shifts=generateShifts(NDIM,-1,1,1);
};

NeighborList::NeighborList(const NeighborList &source){
	neighborList=source.neighborList;
};

NeighborList::~NeighborList(){
};

NeighborList &NeighborList::operator=(const NeighborList &source){
	neighborList=source.neighborList;
	return *this;
};





void NeighborList::build(AbstractSystem &system, std::map<std::pair<int,int>,double> &cutoffs){

	int nAtoms=system.getNAtoms();


	bc=system.getCell();

	std::array<double,NDIM> rPosition;
	std::array<double,NDIM> sPosition;
	std::array<int,NDIM> sc;
	std::array<int,NDIM> cg;
	std::array<int,NDIM> cg2;
	std::array<int,NDIM> nCells;
	std::array<double,NDIM> cellSizeS;
	std::array<double,NDIM> cellSizeR;

	/*
	   //compute the projectors here.
	   MatrixXd rspe(NDIM,NDIM);
	   MatrixXd cspe(NDIM,NDIM);

	   for(int i=0; i<NDIM; i++) {
	        bc.periodic[i]=system.getPeriodic(i);
	        for(int j=0; j<NDIM; j++) {
	                bc.rsp[i][j]=system.getBox(i,j);
	                rspe(i,j)=system.getBox(i,j);
	        }
	   }
	   cspe=rspe.inverse();
	   for(int i=0; i<NDIM; i++) {
	        for(int j=0; j<NDIM; j++) {
	                bc.csp[i][j]=cspe(i,j);
	        }
	   }
	 */

	{
		//boost::timer::auto_cpu_timer t;
		//compute the length of the edges of the cell
		std::array<double,NDIM> edgeLength;
		for(int i=0; i<NDIM; i++) {
			edgeLength[i]=0;
			for(int j=0; j<NDIM; j++) {
				edgeLength[i]+=bc.rsp[j][i]*bc.rsp[j][i];
			}
			edgeLength[i]=sqrt(edgeLength[i]);
			//std::cout<<"EDGE LENGTH("<<i<<")= "<<edgeLength[i]<<std::endl;
		}


		//identify the maximum cutoff in order to build a grid of cells
		double maxCut=0;
		for(auto it=cutoffs.begin(); it!=cutoffs.end(); it++) {
			maxCut=std::max((*it).second,maxCut);
		}



		std::array<double,NDIM> minCellSize;
		std::array<double,NDIM> bv;
		std::array<double,NDIM> bvO;
		double bvdotb;
		double corrFact;
		//for every basis vector
		for(int i=0; i<NDIM; i++) {
			for(int j=0; j<NDIM; j++) {
				bv[j]=bc.rsp[j][i];
				bvO[j]=bv[j];
			}

			//orthogonalize against all other basis vectors
			for(int j=0; j<NDIM; j++) {
				if(i!=j) {
					bvdotb=0;
					for(int k=0; k<NDIM; k++) {
						bvdotb+=bv[k]*bc.rsp[k][j];
					}
					bvdotb/=pow(edgeLength[j],2);
					for(int k=0; k<NDIM; k++) {
						bv[k]-=bvdotb*bc.rsp[k][j];
					}
					//cout<<i<<" "<<j<<" "<<bvdotb<<endl;
				}
			}
			//correct  to be sure that neighbors can only be in the first layer of cells even for strongly non-cubic cells
			//for that, the perpendicular width of the cell must be larger than minCellSize.
			corrFact=edgeLength[i]/sqrt( std::inner_product(bvO.begin(), bvO.end(), bv.begin(), 0.) );
			//cout<<"CORRECTION FACTOR FOR VECTOR "<<i<<" = "<<corrFact<<endl;
			minCellSize[i]=(maxCut*corrFact);
			//cout<<"MINIMUM CELL SIZE "<<i<<" = "<<minCellSize(i)<<endl;
		}



		//Compute the cell size in every direction
		for(int j=0; j<NDIM; j++) {
			nCells[j]=int(edgeLength[j]/minCellSize[j]);
			nCells[j]=std::max(nCells[j],1);
			cellSizeR[j]=(edgeLength[j]/nCells[j]);
			cellSizeS[j]=1.0/nCells[j];
			//std::cout<<j<<" NCELLS "<<nCells[j]<<" "<<cellSizeR[j]<<std::endl;
		}


		cells.clear();


		//compute the cell indices for each atom
		for(int i=0; i<nAtoms; i++) {
			//Get the s-coordinates in the primary cell
			//gPrimaryCellPosition(origin,periodicity,rsp,csp,position,i,rPosition,sPosition);
			auto p=primaryCellPosition(i, system);
			rPosition=p.first;
			sPosition=p.second;

			//get the cell coordinates in grid-space
			for(int j=0; j<NDIM; j++) {
				sc[j]=floor(double(sPosition[j])/cellSizeS[j]);
			}

			int cellIndex=sc[0]*nCells[1]*nCells[2]+sc[1]*nCells[2]+sc[2];

			//std::cout<<cellIndex<<" "<<sc[0]<<" "<<sc[1]<<" "<<sc[2]<<" "<<sPosition[0]<<" "<<sPosition[1]<<" "<<sPosition[2]<<std::endl;
			cells[cellIndex].push_back(i);
		}

		//std::cout<<"CELLS: "<<cells.size()<<std::endl;
	}

	{

		//boost::timer::auto_cpu_timer t;
		neighborList.clear();


		//Build the neighbor list
		size_t nin=neighborList.size();
		for(int i=0; i<nin; i++) {
			neighborList[i].clear();
			neighborList[i].reserve(70);
		}
		for(size_t i=nin; i<nAtoms; i++) {
			neighborList.push_back(std::vector<int>());
			neighborList[i].reserve(70);
		}

	}



	{
		//boost::timer::auto_cpu_timer t;


		{
			//for each cell
			//for(int iCell=0;iCell<nCellsT;iCell++)
			for(auto itCell1=cells.begin(); itCell1!=cells.end(); itCell1++) {
				int iCell1=itCell1->first;

				if( cells[iCell1].size()>0) {
					//compute the cell coordinate
					int ict=iCell1;
					cg[2]=ict%nCells[2];
					ict-=cg[2];
					ict/=nCells[2];
					cg[1]=ict%nCells[1];
					ict-=cg[1];
					ict/=nCells[1];
					cg[0]=ict;

					//make a list of cells we need to scan
					std::set<int> scannedCells;
					for(int si=0; si<shifts.size(); si++) {
						for(int sj=0; sj<NDIM; sj++) {
							cg2[sj]=cg[sj]+shifts[si][sj];
							//fold the grid coordinates
							cg2[sj]-=nCells[sj]*floor(double(cg2[sj])/nCells[sj]);
						}

						int iCell2=cg2[0]*nCells[1]*nCells[2]+cg2[1]*nCells[2]+cg2[2];
						if(cells[iCell2].size()>0 and iCell1 <=iCell2) {
							scannedCells.insert(iCell2);
						}
					}

					//std::cout<<scannedCells.size()<<" "<<shifts.size()<<std::endl;
					//loop on the neighboring cells
					for(std::set<int>::iterator itCell2=scannedCells.begin(); itCell2!=scannedCells.end(); itCell2++) {
						int iCell2=*itCell2;


						for(auto itAtom1=cells[iCell1].begin(); itAtom1!=cells[iCell1].end(); itAtom1++ ) {
							for(auto itAtom2=cells[iCell2].begin(); itAtom2!=cells[iCell2].end(); itAtom2++ ) {
								double dr=nearestImageDistance(*itAtom1, *itAtom2, system, bc);

								if(dr<cutoffs[std::make_pair(system.getSpecies(*itAtom1),system.getSpecies(*itAtom2))]) {
									if(*itAtom1!=*itAtom2) {
										neighborList[*itAtom1].push_back(*itAtom2);
										if(iCell1!=iCell2) {
											neighborList[*itAtom2].push_back(*itAtom1);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
};



void NeighborList::buildJoint(AbstractSystem &reference, AbstractSystem &system, std::map<std::pair<int,int>,double> &cutoffs){
	int nAtomsRef=reference.getNAtoms();
	int nAtoms=system.getNAtoms();

	bc=system.getCell();


	//compute the length of the edges of the cell
	std::array<double,NDIM> edgeLength;
	for(int i=0; i<NDIM; i++) {
		edgeLength[i]=0;
		for(int j=0; j<NDIM; j++) {
			edgeLength[i]+=bc.rsp[j][i]*bc.rsp[j][i];
		}
		edgeLength[i]=sqrt(edgeLength[i]);
		//cout<<"EDGE LENGTH("<<i<<")= "<<edgeLength(i)<<endl;
	}


	//identify the maximum cutoff in order to build a grid of cells
	double maxCut=0;
	for(auto it=cutoffs.begin(); it!=cutoffs.end(); it++) {
		maxCut=std::max((*it).second,maxCut);
	}



	std::array<double,NDIM> minCellSize;
	std::array<double,NDIM> bv;
	std::array<double,NDIM> bvO;
	double bvdotb;
	double corrFact;
	//for every basis vector
	for(int i=0; i<NDIM; i++) {
		for(int j=0; j<NDIM; j++) {
			bv[j]=bc.rsp[j][i];
			bvO[j]=bv[j];
		}

		//orthogonalize against all other basis vectors
		for(int j=0; j<NDIM; j++) {
			if(i!=j) {
				bvdotb=0;
				for(int k=0; k<NDIM; k++) {
					bvdotb+=bv[k]*bc.rsp[k][j];
				}
				bvdotb/=pow(edgeLength[j],2);
				for(int k=0; k<NDIM; k++) {
					bv[k]-=bvdotb*bc.rsp[k][j];
				}
				//cout<<i<<" "<<j<<" "<<bvdotb<<endl;
			}
		}
		//correct  to be sure that neighbors can only be in the first layer of cells even for strongly non-cubic cells
		//for that, the perpendicular width of the cell must be larger than minCellSize.
		corrFact=edgeLength[i]/sqrt( std::inner_product(bvO.begin(), bvO.end(), bv.begin(), 0.) );
		//cout<<"CORRECTION FACTOR FOR VECTOR "<<i<<" = "<<corrFact<<endl;
		minCellSize[i]=(maxCut*corrFact);
		//cout<<"MINIMUM CELL SIZE "<<i<<" = "<<minCellSize(i)<<endl;
	}

	std::array<int,NDIM> nCells;
	std::array<double,NDIM> cellSizeS;
	std::array<double,NDIM> cellSizeR;

	//Compute the cell size in every direction
	for(int j=0; j<NDIM; j++) {
		nCells[j]=int(edgeLength[j]/minCellSize[j]);
		nCells[j]=std::max(nCells[j],1);
		cellSizeR[j]=(edgeLength[j]/nCells[j]);
		cellSizeS[j]=1.0/nCells[j];
		//cout<<j<<" NCELLS "<<nCells(j)<<" "<<cellSizeR(j)<<endl;
	}


	std::array<double,NDIM> rPosition;
	std::array<double,NDIM> sPosition;
	std::array<int,NDIM> sc;
	std::array<int,NDIM> cg;
	std::array<int,NDIM> cg2;


	cells.clear();


	std::map<int, std::list<int> > cellsRef;
	std::map<int, std::list<int> > cellsSys;

	//compute the cell indices for each atom
	for(int i=0; i<nAtomsRef; i++) {
		//Get the s-coordinates in the primary cell
		//gPrimaryCellPosition(origin,periodicity,rsp,csp,position,i,rPosition,sPosition);
		auto p=primaryCellPosition(i, reference);
		rPosition=p.first;
		sPosition=p.second;

		//get the cell coordinates in grid-space
		for(int j=0; j<NDIM; j++) {
			sc[j]=floor(double(sPosition[j])/cellSizeS[j]);
		}

		int cellIndex=sc[0]*nCells[1]*nCells[2]+sc[1]*nCells[2]+sc[2];
		cellsRef[cellIndex].push_back(i);
	}

	for(int i=0; i<nAtoms; i++) {
		//Get the s-coordinates in the primary cell
		//gPrimaryCellPosition(origin,periodicity,rsp,csp,position,i,rPosition,sPosition);
		auto p=primaryCellPosition(i, system);
		rPosition=p.first;
		sPosition=p.second;

		//get the cell coordinates in grid-space
		for(int j=0; j<NDIM; j++) {
			sc[j]=floor(double(sPosition[j])/cellSizeS[j]);
		}

		int cellIndex=sc[0]*nCells[1]*nCells[2]+sc[1]*nCells[2]+sc[2];
		cellsSys[cellIndex].push_back(i);
	}




	neighborList.clear();


	//Build the neighbor list
	size_t nin=neighborList.size();
	for(int i=0; i<nin; i++) {
		neighborList[i].clear();
		neighborList[i].reserve(70);
	}
	for(size_t i=nin; i<nAtoms; i++) {
		neighborList.push_back(std::vector<int>());
		neighborList[i].reserve(70);
	}



	{
		//for each cell
		//for(int iCell=0;iCell<nCellsT;iCell++)
		//loop on all cells from system
		for(auto itCell1=cellsSys.begin(); itCell1!=cellsSys.end(); itCell1++) {
			int iCell1=itCell1->first;

			if( cellsSys[iCell1].size()>0) {
				//compute the cell coordinate
				int ict=iCell1;
				cg[2]=ict%nCells[2];
				ict-=cg[2];
				ict/=nCells[2];
				cg[1]=ict%nCells[1];
				ict-=cg[1];
				ict/=nCells[1];
				cg[0]=ict;

				//make a list of cells we need to scan
				std::set<int> scannedCells;
				for(int si=0; si<shifts.size(); si++) {
					for(int sj=0; sj<NDIM; sj++) {
						cg2[sj]=cg[sj]+shifts[si][sj];
						//fold the grid coordinates
						cg2[sj]-=nCells[sj]*floor(double(cg2[sj])/nCells[sj]);
					}

					int iCell2=cg2[0]*nCells[1]*nCells[2]+cg2[1]*nCells[2]+cg2[2];
					//we only need to bother if the corresponding cell is occupied in reference
					if(cellsRef[iCell2].size()>0 and iCell1 <=iCell2) {
						scannedCells.insert(iCell2);
					}
				}

				//loop on the neighboring cells
				for(std::set<int>::iterator itCell2=scannedCells.begin(); itCell2!=scannedCells.end(); itCell2++) {
					int iCell2=*itCell2;

					for(auto itAtom1=cellsSys[iCell1].begin(); itAtom1!=cellsSys[iCell1].end(); itAtom1++ ) {
						std::array<double,NDIM> positions;
						for(int kk=0; kk<NDIM; kk++) {
							positions[kk]=system.getPosition(*itAtom1,kk);
						}
						for(auto itAtom2=cellsRef[iCell2].begin(); itAtom2!=cellsRef[iCell2].end(); itAtom2++ ) {
							//double dr=nearestImageDistance(*itAtom1, *itAtom2, system, bc);

							double dr=nearestImageDistance(positions, *itAtom2, reference, bc);

							if(dr<cutoffs[std::make_pair(reference.getSpecies(*itAtom1),system.getSpecies(*itAtom2))]) {
								if(*itAtom1!=*itAtom2) {
									neighborList[*itAtom1].push_back(*itAtom2);
									if(iCell1!=iCell2) {
										neighborList[*itAtom2].push_back(*itAtom1);
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

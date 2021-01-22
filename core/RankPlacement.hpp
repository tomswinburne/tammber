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

#ifndef RankPlacement_hpp
#define RankPlacement_hpp

#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <algorithm>
#include <random>

struct RankPlacer {
public:

	int seed;
	int rank;
	int taskGroup;
	int nNodes;
	int nSlotsPerNode;
	int nSlotsPerWorker;
	int nWorkManagers;
	bool padNodes;

	//rank, taskGroup,parent,task,seed
	std::map<int,std::tuple<int,int,int,std::string,int> > layout;


	void allocate(int nRanks,std::string task, int parent){
		std::default_random_engine generator;
		generator.seed(seed+taskGroup);
		std::uniform_int_distribution<int> distribution(1,100000);
		int s=distribution(generator);
		for(int i=0; i<nRanks; i++) {
			layout[rank]=std::make_tuple(rank,taskGroup,parent,task,s);
			rank++;
		}
		taskGroup++;
	};

	void pad(){
		if(padNodes) {
			allocate(nSlotsPerNode-rank%nSlotsPerNode,"Nothing", -1);
		}
	};

	RankPlacer(int nNodes_, int nSlotsPerNode_, int nSlotsPerWorker_,  int nWorkManagers_, bool padNodes_, int seed_){
		seed=seed_;
		nNodes=nNodes_;
		nSlotsPerNode=nSlotsPerNode_;
		nSlotsPerWorker=nSlotsPerWorker_;
		nWorkManagers=nWorkManagers_;
		padNodes=padNodes_;


		rank=0;
		taskGroup=0;
		int parent=-1;

		allocate(1,"Splicer",parent);
		pad();
		allocate(1,"PersistentDB",parent);
		pad();
		allocate(1,"InMemoryDB",parent);
		pad();

		if(padNodes) {
			int nNodesLeft=nNodes-3;
			int nNodesPerWM=nNodesLeft/nWorkManagers;
			int nc=nNodesLeft % nWorkManagers;

			for(int i=0; i<nWorkManagers; i++) {
				parent=rank;
				int np=(i<nc ? 1 : 0);

				int nn=nNodesPerWM+np;
				int budget=nn*nSlotsPerNode;
				int r0=rank;

				allocate(1,"WorkManager",-1);
				pad();

				while(budget-(rank-r0) >= nSlotsPerWorker) {
					allocate(nSlotsPerWorker,"Worker",parent);
				}

				int nl=budget-(rank-r0);
				if(nl>0) {
					allocate(nl,"Nothing",-1);
				}
			}
		}
		else{
			int nRanksLeft=nNodes*nSlotsPerNode-3;
			int nRanksPerWM=nRanksLeft/nWorkManagers;
			int nc=nRanksLeft % nWorkManagers;

			for(int i=0; i<nWorkManagers; i++) {
				parent=rank;
				int np=(i<nc ? 1 : 0);

				int nr=nRanksPerWM+np;
				int budget=nr;
				int r0=rank;

				allocate(1,"WorkManager",-1);
				pad();
				while(budget-(rank-r0) >= nSlotsPerWorker) {
					allocate(nSlotsPerWorker,"Worker",parent);
				}
				int nl=budget-(rank-r0);
				if(nl>0) {
					allocate(nl,"Nothing",-1);
				}
			}

		}
		if(nNodes*nSlotsPerNode-rank>0) {
			allocate(nNodes*nSlotsPerNode-rank,"Nothing",-1);
		}
	}
};




#endif

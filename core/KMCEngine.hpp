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



#ifndef KMCEngine_hpp
#define KMCEngine_hpp

#include <thread>
#include <stdio.h>
#include <math.h> /* exp */
#include <boost/algorithm/string.hpp>

#include "AbstractEngine.hpp"
#include "Constants.hpp"
#include "Types.hpp"

#define wct_per_forcecall_us 1   /* Sleep this many us for each approximated "force-call" */
#define forcecalls_per_min   10  /* Assume this many force calls for each quench */


#if 0
template <class System, class Labeler > class KMCEngine : public AbstractEngine<System,Labeler> {
public:

KMCEngine(boost::property_tree::ptree &config,MPI_Comm localComm_, int seed_) : AbstractEngine<System,Labeler>(config, localComm_,seed_) {

	boost::random::random_device rd;
	rng.seed(rd());

	/*
	   //dummy two-state model (Danny's Original Test)
	    boost::random::random_device rd;
	    rng.seed(rd());

	    for(int i=0; i<10; i++) {
	            markovModel.update(1,1);
	            markovModel.update(2,2);
	    }
	    markovModel.update(1,2);
	    markovModel.update(2,1);
	 */
};

private:
virtual void md_impl(Task<System> &task){

	/* Variables */
	double k_intra, k_inter, k_total;

	/* Constants for kMC model */
	static int n_sprbasn_states     = 12;
	static int n_sprbasn            = 100;
	static double Ea_inter_sprbasn  = 0.6;
	static double Ea_intra_sprbasn  = 0.3;
	static int n_path_inter_sprbasn = 2;
	static int n_path_intra_sprbasn = 2;
	static double nu_inter_sprbasn  = 2.59852e12;
	static double nu_intra_sprbasn  = 2.59852e12;
	static double kb                = 8.617e-5;
	int n_global_states = n_sprbasn_states * n_sprbasn;

	/* Determine legnth of md block */
	double dt   = boost::lexical_cast<double>(task.parameters["Timestep"]);
	dt = dt * 1e-12;
	double time = boost::lexical_cast<double>(task.parameters["BlockTime"]);
	time = time * 1e-12;
	int nsteps   = static_cast<int> (time / dt);
	double T = boost::lexical_cast<double>(task.parameters["Temperature"]);

	/* Determine transition rates */
	k_intra = n_path_intra_sprbasn * nu_intra_sprbasn * exp( -Ea_intra_sprbasn / ( kb * T ) ); // Rate to move within superbasin
	k_inter = n_path_inter_sprbasn * nu_inter_sprbasn * exp( -Ea_inter_sprbasn / ( kb * T ) ); // Rate to more out of superbasin
	k_total = k_intra + k_inter;                                                            // Total Transition Rate

	/* Generate random number */
	double r = uniform(rng);

	/* Check for transition */
	int linitial   = (int) task.systems[0].label;
	int linitial_l = linitial %  n_sprbasn_states;
	int lfinal     = linitial;
	if (r <= (1.0 - exp(-k_total * time)))
	{
		/* Yes, transition */
		/* Determine next state using "ring of superbasins" */
		r = uniform(rng);
		double prob_stay = k_intra / k_total;
		//double prob_out1 = prob_stay + (1.0 - prob_stay) / 2.0;
		double prob_out1 = 1- prob_stay;
		int lref = (linitial / n_sprbasn_states)*n_sprbasn_states;

		if (r < prob_stay)
		{
			r = uniform(rng);
			/* Stay in superbasin */
			if (r < 0.5)
			{
				lfinal = linitial_l + 1;
				if (lfinal >= n_sprbasn_states) lfinal = 0;
			}
			else
			{
				lfinal = linitial_l - 1;
				if (lfinal < 0) lfinal = n_sprbasn_states-1;
			}
			lfinal += lref;
		}
		else
		{
			r = uniform(rng);
			/* Leave superbasin */
			if (r < 0.5)
			{
				lfinal = linitial + n_sprbasn_states;
				if (lfinal >= n_global_states) lfinal -= n_global_states;
			}
			else
			{
				lfinal = linitial - n_sprbasn_states;
				if (lfinal < 0) lfinal += n_global_states;
			}
		}
		std::cout<<"linitial: "<<linitial<<" lfinal: "<<lfinal<<std::endl;
	}
	task.systems[0].label = (uint64_t) lfinal;

	/* Sleep to approximate corce-calls */
	int n_sleep = wct_per_forcecall_us * nsteps;
	std::this_thread::sleep_for(std::chrono::microseconds(n_sleep));

	/*
	   //dummy two-state model (Danny's Original Test)
	    uint64_t initial=task.systems[0].label;
	    uint64_t final=markovModel.statistics[initial].sampleSegmentBKL(rng,uniform);
	    task.systems[0].label=final;
	    std::cout<<"md: "<<initial<<" "<<final<<std::endl;
	    //std::this_thread::sleep_for(std::chrono::seconds(1));
	 */

};

virtual void min_impl(Task<System> &task){

	/* Sleep to approximate corce-calls */
	int n_sleep = wct_per_forcecall_us * forcecalls_per_min;
	std::this_thread::sleep_for(std::chrono::microseconds(n_sleep));

};

virtual void forces_impl(Task<System> &task){

	/* Sleep to approximate corce-call */
	int n_sleep = wct_per_forcecall_us;
	std::this_thread::sleep_for(std::chrono::microseconds(n_sleep));

};

virtual void file_init_impl(Task<System> &task){
/*
        System s;
        s.label=1;
        task.systems.clear();
        task.systems.push_back(s);
 */

	std::string filename=boost::lexical_cast<std::string>(task.parameters["Filename"]);
	boost::trim(filename);

	task.systems.clear();

	std::ifstream in;
	in.open(filename.c_str(), std::ios::in);
	Label lb;

	std::cout<<"FILENAME :"<<filename<<" "<<in.good()<<std::endl;
	while(in>>lb) {
		System s;
		s.label=lb;
		s.type=TYPE::OTHER;
		task.systems.push_back(s);
		std::cout<<"SYSTEM :"<<lb<<std::endl;
	}
	in.close();
};

virtual void init_velocities_impl(Task<System> &task){
	/* DO NOTHING */
};

virtual void file_write_impl(Task<System> &task) {
	std::string filename=boost::lexical_cast<std::string>(task.parameters["Filename"]);
	boost::trim(filename);
	std::ofstream out;
	out.open(filename.c_str(), std::ios::out);
	for(int i=0; i<task.systems.size(); i++) {
		out<<task.systems[i].label<<"\n\n";
	}
	out.close();
}

TransitionStatistics markovModel;
boost::random::mt11213b rng;
boost::random::uniform_01<> uniform;

};

#endif

#endif /* KMCEngine_hpp */

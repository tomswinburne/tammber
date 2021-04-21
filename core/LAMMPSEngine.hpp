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
#ifndef LAMMPSEngine_h
#define LAMMPSEngine_h

#include <mpi.h>
#include <stdio.h>
#include <inttypes.h>
#include <random>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>


//#include "AbstractEngine.hpp"
//#include "MDEngine.hpp" not included here (if used, included in CustomTypes)
#include "Log.hpp"
#include "LAMMPSSystem.hpp"
#include "lammps/lammps.h"
#include "lammps/library.h"

using namespace LAMMPS_NS;

/*
	Templated out base class, to allow more than one abstract engine type to be adapted to lammps
	Idea is that we can derive modified abstract classes from MDEngine (c.f. TADEngine)
	*before* the final lammps functionality is included, without any duplication
	In this spirit we do not have a LAMMPSTaskMapper as LAMMPS doesn't determine the tasks
*/

//class LAMMPSEngine : public MDEngine<LAMMPSSystem,MDTaskMapper>  {
template<class MDBaseEngine>
class LAMMPSEngine : public MDBaseEngine {
public:

friend class LAMMPSSystem;
typedef LAMMPSSystem System;
//typedef MDEngine<LAMMPSSystem,LAMMPSTaskMapper> MDBaseEngine;

LAMMPSEngine(boost::property_tree::ptree &config, MPI_Comm localComm_, int seed_) : MDBaseEngine(config,localComm_,seed_) {

	parser.seed(seed_);
	// TODO: how to pass command-line args to LAMMPS at some point?

	log_lammps = config.get<bool>("Configuration.LAMMPSEngine.LogLammps");

	std::string logfile="none";
	if(log_lammps) logfile="log_"+std::to_string(MDBaseEngine::local_rank)+"_"+std::to_string(seed_)+".lammps";

	int argc=5;
	char *lammps_argv[]={(char *)"tammber",(char *)"-screen",(char *)"none",(char *)"-log",(char *)logfile.c_str()};

	LOGGER("Trying to open lammps worker; local rank = "<<MDBaseEngine::local_rank)

	lmp = NULL;
	lmp = new LAMMPS(argc,lammps_argv,localComm_);

	LOGGER("Opened lammps worker ; local rank = "<<MDBaseEngine::local_rank)




	// error check on integer sizes, for storage in System
	// NOTE: avoid error check for int sizes?

	int bigint_bytes = lammps_extract_setting(lmp,(char *) "bigint");
	int tagint_bytes = lammps_extract_setting(lmp,(char *) "tagint");
	int imageint_bytes = lammps_extract_setting(lmp,(char *) "imageint");

	if (bigint_bytes != 8)
		error("LAMMPS must be compiled for 8-byte big integers");
	if (tagint_bytes != 4)
		error("LAMMPS must be compiled for 4-byte atom IDs");
	if (imageint_bytes != 4)
		error("LAMMPS must be compiled for 4-byte image flags");


	//read the LAMMPS scripts

	//bootstrapScript=config.get<std::string>("Configuration.LAMMPSEngine.BootstrapScript");
	mdScript=config.get<std::string>("Configuration.LAMMPSEngine.MDScript");
	minScript=config.get<std::string>("Configuration.LAMMPSEngine.MinScript");
	writeScript=config.get<std::string>("Configuration.LAMMPSEngine.WriteScript");
	initScript=config.get<std::string>("Configuration.LAMMPSEngine.InitScript");
	postInitScript=config.get<std::string>("Configuration.LAMMPSEngine.PostInitScript");
	velocityInitScript=config.get<std::string>("Configuration.LAMMPSEngine.VelocityInitScript");

	int nGPU=config.get<int>("Configuration.LAMMPSEngine.NAccelerators",1);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(nGPU<0) {nGPU=1;}
	int gpuID=rank%nGPU;

	std::unordered_map<std::string, std::string> p;
	p["AcceleratorID"]=boost::str(boost::format("%1%" ) % gpuID );
	p["NAccelerators"]=boost::str(boost::format("%1%" ) % nGPU );

	//bootstrap the driver here
	//this should setup the base LAMMPS parameters (potentials, species, etc.)
	std::string ICstr =	config.get<std::string>("Configuration.InitialConfigurations");

	std::vector<std::string> initialConfs;
	boost::split(initialConfs,ICstr,boost::is_any_of("\n, "));
	for(auto initialConf : initialConfs) {
		boost::trim(initialConf);
		if(initialConf.length()>0) {
			LOGGER("BOOTSTRAP: "<<initialConf)
			p["Filename"] = initialConf;
			break;
		}
	}
	bootstrapScript = parser.parse(initScript,p);
	std::vector<std::string> cmdVector=parser.splitLines(bootstrapScript);
	std::string cmd;
	for(int i=0; i<cmdVector.size(); i++) {
		LOGGER("LAMMPSCommand: "<<cmdVector[i])
		lammps_commands_string(lmp,(char *) cmdVector[i].c_str());
	}

	// Any overwritten functions must be reassigned to dispatchor
	MDBaseEngine::BaseEngine::impls["TASK_DIE"] = LAMMPSEngine::die_impl;
	MDBaseEngine::BaseEngine::impls["TASK_MD"] = LAMMPSEngine::md_impl;
	MDBaseEngine::BaseEngine::impls["TASK_INIT_VELOCITIES"] = LAMMPSEngine::init_velocities_impl;
	MDBaseEngine::BaseEngine::impls["TASK_MIN"] = LAMMPSEngine::min_impl;
	MDBaseEngine::BaseEngine::impls["TASK_INIT_FROM_FILE"] = LAMMPSEngine::file_init_impl;
	MDBaseEngine::BaseEngine::impls["TASK_WRITE_TO_FILE"] = LAMMPSEngine::file_write_impl;
	MDBaseEngine::BaseEngine::impls["TASK_FORCES"] = LAMMPSEngine::forces_impl;
	MDBaseEngine::BaseEngine::impls["TASK_CENTRO"] = LAMMPSEngine::centro_impl;
};

virtual bool failed(){
	if(bool(lammps_has_error(lmp))) {
		char error_message[2048];
		int error_type = lammps_get_last_error_message(lmp,error_message,2048);
		if(MDBaseEngine::local_rank==0) LOGGERA("LAMMPS ERROR! type:"<<error_type<<" msg:"<<error_message)
		return true;
	} else return false;
};


~LAMMPSEngine() {
	lammps_close(lmp);
};

private:
std::function<void(GenericTask&)> die_impl = [this](GenericTask &task) {
	lammps_close(lmp);
};

std::function<void(GenericTask&)> md_impl = [this](GenericTask &task) {

std::unordered_map<std::string,std::string> parameters=extractParameters(task.type,task.flavor,MDBaseEngine::defaultFlavor,MDBaseEngine::taskParameters);

	LAMMPSSystem s;
	extract("State",task.inputData,s);
	LAMMPSSystem *sys = &s;

	int natomsBegin;
	int natomsEnd;

	do {
		transferSystemToLammps(*sys, parameters);
		natomsBegin = (int) *( (int64_t *) lammps_extract_global(lmp,(char *) "natoms"));

		double temperature;
		if(extract("Temperature",task.arguments,temperature))
			parameters["Temperature"] = boost::lexical_cast<std::string>(temperature);

		//figure out how many steps to do
		double dt=safe_extractor<double>(parameters,"Timestep",0.001);


		double time;
		if(!extract("BlockTime",task.arguments,time))
			time = safe_extractor<double>(parameters,"BlockTime",1.0);


		int nsteps = static_cast<int> (time / dt);
		parameters["Nsteps"]=boost::str(boost::format("%1%" ) % nsteps );
		LOGGER("Setting up MD at temperature of "<<temperature<<"K for "<<nsteps<<" steps")

		//parse the command string
		std::string rawCmd = mdScript; //task.parameters["MDScript"];
		std::string parsedCmd=parser.parse(rawCmd, parameters);
		std::vector<std::string> cmdVector=parser.splitLines(parsedCmd);

		//execute the command string
		std::string cmd;
		for(int i=0; i<cmdVector.size(); i++) {
			cmd+=cmdVector[i]+"\n";
		}
		lammps_commands_string(lmp,(char *) cmd.c_str());

		natomsEnd = (int) *((int64_t *) lammps_extract_global(lmp,(char *) "natoms"));
		if(natomsEnd!=natomsBegin && MDBaseEngine::local_rank==0) LOGGERA("ERROR: LAMMPS LOST ATOMS. RESTARTING TASK")
	} while(natomsBegin != natomsEnd );

	transferAtomsFromLammps(s);
	//task.clearOutputs();
	//insert("State",task.outputData,s);

	auto it=task.inputData.find("State");
	Label lb=0;
	extract("Label",task.arguments,lb);
	task.outputData.clear();
	insert("State",lb,it->second.location, it->second.shared,task.outputData,s);
};

std::function<void(GenericTask&)> forces_impl = [this](GenericTask &task) {
	std::unordered_map<std::string,std::string> parameters =
		extractParameters(task.type,task.flavor,MDBaseEngine::defaultFlavor,MDBaseEngine::taskParameters);

	bool reset; if(!extract("Reset",task.arguments,reset)) reset = true;

	bool prepost; if(!extract("PrePost",task.arguments,prepost)) prepost = true;

	std::list<System> slist;
	task.clearOutputs();
	extract("State",task.inputData,slist);

	for(auto &ns: slist) {
		if(reset) {
			LAMMPSSystem *sys = &ns;
			transferSystemToLammps(*sys, parameters);
			reset=false;
		}
		singleForceEnergyCall(ns,false,prepost);
		insert("State",task.outputData,ns);
	}
	task.clearInputs();
};

std::function<void(GenericTask&)> centro_impl = [this](GenericTask &task) {
	std::unordered_map<std::string,std::string> parameters =
		extractParameters(task.type,task.flavor,MDBaseEngine::defaultFlavor,MDBaseEngine::taskParameters);
		int nn; if(!extract("CentroNeighbors",task.arguments,nn)) nn = 8;

		LOGGER("CSNN (C): "<<nn)

		task.clearOutputs();
		System s;
		if (extract("State",task.inputData,s)) {
			LAMMPSSystem *sys = &s;
			transferSystemToLammps(*sys, parameters);
			std::vector<double> csl = calculateCentroSymmetry(s,nn);
			insert("CentroSymmetry",task.outputData,csl);
		}
		task.clearInputs();
};
std::function<void(GenericTask&)> init_velocities_impl = [this](GenericTask &task){
	std::unordered_map<std::string,std::string> parameters=extractParameters(task.type,task.flavor,MDBaseEngine::defaultFlavor,MDBaseEngine::taskParameters);

	LAMMPSSystem s;
	extract("State",task.inputData,s);
	LAMMPSSystem *sys = &s;

	transferSystemToLammps(*sys, parameters);

	double InitTemperature;
	if(extract("InitTemperature",task.arguments,InitTemperature))
		parameters["InitTemperature"]=boost::lexical_cast<std::string>(InitTemperature);

	//parse the command string
	std::string rawCmd = velocityInitScript;//task.parameters["MinScript"];
	std::string parsedCmd=parser.parse(rawCmd, parameters);
	std::vector<std::string> cmdVector=parser.splitLines(parsedCmd);

	//execute the command string

	std::string cmd;
	for(int i=0; i<cmdVector.size(); i++) cmd+=cmdVector[i]+"\n";
	lammps_commands_string(lmp,(char *) cmd.c_str());

	//transfer back
	Label lb=0;
	extract("Label",task.arguments,lb);
	auto it=task.inputData.find("State");
	insert("State",lb,it->second.location, it->second.shared,task.outputData,s);

};

std::function<void(GenericTask&)> min_impl = [this](GenericTask &task){

	std::unordered_map<std::string,std::string> parameters=\
		extractParameters(task.type,task.flavor,MDBaseEngine::defaultFlavor,MDBaseEngine::taskParameters);

	LAMMPSSystem s;
	extract("State",task.inputData,s);
	LAMMPSSystem *sys = &s;


	transferSystemToLammps(*sys, parameters);

	//parse the command string
	std::string rawCmd = minScript;//task.parameters["MinScript"];
	std::string parsedCmd=parser.parse(rawCmd, parameters);
	std::vector<std::string> cmdVector=parser.splitLines(parsedCmd);

	//execute the command string
	std::string cmd;
	for(int i=0; i<cmdVector.size(); i++) {
		cmd+=cmdVector[i]+"\n";
	}
	lammps_commands_string(lmp,(char *) cmd.c_str());

	// get energy
	//singleForceEnergyCall(s,false,true); // options are "no forces" and (re)compute neighbor list
	transferAtomsFromLammps(s);

	//task.clearOutputs();
	//insert("State",task.outputData,s);

	//transfer back
	Label lb=0;
	extract("Label",task.arguments,lb);

	auto it=task.inputData.find("State");
	insert("State",lb,it->second.location, it->second.shared,task.outputData,s);


};

std::function<void(GenericTask&)> file_init_impl = [this](GenericTask &task){

	std::string filename;
	std::unordered_map<std::string,std::string> parameters=extractParameters(task.type,task.flavor,MDBaseEngine::defaultFlavor,MDBaseEngine::taskParameters);


	if( extract("Filename",task.arguments,filename) ) {
		parameters["Filename"]=filename;
		LOGGER("Filename: "<<filename)

		//parse the command string
		std::string rawCmd = initScript;// task.parameters["InitScript"];
		std::string parsedCmd=parser.parse(rawCmd, parameters);
		std::vector<std::string> cmdVector=parser.splitLines(parsedCmd);


		//execute the command string
		std::string cmd;
		for(int i=0; i<cmdVector.size(); i++) {

			LOGGER(cmdVector[i])
			cmd=cmdVector[i];
			lammps_command(lmp,(char *)cmd.c_str());
		}
		//cmd+=cmdVector[i]+"\n";
		//lammps_commands_string(lmp,(char *) cmd.c_str());



	} else error("LAMMPS file init: unrecognized parameters");


	//initialize the system from the LAMMPS instance
	LAMMPSSystem s;
	transferSystemFromLammps(s);

	Label lb=0;
	extract("Label",task.arguments,lb);

	auto lit=task.returns.erase("Label");
	insert("Label", task.returns, lb);

	auto it=task.inputData.find("State");
	insert("State",lb,it->second.location, it->second.shared,task.outputData,s);
};

std::function<void(GenericTask&)> file_write_impl = [this](GenericTask &task) {

	std::unordered_map<std::string,std::string> parameters=extractParameters(task.type,task.flavor,MDBaseEngine::defaultFlavor,MDBaseEngine::taskParameters);

	LAMMPSSystem s;
	bool gotState=extract("State",task.inputData,s);
	LAMMPSSystem *sys = &s;
	transferSystemToLammps(*sys, parameters);

	std::string filename;
	if(extract("Filename",task.arguments,filename) ) {
		filename=parser.parse(filename, parameters);
		boost::filesystem::path p {filename};
		boost::filesystem::path dir=p.parent_path();
		try {
			boost::filesystem::create_directories(dir);
		}
		catch(...) {};

		parameters["Filename"]=filename;

		//parse the command string
		std::string rawCmd = writeScript;//task.parameters["WriteScript"];
		std::string parsedCmd=parser.parse(rawCmd, parameters);

		std::vector<std::string> cmdVector=parser.splitLines(parsedCmd);


		//execute the command string
		std::string cmd;
		for(int i=0; i<cmdVector.size(); i++) cmd+=cmdVector[i]+"\n";
		lammps_commands_string(lmp,(char *) cmd.c_str());

	} else error("LAMMPS file init: unrecognized parameters");

};

void transferAtomsFromLammps(System &s){
	lammps_command(lmp,(char *) "run 0");
	LAMMPSSystem *sys = &s;
	lammps_gather(lmp,(char *) "id",LAMMPS_INT,1,&sys->id[0]);
	lammps_gather(lmp,(char *) "type",LAMMPS_INT,1,&sys->species[0]);
	lammps_gather(lmp,(char *) "x",LAMMPS_DOUBLE,3,&sys->x[0]);
	lammps_gather(lmp,(char *) "v",LAMMPS_DOUBLE,3,&sys->v[0]);
	if (sys->qflag) lammps_gather(lmp,(char *) "q",LAMMPS_DOUBLE,1,&sys->q[0]);
	if (sys->fflag) lammps_gather(lmp,(char *) "f",LAMMPS_DOUBLE,3,&sys->f[0]);

	void * lmpE = lammps_extract_compute(lmp,(char *) "pe",0,0);
	if (lmpE == NULL) {
		lammps_command(lmp,(char *) "compute pe all pe");
		lammps_command(lmp,(char *) "run 0");
		lmpE = lammps_extract_compute(lmp,(char *) "pe",0,0);
	}
	s.setEnergy(*((double *) lmpE));
};

void transferSystemFromLammps(System &s){
	LAMMPSSystem *sys = &s;
	int triclinic = *((int *) lammps_extract_global(lmp,(char *) "triclinic"));
	int qflag = *((int *) lammps_extract_global(lmp,(char *) "q_flag"));
	LOGGER("LAMMPSEngine::transferSystemFromLammps : q_flag="<<qflag)
	sys->setNTimestep(0);
	sys->qflag = qflag;
	int natoms = (int) *((int64_t *)
	                     lammps_extract_global(lmp,(char *) "natoms"));
	sys->setNAtoms(natoms);

	// box info

	double *boxlo = (double *) lammps_extract_global(lmp,(char *) "boxlo");
	double *boxhi = (double *) lammps_extract_global(lmp,(char *) "boxhi");
	sys->xy = *((double *) lammps_extract_global(lmp,(char *) "xy"));
	sys->xz = *((double *) lammps_extract_global(lmp,(char *) "xz"));
	sys->yz = *((double *) lammps_extract_global(lmp,(char *) "yz"));

	int *periodicity = (int *) lammps_extract_global(lmp,(char *) "periodicity");
	for(int i=0; i<3; i++) {
		sys->boxlo[i]=boxlo[i];
		sys->boxhi[i]=boxhi[i];
		sys->periodicity[i]=periodicity[i];
	}
	sys->updateCell();
	transferAtomsFromLammps(s);
};

void transferSystemToLammps(System &sys, std::unordered_map<std::string, std::string> &parameters){
	int natoms = sys.getNAtoms();
	char cmdc[512];
	lammps_command(lmp,(char *) "delete_atoms group all");

	bool triclinic = *((int *) lammps_extract_global(lmp,(char *) "triclinic"));
	if (triclinic)
		sprintf(cmdc,"change_box all "
		        "x final %g %g y final %g %g z final %g %g "
		        "xy final %g xz final %g yz final %g units box",
		        sys.boxlo[0],sys.boxhi[0],
		        sys.boxlo[1],sys.boxhi[1],
		        sys.boxlo[2],sys.boxhi[2],
		        sys.xy,sys.xz,sys.yz);
	else
		sprintf(cmdc,"change_box all "
		        "x final %g %g y final %g %g z final %g %g units box",
		        sys.boxlo[0],sys.boxhi[0],
		        sys.boxlo[1],sys.boxhi[1],
		        sys.boxlo[2],sys.boxhi[2]);
	lammps_command(lmp,cmdc);

	// NOTE: should I also (re)set periodicity?   or check has not changed?
	// NOTE: periodicity is a double within System?
	// NOTE: what about remapping of atom back into box with image change?


	lammps_create_atoms(lmp,natoms,&sys.id[0],&sys.species[0],&sys.x[0],&sys.v[0],NULL,1);
	//if(sys.qflag) lammps_scatter_subset(lmp,(char *) "q",LAMMPS_DOUBLE,1,sys.id.size(),&sys.id[0],&sys.q[0]);
	LOGGER("LAMMPSEngine::transferSystemToLammps : q_flag="<<sys.qflag)
	if(sys.qflag) lammps_scatter(lmp,(char *) "q",LAMMPS_DOUBLE,1,&sys.q[0]);

	//parse the command string
	std::string rawCmd = postInitScript; //parameters["PostInitScript"];
	std::string parsedCmd=parser.parse(rawCmd, parameters);
	std::vector<std::string> cmdVector=parser.splitLines(parsedCmd);

	//execute the command string
	std::string cmd;
	for(int i=0; i<cmdVector.size(); i++) cmd+=cmdVector[i]+"\n";
	lammps_commands_string(lmp,(char *) cmd.c_str());
};

void singleForceEnergyCall(System &s, bool noforce=false,bool prepost=true) {

	//lammps_scatter_subset(lmp,(char *) "x",LAMMPS_DOUBLE,3,s.id.size(),&s.id[0],&s.x[0]);
	lammps_scatter(lmp,(char *) "x",LAMMPS_DOUBLE,3,&s.x[0]);
	if(prepost) lammps_command(lmp,(char *) "run 0"); // annoying but need to build nlist
	else lammps_command(lmp,(char *) "run 1 pre no post no"); // OK with check..
	int nAtoms = s.getNAtoms();
	if(!noforce) {
		if(s.fflag==0) s.fflag=1;
		if(s.f.size() != NDIM*nAtoms) s.f.resize(NDIM*nAtoms,0.0);
		lammps_gather_subset(lmp,(char *)"f",LAMMPS_DOUBLE,3,s.id.size(),&s.id[0],&s.f[0]);
	}

	// if a compute called 'pe' does not exists we create it- should only happen once
	void * lmpE = lammps_extract_compute(lmp,(char *) "pe",0,0);
	if (lmpE == NULL) {
		lammps_command(lmp,(char *) "compute pe all pe");
		lammps_command(lmp,(char *) "run 0");
		lmpE = lammps_extract_compute(lmp,(char *) "pe",0,0);
	}
	s.setEnergy(*((double *) lmpE));
};


std::vector<double> calculateCentroSymmetry(System &s,int nn=8) {
	//lammps_scatter_subset(lmp,(char *) "x",LAMMPS_DOUBLE,3,s.id.size(),&s.id[0],&s.x[0]);
	lammps_scatter(lmp,(char *) "x",LAMMPS_DOUBLE,3,&s.x[0]);

	lammps_command(lmp,(char *) "run 0"); // annoying but need to build nlist
	std::string cc="compute calcCS all centro/atom "+std::to_string(nn);
	lammps_command(lmp,(char *) cc.c_str());
	lammps_command(lmp,(char *) "run 0");
	std::vector<double> res(s.id.size(),0.0);
	lammps_gather(lmp,(char *)"c_calcCS",LAMMPS_DOUBLE,1,&res[0]);
	lammps_command(lmp, (char *)"uncompute calcCS");
	return res;
};



// LAMMPS specific variables
int me,nprocs;               // MPI info
LAMMPS *lmp;                  // instance of LAMMPS
void error(const char *str){
	if (me == 0) printf("ERROR: %s\n",str);
	MPI_Abort(MPI_COMM_WORLD,1);
};

private:
ParameterParser parser;
std::size_t identityHash;
std::size_t coordinatesHash;

bool log_lammps;
std::string bootstrapScript;
std::string mdScript;
std::string minScript;
std::string writeScript;
std::string initScript;
std::string postInitScript;
std::string velocityInitScript;
};



#endif /* LAMMPSEngine_h */

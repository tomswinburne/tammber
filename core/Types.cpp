#include "Types.hpp"
#include "Log.hpp"

////////////////////////////////////////////////////////////////////////////////////////////

StateStatistics::StateStatistics(){
	nSegments=0;
	nTransitions=0;
};

StateStatistics::StateStatistics(uint64_t label_){
	nSegments=0;
	nTransitions=0;
	label=label_;
};

void StateStatistics::update(uint64_t finalState){
	nSegments++;
	if(finalState!=label) {
		counts[finalState]+=1;
		nTransitions+=1;
	}

	//std::cout<<"UPDATE: "<<label<<" "<<finalState<<" "<<nSegments<<" "<<nTransitions<<std::endl;
};

void StateStatistics::update(StateStatistics s){
	nSegments+=s.nSegments;
	nTransitions+=s.nTransitions;
	for(auto it=s.counts.begin(); it!=s.counts.end(); it++) {
		counts[it->first]+=s.counts[it->first];
	}
};



void StateStatistics::sampleSegmentBKL(double unknownPrior, boost::random::mt11213b &rand, boost::random::uniform_01<> &uniform,Label &final,  bool &absorbed){
	unknownPrior=std::max(unknownPrior,1e-10);
	double r1=uniform(rand);
	final=label;

	double rN=r1*(nSegments+unknownPrior+1);

	double rAccum=unknownPrior;

	//this was absorbed into the unknown
	if(rAccum>=rN) {
		absorbed=true;
	}

	for(auto it=counts.begin(); it!=counts.end(); it++) {
		rAccum+=it->second;
		//this is a real transition
		if(rAccum>=rN) {
			final=it->first;
			break;
		}
	}
	absorbed=false;
};

void StateStatistics::sampleEscapeBKL(double unknownPrior, boost::random::mt11213b &rand, boost::random::uniform_01<> &uniform, Label &final,  bool &absorbed,int &nConsumed){

	unknownPrior=std::max(unknownPrior,1e-10);

	double pstay=unknownPrior;
	for(auto it=counts.begin(); it!=counts.end(); it++) {
		pstay+=it->second;
		//std::cout<<it->second<<" "<<pstay<<std::endl;
	}
	pstay/=(nSegments+unknownPrior+1);
	pstay=1-pstay;


	//std::cout<<"pstay: "<<pstay<<" "<<nSegments<<" "<<unknownPrior <<std::endl;
	//sample a number of self transition before an escape
	boost::random::negative_binomial_distribution<int,double> nbino(1,1-pstay);
	nConsumed=nbino(rand);

	//std::cout<<"pstay: "<<pstay<<" "<<nConsumed <<std::endl;

	//sample a transition
	double r1=uniform(rand);
	double rN=r1*(nTransitions+unknownPrior);
	double rAccum=unknownPrior;

	final=label;

	//this was absorbed into the unknown
	if(rAccum>=rN) {
		absorbed=true;
	}

	for(auto it=counts.begin(); it!=counts.end(); it++) {
		rAccum+=it->second;
		//this is a real transition
		if(rAccum>=rN) {
			final=it->first;
			break;
		}
	}
	absorbed=false;

};



void StateStatistics::clear(){
	nSegments=0;
	counts.clear();
};

////////////////////////////////////////////////////////////////////////////////////////////

TransitionStatistics::TransitionStatistics(){
	boost::random::random_device rd;
	rng.seed(rd());
};

void TransitionStatistics::clear(){
	statistics.clear();
};

void TransitionStatistics::update(uint64_t initialState, uint64_t finalState){

	if(statistics.count(initialState)==0) {
		StateStatistics s(initialState);
		statistics.insert(std::make_pair(initialState,s));
	}
	statistics.find(initialState)->second.update(finalState);
};

void TransitionStatistics::assimilate(TransitionStatistics &s){
	for(auto it=s.statistics.begin(); it!=s.statistics.end(); it++) {
		if(statistics.count(it->first)==0) {
			StateStatistics s(it->first);
			statistics.insert(std::make_pair(it->first,s));
		}
		statistics.find(it->first)->second.update(it->second);
	}
	s.clear();
};

void TransitionStatistics::sampleSegmentBKL(Label &lb,double unknownPrior, Label &final, bool &absorbed){
	//unknown state
	final=lb;
	absorbed=true;
	if( statistics.count(lb)==0 ) {
		StateStatistics s(lb);
		statistics.insert(std::make_pair(lb,s));
	}
	statistics.find(lb)->second.sampleSegmentBKL(unknownPrior,rng,uniform,final,absorbed);

};

void TransitionStatistics::sampleEscapeBKL(Label &lb,double unknownPrior, Label &final, bool &absorbed, int &nSegments){
	//unknown state
	final=lb;
	absorbed=true;
	nSegments=0;

	if( statistics.count(lb)==0 ) {
		StateStatistics s(lb);
		statistics.insert(std::make_pair(lb,s));
	}
	statistics.find(lb)->second.sampleEscapeBKL(unknownPrior,rng,uniform,final,absorbed,nSegments);

};



int TransitionStatistics::size(){
	return statistics.size();
};

////////////////////////////////////////////////////////////////////////////////////////////


Trajectory::Trajectory(){
	(*this).length=0;
	(*this).overhead_=0;
	(*this).nSplice=1;
};

void Trajectory::print(){
	std::cout<<"============"<<std::endl;

	for(auto it=visits.begin(); it!=visits.end(); it++) {
		std::cout<<it->label<<" "<<it->duration<<std::endl;
	}
};

void Trajectory::log(){
	#ifdef USE_BOOST_LOG
	boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
	BOOST_LOG_SEV(lg, boost::log::trivial::info) <<"============";
	for(auto it=visits.begin(); it!=visits.end(); it++) {
		//std::cout<<it->label<<" "<<it->duration<<std::endl;
		BOOST_LOG_SEV(lg, boost::log::trivial::info) <<it->label<<" "<<it->duration;
	}
	#else
	std::cout<<"============\n";
	for(auto it=visits.begin(); it!=visits.end(); it++) {
		std::cout<<it->label<<" "<<it->duration<<std::endl;
	}
	#endif
};


void Trajectory::clear(){
	visits.clear();
	length=0;
	index=0;
	nSplice=0;
	overhead_=0;
};

//append a visit to a trajectory. Extends the last visit if possible
void Trajectory::appendVisit(Visit &v, bool extend){
	(*this).length+=v.duration;

	if( extend and !(*this).empty() and (*this).back().label==v.label) {
		//extend the last visit
		(*this).back().duration+=v.duration;
		return;
	}
	//new visit
	(*this).visits.push_back(v);
};


//splice two trajectories
bool Trajectory::splice(Trajectory &t){

	//do not splice empty trajectories
	if( (*this).empty() or t.empty() ) {
		return false;
	}

	//splice only if the beginning and end match
	if(  (*this).back().label == t.front().label ) {
		(*this).back().duration+=t.front().duration;
		t.visits.pop_front();
		(*this).visits.splice((*this).visits.end(), t.visits);
		(*this).length+=t.length;
		(*this).overhead_+=t.overhead_;
		//consume t
		t.clear();

		(*this).nSplice++;

		return true;
	}
	return false;
};


void Trajectory::truncate(int maxDuration, Trajectory &leftover){
	leftover.clear();
	//we don't need to do anything
	if( (*this).duration()<=maxDuration) {
		return;
	}


	int dt=(*this).duration();

	leftover=*this;
	(*this).clear();

	int nUsed=0;
	while(nUsed<maxDuration) {
		Visit front=leftover.front();
		if(nUsed+front.duration<=maxDuration) {
			//use the whole visit
			(*this).appendVisit(front);
			nUsed+=front.duration;
			leftover.pop_front();
		}
		else{
			//take only a fraction
			Visit v;
			v.label=front.label;
			v.duration=maxDuration-nUsed;
			(*this).appendVisit(v);
			leftover.front().duration-=(maxDuration-nUsed);
			leftover.length-=(maxDuration-nUsed);
			nUsed=maxDuration;
		}
	}

	if(dt!=(*this).duration()+leftover.duration() ) {
		std::cout<<"ERROR!!! LOST SOME BLOCKS"<<std::endl;
	}
};


uint64_t Trajectory::duration(){
	return length;
};

uint64_t& Trajectory::overhead(){
	return overhead_;
};

bool Trajectory::empty(){
	return visits.size()==0;
};

Visit& Trajectory::back(){
	return visits.back();
};

Visit& Trajectory::front(){
	return visits.front();
};

void Trajectory::pop_back(){
	if(visits.size()>0) {
		length-=visits.back().duration;
		visits.pop_back();
	}
};

void Trajectory::pop_front(){
	if(visits.size()>0) {
		length-=visits.front().duration;
		visits.pop_front();
	}
};

////////////////////////////////////////////////////////////////////////////////////////////

void SegmentDatabase::clear(){
	db.clear();
};

int SegmentDatabase::size(){
	return db.size();
};

int SegmentDatabase::count(Label lb, int flavor){
	auto key=std::make_pair(lb,flavor);
	if(db.count(key)) {
		return db[key].size();
	}
	return 0;
};

bool SegmentDatabase::front(Label lb, int flavor, Trajectory &t){
	if(count(lb,flavor)==0) {
		return false;
	}
	auto key=std::make_pair(lb,flavor);
	t=db[key].front();
	return true;
};

void SegmentDatabase::pop_front(Label lb, int flavor){
	if(count(lb,flavor)!=0) {
		auto key=std::make_pair(lb,flavor);
		db[key].pop_front();
		if(db[key].empty()) {
			db.erase(key);
		}
	}
};

void SegmentDatabase::add(int flavor, Trajectory &t){
	if(not t.empty() and not (t.duration()==0)) {
		Label lb=t.front().label;
		auto key=std::make_pair(lb,flavor);
		//if trajectory exist in that bin, try to splice at back. Otherwise, add
		if(db.count(key)>0) {
			if(db[key].empty() or (not db[key].back().splice(t))) {
				db[key].push_back(t);
			}
		}
		else{
			db[key]=std::deque<Trajectory>();
			db[key].push_back(t);
		}
	}
	else{
		std::cout<<"SegmentDatabase::add WARNING! Adding an empty trajectory"<<std::endl;
	}
};


void SegmentDatabase::log(){
	#ifdef USE_BOOST_LOG
	boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
	BOOST_LOG_SEV(lg, boost::log::trivial::info) <<"SEGMENT DB ";
	for(auto it=db.begin(); it!=db.end(); it++) {
		//std::cout<<it->first.first<<" "<<it->first.second<<std::endl;
		BOOST_LOG_SEV(lg, boost::log::trivial::info) <<it->first.first<<" "<<it->first.second<<" "<<it->second.size();
		for(auto itt=it->second.begin(); itt!=it->second.end(); itt++) {
			itt->log();
		}
		BOOST_LOG_SEV(lg, boost::log::trivial::info) <<" - ";
	}
	#else
	std::cout<<"SEGMENT DB\n";
	for(auto it=db.begin(); it!=db.end(); it++) {
		std::cout<<it->first.first<<" "<<it->first.second<<" "<<it->second.size()<<std::endl;
		for(auto itt=it->second.begin(); itt!=it->second.end(); itt++) itt->log();
		std::cout<<"-\n";
	}
	#endif
};


void SegmentDatabase::print(){

	for(auto it=db.begin(); it!=db.end(); it++) {
		std::cout<<it->first.first<<" "<<it->first.second<<std::endl;
		for(auto itt=it->second.begin(); itt!=it->second.end(); itt++) {
			itt->print();
		}
	}
};

int SegmentDatabase::duration(){
	int dura=0;
	for(auto it=db.begin(); it!=db.end(); it++) {
		//std::cout<<it->first.first<<" "<<it->first.second<<std::endl;
		for(auto itt=it->second.begin(); itt!=it->second.end(); itt++) {
			dura+=itt->duration();
		}
	}
	return dura;
};

int SegmentDatabase::duration(Label lb, int flavor){
	int dura=0;
	auto key=std::make_pair(lb,flavor);

	for(auto it=db[key].begin(); it!=db[key].end(); it++) {
		dura+=it->duration();
	}
	return dura;
};



void SegmentDatabase::merge(SegmentDatabase &supp){
	//loop over initial states
	for(auto it=supp.db.begin(); it!=supp.db.end(); it++) {
		//loop over trajectories
		int flavor=it->first.second;
		for(auto itt=it->second.begin(); itt!=it->second.end(); itt++) {
			Trajectory &t=*itt;
			this->add(flavor,t);
		}
	}
	supp.clear();
};

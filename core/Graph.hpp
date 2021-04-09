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



#ifndef Graph_hpp
#define Graph_hpp

#include <iostream>
#include <stdio.h>
#include <map>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/function.hpp>
#include <boost/functional/factory.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/optional/optional.hpp>
#include <Eigen/Dense>
#include <boost/timer/timer.hpp>

#include "AbstractSystem.hpp"
#include "NeighborList.hpp"
#include "Hash.hpp"
#include "Log.hpp"



//using namespace boost;

extern "C"
{
#include "nauty/traces.h"
#include "nauty/nausparse.h"
}


typedef std::list<std::pair< std::map<int,int>, std::map<int,int> > >  StateToStateMappings;

typedef std::pair< std::map<int,int>, std::map<int,int> > StateToStateMapping;



struct vertex_p {
	vertex_p(){
		specie=0;
		uniqueIndex=-1;
		consecutiveIndex=-1;
		color=0;
	};
	int specie;
	int uniqueIndex;
	int consecutiveIndex;
	int color;
};

struct edge_p {
	edge_p(){
		type=0;
	};
	int type;
};

struct Distance {
	int i;
	double d;

	bool operator<(Distance other) const
	{
		return d < other.d;
	}

};





struct Interval {
	int i;
	int j;
	double di;
	double dj;
	double dd;
	double score;
	double score1;
	double score2;
	bool operator<(Interval other) const
	{
		return dd < other.dd;
	}
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, vertex_p,edge_p> graph_type;

class AbstractStateLabeler {
public:
virtual void initialize(boost::property_tree::ptree &config)=0;
virtual uint64_t hash(AbstractSystem &s, bool canonical)=0;
virtual int connectedComponents(AbstractSystem &s, std::vector<int> &cluster_occ)=0;
virtual void isomorphicMap(AbstractSystem &s, AbstractSystem &t, std::map<int,int> &map)=0;
virtual void isomorphicSelfMaps(AbstractSystem &s, std::list<std::map<int,int>> &maps)=0;
virtual void canonicalMap(AbstractSystem &s, std::map<int,int> &map, uint64_t &label){
	map.clear();
};

};


// Custom callback for vf2_graph_iso
template <typename Graph1,
					typename Graph2>
struct vf2_map_builder {

	vf2_map_builder(const Graph1& graph1, const Graph2& graph2)
		: graph1_(graph1), graph2_(graph2) {}

	template <typename CorrespondenceMap1To2,
						typename CorrespondenceMap2To1>
	bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1 g) {
		std::map<int,int> map;
		BGL_FORALL_VERTICES_T(v, graph1_, Graph1)
			map.insert(std::make_pair(graph1_[v].uniqueIndex,graph2_[get(f,v)].uniqueIndex));
		//get(boost::vertex_index_t(), graph1_, v),get(boost::vertex_index_t(), graph2_, get(f, v))));
		maps->push_back(map);
		return true;
	}
public:
	std::list< std::map<int,int> > *maps;
private:
	const Graph1& graph1_;
	const Graph2& graph2_;
};


class ConnectivityGraphStateLabeler : public AbstractStateLabeler {
public:
std::map<std::pair<int,int>,double> cutoffs;
std::map<std::pair<int,int>,double> mcutoffs;
std::set<int> distinguishableSpecies;
std::map<int,int> TypeMap; // between MDEngineTypes
bool canonical;
double skin;
double cut;

virtual void initialize(boost::property_tree::ptree &config) {
	std::string ts=config.get<std::string>("Configuration.StateLabeler.DistinguishableSpecies", "");
	boost::trim(ts);
	if(ts.size()>0) {
		std::vector<std::string> strs;
		boost::split(strs,ts,boost::is_any_of(","));
		for(auto it=strs.begin(); it!=strs.end(); it++) {
			std::string s=*it;
			boost::trim(s);
			int sp=boost::lexical_cast<int>(s);

			distinguishableSpecies.insert(sp);
		}
	}
	std::set<int> bond_set; // to replace typemaps
	BOOST_FOREACH(boost::property_tree::ptree::value_type &v, config.get_child("Configuration.StateLabeler.Bonds")) {
		std::string s=v.second.get<std::string>("Between");
		boost::trim(s);
		std::vector<std::string> sp;
		boost::split(sp, s, boost::is_any_of("\t "), boost::token_compress_on);
		int i0=boost::lexical_cast<int>(sp[0]);
		int i1=boost::lexical_cast<int>(sp[1]);
		bond_set.insert(i0);
		bond_set.insert(i1);
		double cut=v.second.get<double>("Cutoff");
		mcutoffs[std::make_pair(i0,i1)]=cut;
		mcutoffs[std::make_pair(i1,i0)]=cut;
	}

	// Mapping between species and type for MD ENGINE- many species for a single type
	TypeMap.clear(); // can't be a vector as may not be sequential
	boost::optional< boost::property_tree::ptree& >
		has_type_maps = config.get_child_optional("Configuration.StateLabeler.TypeMaps");

	if( has_type_maps ) {
		BOOST_FOREACH(boost::property_tree::ptree::value_type &v, config.get_child("Configuration.StateLabeler.TypeMaps")) {
			std::string s=v.second.data();
			boost::trim(s);
			std::vector<std::string> sp;
			boost::split(sp, s, boost::is_any_of("\t "), boost::token_compress_on);
			int i0=boost::lexical_cast<int>(sp[0]);
			int i1=boost::lexical_cast<int>(sp[1]);
			TypeMap.insert(std::make_pair(i0,i1));
			LOGGER("TypeMap: "<<i0<<" -> "<<i1)
		}
	}
	if(TypeMap.size()==0) for(auto i: bond_set)
		TypeMap.insert(std::make_pair(i,i)); // backup option

	// need to inverse TypeMap here- so go through double loop
	cutoffs.clear();
	for(auto &s0: TypeMap) for(auto &s1: TypeMap) {
		cutoffs[std::make_pair(s0.first,s1.first)] = mcutoffs[std::make_pair(s0.second,s1.second)];
		LOGGER(s0.second<<","<<s1.second<<" -> "<<s0.first<<","<<s1.first<<" = "<<mcutoffs[std::make_pair(s0.second,s1.second)])
	}
};

virtual uint64_t hash(AbstractSystem &s, bool canonical){
	graph_type g;
	buildGraph(s, cutoffs, g);

	std::map<int,int> map;
	uint64_t hash=hashGraph(g, canonical, map);
	//hash it
	uint32_t lh=fmix32(s.getEpoch());
	//add it to the hash
	hash = hash ^ lh;
	return hash;
};

virtual int connectedComponents(AbstractSystem &s, std::vector<int> &cluster_occ){
	graph_type g;
	buildGraph(s, cutoffs, g);
	cluster_occ.clear();
	cluster_occ.resize(num_vertices(g),0);
	int clusters = boost::connected_components(g,&cluster_occ[0]);
	return clusters;
};



void buildGraph(AbstractSystem &system,std::map<std::pair<int,int>,double> &cutoffs, graph_type &graph ){
	//boost::timer::auto_cpu_timer t;

	NeighborList nList;
	AbstractSystem *s=dynamic_cast<AbstractSystem*>(&system);
	nList.build(*s,cutoffs);

	//nList.build(system,cutoffs);

	graph.clear();

	boost::graph_traits<graph_type>::vertex_descriptor v;
	boost::graph_traits<graph_type>::vertex_descriptor u;
	boost::vertex_bundle_type<graph_type>::type vp;
	boost::edge_bundle_type<graph_type>::type ep;
	std::map<int, boost::graph_traits<graph_type>::vertex_descriptor > vertexMap;

	int nAtoms=system.getNAtoms();

	for(int i=0; i<nAtoms; i++) {
		vp.consecutiveIndex=i;
		vp.uniqueIndex=system.getUniqueID(i);
		vp.specie=TypeMap[system.getSpecies(i)];

		int color;
		if(distinguishableSpecies.count(vp.specie)>0) {
			color=(vp.uniqueIndex+nAtoms+666);
		}
		else{
			color=vp.specie;
		}

		vp.color=color;
		v=add_vertex(vp,graph);
		vertexMap[i]=v;
	}


	for(int i=0; i<nAtoms; i++) {
		for(int j=0; j<nList.getNumberOfNeighbors(i); j++) {
			int ij=nList.getNeighbor(i,j);
			if(i>ij) {
				ep.type=0;
				v=vertexMap[i];
				u=vertexMap[ij];
				add_edge(v,u,ep,graph);
			}
		}
	}


};

std::map<int,int> canonicalSort(graph_type &ga){
	std::map<int,int> canonicalMap;


	//map boost iterator order to consecutive order
	std::map<int, int > vertexMap;
	std::map<int, int > vertexMapUnique;
	std::map<int, int > ivertexMap;

	std::map<int, std::list<int> > colorMap;
	std::map<int, int > cMap;

	int i=0;
	auto vi= boost::vertices(ga);
	for(auto itv=vi.first; itv!=vi.second; itv++) {
		//vertex color
		int color=ga[*itv].color;
		colorMap[color].push_back(i);
		//index map
		vertexMap[i]=ga[*itv].consecutiveIndex;
		vertexMapUnique[i]=ga[*itv].uniqueIndex;
		ivertexMap[ga[*itv].consecutiveIndex]=i;
		cMap[i]=ga[*itv].color;
		i++;
	}

	unsigned long nVerts = num_vertices(ga);
	unsigned long nEdges=num_edges(ga);



	//convert to TRACES format
	static DEFAULTOPTIONS_TRACES(options);
	TracesStats stats;
	SG_DECL(sg1);
	SG_DECL(cg1);
	DYNALLSTAT(int,lab1,lab1_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);

	options.getcanon = true;
	options.defaultptn = false;

	SG_ALLOC(sg1,nVerts,2*nEdges,"malloc");
	DYNALLOC1(int,lab1,lab1_sz,nVerts,"malloc");
	DYNALLOC1(int,ptn,ptn_sz,nVerts,"malloc");
	DYNALLOC1(int,orbits,orbits_sz,nVerts,"malloc");

	sg1.nv = int(nVerts); /* Number of vertices */
	sg1.nde = int(2*nEdges); /* Number of directed edges */

	//colors
	i=0;
	for(auto itc=colorMap.begin(); itc!=colorMap.end(); itc++) {
		for(auto ita=itc->second.begin(); ita!=itc->second.end(); ita++) {
			lab1[i]=*ita;
			ptn[i]=1;
			i++;
		}
		ptn[i-1]=0;
	}

	i=0;
	int ie=0;
	//create the edges in traces ordering

	for(auto itv=vi.first; itv!=vi.second; itv++) {

		unsigned long na=boost::out_degree(*itv, ga);
		sg1.d[i]=int(na);
		sg1.v[i]=ie;


		auto va=boost::adjacent_vertices(*itv,ga);
		int j=0;
		for(auto itva=va.first; itva!=va.second; itva++) {
			sg1.e[sg1.v[i]+j]=ivertexMap[ga[*itva].consecutiveIndex];
			j++;
		}
		ie+=j;
		i++;
	}

	//call traces
	Traces(&sg1,lab1,ptn,orbits,&options,&stats,&cg1);



	//convert to map unique indices
	for(int i=0; i<nVerts; i++) {

		//this is a distinguishable specie. We should not remap
		if(cMap[i]>666) {
			int tracesIndex=i+1;
			int uniqueIndex=vertexMapUnique[tracesIndex];
			canonicalMap[uniqueIndex]=uniqueIndex;
		}
		else{
			int tracesIndex=i+1;
			int uniqueIndex=vertexMapUnique[tracesIndex];

			int remappedTraceIndex=lab1[i];
			int remappedUniqueIndex=vertexMapUnique[remappedTraceIndex];

			canonicalMap[remappedUniqueIndex]=tracesIndex;
		}
	}

	/*

  //keep a set of indices that are available for remapping.
  std::set<int> flexibleIndices;
  i=0;
  for(auto itv=vi.first; itv!=vi.second; itv++) {
    int color=ga[*itv].color;
    int uniqueIndex=vertexMapUnique[i];
    if(color<=666) flexibleIndices.insert(uniqueIndex);
    i++;
  }

	i=0;
	std::map<int,int> uniqueC;
	//proposal remapping. Has to be corrected for constraints on which atom is allowed to be remapped
	for(auto itv=vi.first; itv!=vi.second; itv++) {
		int color=ga[*itv].color;
		//int tracesIndex=i+1;
		int remappedTraceIndex=lab1[i];
		int uniqueIndex=vertexMapUnique[i];
		int remappedUniqueIndex=vertexMapUnique[remappedTraceIndex];
		uniqueC[remappedUniqueIndex]=color;
		canonicalMap[remappedUniqueIndex]=uniqueIndex;
		i++;
	}

	std::set<int> io;
	std::set<int> ir;
	std::map<int,int> constrainedCanonicalMap;
	int cRemap;
	for(auto it=canonicalMap.begin(); it!=canonicalMap.end(); it++) {
		int color=uniqueC[it->first];
		//This is a distinguishable particle. it should not change index.
		if(color>666) cRemap=it->second;
		//this one can change index. Pick the first one available.
		else {
			cRemap=*flexibleIndices.begin();
			flexibleIndices.erase(flexibleIndices.begin());
		}
		constrainedCanonicalMap[cRemap]=it->second;
		io.insert(it->second);
		ir.insert(cRemap);
	}
	canonicalMap=constrainedCanonicalMap;
	assert(io==ir);
	*/
	//free traces memory
	SG_FREE(sg1);
	SG_FREE(cg1);
	DYNFREE(lab1,lab1_sz);
	DYNFREE(ptn,ptn_sz);
	DYNFREE(orbits,orbits_sz);

	return canonicalMap;
};

void canonicalMap(AbstractSystem &s, std::map<int,int> &map, uint64_t &label){
	graph_type g;
	buildGraph(s, cutoffs, g);
	map.clear();
	label=hashGraph(g,true,map);
};

virtual void isomorphicMap(AbstractSystem &s, AbstractSystem &t, std::map<int,int> &map) {
	std::map<int,int> map1,map2,imap1;
	graph_type g,h;
	map.clear();
	buildGraph(s, cutoffs, g);
	uint64_t label1 = hashGraph(g,true,map1);
	buildGraph(t, cutoffs, h);
	uint64_t label2 = hashGraph(h,true,map2);
	if (label1==label2) {
		for(auto &map_el: map1) imap1[map_el.second] = map_el.first; // C -> 1
		for(auto &map_el: map2) map[map_el.first] = imap1[map_el.second]; // 2 -> C -> 1
	}
};

virtual void isomorphicSelfMaps(AbstractSystem &s, std::list<std::map<int,int>> &maps) {
	graph_type g;
	buildGraph(s, cutoffs, g);
	vf2_map_builder<graph_type, graph_type> callback(g,g);
	callback.maps = &maps;
	boost::vf2_graph_iso(g,g,callback);
};

uint64_t hashGraph(graph_type &ga, bool canonical, std::map<int,int> &canonicalMap){
	//boost::timer::auto_cpu_timer t;
	canonicalMap.clear();

	//find the canonical labeling if necessary
	if(canonical) {
		canonicalMap=canonicalSort(ga);
	}
	else{
		//identity mapping
		auto vi= boost::vertices(ga);
		for(auto itv=vi.first; itv!=vi.second; itv++) {
			int uniqueIndex=ga[*itv].uniqueIndex;
			int remappedUniqueIndex=uniqueIndex;
			canonicalMap[uniqueIndex]=remappedUniqueIndex;
		}
	}


	//hash
	uint64_t hash=0;


	//hash the vertices
	auto vi= boost::vertices(ga);
	for(auto itv=vi.first; itv!=vi.second; itv++) {

		//create a 64 bits vertex index with the unique index and the color

		uint32_t i=uint32_t(canonicalMap[ga[*itv].uniqueIndex]);
		uint32_t c=uint32_t(ga[*itv].color);
		i=fmix32(i);
		c=fmix32(c);

		uint64_t l= i + (uint64_t(c) << 32);

		//hash it
		uint64_t lh=fmix64(l);
		//add it to the hash
		hash = hash ^ lh;
	}

	//hash the edges
	BGL_FORALL_EDGES_T(ed, ga, graph_type)
	{
		unsigned long v1, v2;
		v1 = source(ed,ga);
		v2 = target(ed,ga);


		uint32_t unique1=uint32_t(canonicalMap[ga[v1].uniqueIndex]);
		uint32_t unique2=uint32_t(canonicalMap[ga[v2].uniqueIndex]);

		unique1=fmix32(unique1);
		unique2=fmix32(unique2);

		//create a 64 bits edge index with the two vertices
		uint64_t min=uint64_t(std::min(unique1,unique2));
		uint64_t max=uint64_t(std::max(unique1,unique2));
		uint64_t l= max + (min << 32);


		//hash it
		uint64_t lh=fmix64(l);
		//add it to the hash
		hash = hash ^ lh;

	}


	return hash;

};

};


typedef boost::function< std::shared_ptr<AbstractStateLabeler>() > LabelerFactory;



const std::map<std::string, LabelerFactory> labelerFactory={
	//{"KMCStateLabeler", boost::factory<std::shared_ptr<KMCStateLabeler> >() },
	{"ConnectivityGraphStateLabeler", boost::factory<std::shared_ptr<ConnectivityGraphStateLabeler> >() }
};



#endif /* Graph_hpp */

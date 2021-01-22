//
//  LocalStore.cpp

//
//  Created by Danny Perez on 1/9/17.
//  Copyright © 2017 dp. All rights reserved.
//

#include "LocalStore.hpp"

#include "CustomTypedefs.hpp"
#include <boost/filesystem.hpp>
#include <boost/format.hpp>




/*
 * Definition of AbstractLocalDataStore
 */
AbstractLocalDataStore::AbstractLocalDataStore(){
	currentSize=0;
	maxSize=0;
};


unsigned int AbstractLocalDataStore::count(unsigned int dbKey, Label &key){
	return storedKeys[dbKey].count(key);
};

void AbstractLocalDataStore::touch(unsigned int dbKey, Label &key){
	if(count(dbKey,key)>0) {
		mru.touch(std::make_pair(dbKey,key));
	}
};

std::set<Label> AbstractLocalDataStore::availableKeys(unsigned int dbKey){
	if(storedKeys.count(dbKey)==0) {
		return std::set<Label>();
	}
	return storedKeys[dbKey];
};


std::set<unsigned int> AbstractLocalDataStore::availableDbKeys(){

	std::set<unsigned int> dbKeys;
	for(auto it=storedKeys.begin(); it!=storedKeys.end(); it++) {
		dbKeys.insert(it->first);
	}
	return dbKeys;
};

std::unordered_map<unsigned int, std::set<Label> > AbstractLocalDataStore::getStoredKeys(){
	return storedKeys;
};


void AbstractLocalDataStore::purge(){

	while(maxSize>0 and currentSize>maxSize) {
		auto p=mru.oldest();
		if(not p.first) {
			LOGGERA("AbstractLocalDataStore::purge ERROR: "<<currentSize<<" "<<maxSize<< " "<<p.second.first<<" "<<p.second.second)
			break;
		}
		erase(p.second.first,p.second.second);
	}

};



void AbstractLocalDataStore::setMaximumSize(unsigned long maxSize_){
	maxSize=maxSize_;
};






unsigned int PersistentLocalStore::count(unsigned int dbKey, Label &key){
	//return 0;
	return storedKeys[dbKey].count(key);
};

PersistentLocalStore::PersistentLocalStore(){
};

PersistentLocalStore::~PersistentLocalStore(){
	offsets.flush();
	io.flush();
	offsets.close();
	io.close();
};

int PersistentLocalStore::initialize(std::string homeDir, std::string baseName){

	std::string dbName=homeDir+"/"+baseName+".db";
	std::string offsetsName=homeDir+"/"+baseName+".offsets";
	bool restore=boost::filesystem::exists(dbName) && boost::filesystem::exists(offsetsName);

	boost::filesystem::path p(homeDir);
	boost::filesystem::create_directories( p.parent_path().string() );

	offsets.open(offsetsName, std::fstream::in | std::fstream::out | std::fstream::app);
	io.open(dbName, std::fstream::in | std::fstream::out | std::fstream::app | std::fstream::binary);

	std::cout<<"STORE INITIALIZE: "<<offsets.is_open()<<" "<<io.is_open()<<std::endl;

	if(restore) {
		std::cout<<"RESTORING DATABASE"<<std::endl;
		//read in the offsets
		int dbk;
		uint64_t k;
		long int off;
		long int s;
		offsets.seekg(0, offsets.beg);
		while(offsets>>dbk>>k>>off>>s) {
			LOGGER("RESTORING DB: "<<dbk<<" "<<k<<" "<<off<<" "<<s)
			locations[std::make_pair(dbk,k)]=std::pair<std::streampos,std::streamsize>(off,s);
			storedKeys[dbk].insert(k);
		}
		offsets.close();
		offsets.open(offsetsName, std::fstream::in | std::fstream::out | std::fstream::app);
	}
	return 0;
};

int PersistentLocalStore::put(unsigned int dbKey, uint64_t &key, RawDataVector &data){
	auto k=std::make_pair(dbKey,key);
	if(locations.count(k)==0) {
		//insert the item

		//move to the end of the stream
		io.seekp(0, io.end);
		locations[k]=std::pair<std::streampos,std::streamsize>(io.tellp(),data.size());
		offsets.seekp(0, offsets.end);
		offsets<<k.first<<" "<<k.second<<" "<<locations[k].first<<" "<<locations[k].second<<" "<<std::flush;

		io.write(&data[0],data.size());
		io.flush();

		storedKeys[dbKey].insert(key);


		/*
		   {
		        //TEST THE ROUND TRIP
		        RawDataVector dd;
		        get(dbKey,key,dd);
		        std::cout<<"TEST: "<<data.size()<<" "<<dd.size()<<" "<<(data==dd)<<std::endl;
		   }
		 */

	}
	return 0;
};


int PersistentLocalStore::get(unsigned int dbKey, uint64_t &key, RawDataVector &data){
	data.clear();
	auto k=std::make_pair(dbKey,key);
	if(locations.count(k)!=0) {
		auto loc=locations[k];
		io.seekg(loc.first, io.beg);
		data=RawDataVector(loc.second);
		io.read(&data[0],loc.second);
		io.sync();
	}
	else{
		return KEY_NOTFOUND;
	}
	return 0;
};

int PersistentLocalStore::createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet){
	return 0;
};

int PersistentLocalStore::sync(){
	offsets.flush();
	io.flush();
	return 0;
};

int PersistentLocalStore::erase(unsigned int dbKey, uint64_t &key){
	return 0;
};





#ifdef USE_SZ


CompressedPersistentLocalStore::CompressedPersistentLocalStore(){
	treeIndex=1;
	maxBufferSize=1;

	SZ_Init(NULL);
};

CompressedPersistentLocalStore::~CompressedPersistentLocalStore(){

};

int CompressedPersistentLocalStore::initialize(std::string homeDir, std::string baseName){

	compressedStore.initialize(homeDir,baseName+".c");
	uncompressedStore.initialize(homeDir,baseName+".u");
	treeStore.initialize(homeDir,baseName+".tree");

	storedKeys=compressedStore.getStoredKeys();

	return 0;
};

int CompressedPersistentLocalStore::put(unsigned int dbKey, uint64_t &key, RawDataVector &data){
	if(count(dbKey,key)==0) {
		//put in the buffer

		//unpack
		SystemType s;
		unpack(data,s);

		//extract compressible data
		std::vector<double> d;
		d=s.extractCompressibleData();

		auto k=std::make_pair(dbKey,key);
		buffer[k]=d;

		//pack and store the rest of the object
		RawDataVector rd;
		pack(rd,s);
		uncompressedStore.put(dbKey,key,rd);

		storedKeys[dbKey].insert(key);
	}

	//if buffer is full, compress and add to datastore
	if(buffer.size()>=maxBufferSize) {
		compressBuffer();
	}

	return 0;
};

void CompressedPersistentLocalStore::compressBuffer(){
	SZ_Init("sz.config");


	int nbBlocks=0;
	int blockSize=0;


	std::vector<double> rdata;
	std::vector<unsigned char> cdata;

	double realPrecision=0.0001;
	//construct compression buffer
	for(auto it=buffer.begin(); it!=buffer.end(); it++) {
		auto k=std::make_pair(it->first.first,it->first.second);
		blockSize=it->second.size();
		rdata.insert(rdata.end(),it->second.begin(),it->second.end());

	}

	{
		std::ofstream outd;
		outd.open("rd.out",std::ios::out);
		for(int i=0; i<rdata.size(); i++) {
			outd<<rdata[i]<<std::endl;
		}
		outd.close();
	}


	//generate the Huffman tree
	unsigned char* huffmanTree = NULL;
	unsigned int treeByteSize = 0;

	SZ_construct_HuffmanEncoder(SZ_DOUBLE, &(rdata[0]), buffer.size(), blockSize, realPrecision, &huffmanTree, &treeByteSize);

	//save the tree
	RawDataVector tree(huffmanTree,huffmanTree+treeByteSize);
	int ddd=0;
	treeStore.put(ddd,treeIndex,tree);


	int ii=0;
	for(auto it=buffer.begin(); it!=buffer.end(); it++,ii++) {
		auto k=std::make_pair(it->first.first,it->first.second);
		RawDataVector d;
		d.resize(sizeof(double)*blockSize);
		std::vector<unsigned char> cmprBuffer(sizeof(double)*blockSize);
		unsigned int cmprSize=0;
		SZ_compress_double_1D_exaalt_block(&(rdata[ii*blockSize]), blockSize, realPrecision, &(cmprBuffer[0]), &cmprSize);

		cmprBuffer.resize(cmprSize);
		LOGGER("PUT C: compression ratio "<<1.0*blockSize*4/cmprSize)
		pack(d,treeIndex,blockSize,cmprBuffer);
		compressedStore.put(k.first,k.second,d);

	}

	treeIndex++;
	buffer.clear();
	SZ_ReleaseHuffman();
	SZ_Finalize();
};


int CompressedPersistentLocalStore::get(unsigned int dbKey, uint64_t &key, RawDataVector &data){
	auto k=std::make_pair(dbKey,key);
	if(uncompressedStore.count(dbKey,key)==0 ) {
		return KEY_NOTFOUND;
	}

	RawDataVector d;
	uncompressedStore.get(dbKey,key,d);

	SystemType s;
	unpack(d,s);
	std::vector<double> dd;
	std::vector<unsigned char> dc;
	if( buffer.count(k)>0 ) {
		dd=buffer[k];
	}
	else{
		SZ_Init("sz.config");
		RawDataVector dt;
		Label treeLabel;
		int blockSize;
		double realPrecision2 = 0;
		compressedStore.get(dbKey,key,dt);
		unpack(dt,treeLabel,blockSize,dc);
		dd.resize(blockSize);


		RawDataVector tdata;
		int ddd=0;
		treeStore.get(ddd,treeLabel,tdata);

		unsigned char *p=(unsigned char *) &(tdata[0]);
		node huffmanTreeRootNode = SZ_reconstruct_HuffmanEncoder(p, tdata.size(), &realPrecision2);


		size_t cmprBlockSize;
		//decompress
		p=(unsigned char *) &(dc[0]);
		SZ_decompress_double_1D_exaalt_block(huffmanTreeRootNode, realPrecision2, p, dc.size(), &(dd[0]), blockSize);
		SZ_ReleaseHuffman();
		SZ_Finalize();
	}

	//reassemble item
	s.restoreCompressibleData(dd);
	pack(data,s);


	return 0;
};

int CompressedPersistentLocalStore::createDatabase(unsigned int dbKey, bool allowDuplicates, bool eraseOnGet){
	return 0;
};

int CompressedPersistentLocalStore::sync(){
	compressedStore.sync();
	uncompressedStore.sync();
	treeStore.sync();
	return 0;
};

int CompressedPersistentLocalStore::erase(unsigned int dbKey, uint64_t &key){
	return 0;
};
#endif






/*
 * Definition of STLLocalDataStore
 */


STLLocalDataStore::STLLocalDataStore(){
	currentSize=0;
};

int STLLocalDataStore::initialize(std::string homeDir_, std::string baseName_){
	return 0;
};

int STLLocalDataStore::put(unsigned int dbKey, Label &key, RawDataVector &data){

	if( dbm.count(dbKey)==0 ) {
		//*return DBKEY_NOTFOUND;
		dbAttributes[dbKey]=std::make_pair(true, true);

	}

	bool allowDuplicates=dbAttributes[dbKey].first;

	if(allowDuplicates or count(dbKey, key)==0) {
		dbm[dbKey].insert(std::make_pair(key,data));
		storedKeys[dbKey].insert(key);
		currentSize+=data.size();
	}

	mru.touch(std::make_pair(dbKey,key));

	return 0;

};

int STLLocalDataStore::get(unsigned int dbKey, Label &key, RawDataVector &data){

	if( dbm.count(dbKey)==0 ) {
		return DBKEY_NOTFOUND;
	}

	if(dbm[dbKey].count(key)==0) {
		return KEY_NOTFOUND;
	}

	auto it=dbm[dbKey].find(key);
	data=it->second;

	bool eraseOnGet=dbAttributes[dbKey].second;
	if(eraseOnGet) {
		int size;
		size=data.size();
		currentSize-=size;
		dbm[dbKey].erase(it);
	}
	if( dbm[dbKey].count(key)==0) {
		mru.erase(std::make_pair(dbKey,key));
		storedKeys[dbKey].erase(key);
	}
	else{
		mru.touch(std::make_pair(dbKey,key));
	}

	return 0;
};


int STLLocalDataStore::createDatabase(unsigned int dbKey, bool allowDuplicates_, bool eraseOnGet_){
	if(dbm.count(dbKey)==0) {
		dbm[dbKey]=std::unordered_multimap< Label, RawDataVector>();
		dbAttributes[dbKey]=std::make_pair(allowDuplicates_, eraseOnGet_);
	}
	return 0;
};

int STLLocalDataStore::sync(){
	//This does nothing for now, but we could serialize if needed
	return 0;
};


int STLLocalDataStore::erase(unsigned int dbKey, Label &key){
	int size=0;
	if(dbm[dbKey].count(key)!=0) {

		auto it=dbm[dbKey].find(key);
		size=it->second.size();
		dbm[dbKey].erase(it);

		if(dbm[dbKey].count(key)==0) {
			mru.erase(std::make_pair(dbKey,key));
			storedKeys[dbKey].erase(key);
		}
	}
	else{
		std::cout<<"STLLocalDataStore::erase ERROR: no such item "<<dbKey<<" "<<key<<std::endl;
	}
	currentSize-=size;

	return size;
};

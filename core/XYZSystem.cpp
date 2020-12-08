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



#include "XYZSystem.hpp"

#include <cassert>
#include <map>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>


XYZSystem::XYZSystem(){
    nAtoms=0;
    
    for(int i=0;i<NDIM;i++){
        p[i]=1;
        for(int j=0;j<NDIM;j++){
            b[i][j]=0;
        }
    }
};

void XYZSystem::clearAtoms(){
    spe.clear();
    pos.clear();
    id.clear();
    
    nAtoms=0;
};

int XYZSystem::getNAtoms(){
    return nAtoms;
};

void XYZSystem::setNAtoms(int n){
    clearAtoms();
    nAtoms=n;
    spe=std::vector<int>(nAtoms,0);
    pos=std::vector<double>(NDIM*nAtoms,0);
    id=std::vector<int>(nAtoms,0);
    
};

void XYZSystem::setPosition(int iAtom, int iDim, double p){
    int indice=iAtom*NDIM+iDim;
    assert(indice<nAtoms*NDIM);
    pos[indice]=p;
};

double XYZSystem::getPosition(int iAtom, int iDim){
    int indice=iAtom*NDIM+iDim;
    assert(indice<nAtoms*NDIM);
    return pos[indice];
};

void XYZSystem::setSpecies(int iAtom, int species){
    assert(iAtom<nAtoms);
    spe[iAtom]=species;
};

int XYZSystem::getSpecies(int iAtom){
    assert(iAtom<nAtoms);
    return spe[iAtom];
};

void XYZSystem::setUniqueID(int iAtom, int i){
    assert(iAtom<nAtoms);
    id[iAtom]=i;
};

int XYZSystem::getUniqueID(int iAtom){
    assert(iAtom<nAtoms);
    return id[iAtom];
};

double XYZSystem::getBox(int iDim, int jDim){
    return b[iDim][jDim];
};

void XYZSystem::setBox(int iDim, int jDim, double box){
    b[iDim][jDim]=box;
};

int XYZSystem::getPeriodic(int iDim){
    return p[iDim];
};

void XYZSystem::setPeriodic(int iDim, bool periodic){
    p[iDim]=(periodic ? 1 : 0);
};

void XYZSystem::pack(std::vector<char> &v){
    boost::iostreams::back_insert_device<std::vector<char>> sink{v};
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> os{sink};
    boost::archive::binary_oarchive oa(os, boost::archive::no_header);
    oa << *this;
    os.flush();
};

void XYZSystem::unpack(std::vector<char> &v){
    clearAtoms();
    boost::iostreams::array_source source{v.data(), v.size()};
    boost::iostreams::stream<boost::iostreams::array_source> is{source};
    boost::archive::binary_iarchive ia(is, boost::archive::no_header);
    ia >> *this;
};


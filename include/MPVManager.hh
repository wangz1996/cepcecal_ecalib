#ifndef MPVMANAGER_HH
#define MPVMANAGER_HH
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

class MPVManager{
public:
    MPVManager(const std::string& fname,const std::string& tname);
    ~MPVManager();
    void readMPV(const std::string& fname,const std::string& tname);
    float getMPV(int ch){return umap_mpv[ch];}

private:
    TFile *f;
    TTree *t;
    std::unordered_map<int,float> umap_mpv;
};

#endif
#include "MPVManager.hh"

MPVManager::MPVManager(const std::string& fname,const std::string& tname){
    this->readMPV(fname,tname);
}
MPVManager::~MPVManager(){
    
}

void MPVManager::readMPV(const std::string& fname,const std::string& tname){
    f = TFile::Open(TString(fname),"READ");
    t = (TTree*)f->Get(TString(tname));
    double LMPV;
    int CellID;
    t->SetBranchAddress("CellID",&CellID);
    t->SetBranchAddress("LandauMPV",&LMPV);
    for(int i=0;i<t->GetEntries();i++){
        t->GetEntry(i);
        //std::cout<<CellID<<" "<<LMPV<<std::endl;
        umap_mpv[CellID] = LMPV;
    }
}
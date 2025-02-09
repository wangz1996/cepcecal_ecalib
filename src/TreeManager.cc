#include "TreeManager.hh"

TreeManager::TreeManager() {
    wcm = new WCManager();
}
TreeManager::~TreeManager() {
    
}

void TreeManager::readTree(const std::string &filename, const std::string &treename) {
    f = TFile::Open(TString(filename),"READ");
    t = (TTree*)f->Get(TString(treename));
    //Set Branch Address
    t->SetBranchAddress("CellID",&CellID);
    t->SetBranchAddress("Hit_X",&Hit_X);
    t->SetBranchAddress("Hit_Y",&Hit_Y);
    t->SetBranchAddress("Hit_Z",&Hit_Z);
    t->SetBranchAddress("Hit_Energy",&Hit_E);
    t->SetBranchAddress("Hit_ADC",&Hit_ADC);
}

void TreeManager::eventLoop(){
    for(size_t ientry=0; ientry<t->GetEntries(); ientry++){
        t->GetEntry(ientry);
        //clear the wcm
        wcm->clear();
        for(size_t i=0; i<Hit_E->size(); i++){
            int layer = CellID->at(i)/100000;
            wcm->fillHit(layer,Hit_X->at(i),Hit_Y->at(i),Hit_E->at(i));
        }
        //Find the cell which wc falls into
        for(size_t iLayer=0; iLayer<NLayers; iLayer++){
            int cellid=0;
            if(!wcm->layerNoEmpty(iLayer))continue;
            if(iLayer%2==1){//odd layers
                int x_i = round((wcm->getWCX(iLayer)+108.65)/5.3);
                int y_i = round((wcm->getWCY(iLayer)+90.8)/45.4);
                if(y_i!=2){continue;}
                int chnid = posArr[x_i];
                cellid = iLayer*100000+2*10000+chnid;
            }
            else{//even layers
                int x_i = round((wcm->getWCX(iLayer)+90.8)/45.4);
                int y_i = round((wcm->getWCY(iLayer)+108.65)/5.3);
                if(x_i!=2){continue;}
                int chnid = posArr[y_i];
                cellid = iLayer*100000+2*10000+chnid;
            }
            
        }
    }
}
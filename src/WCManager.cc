#include "WCManager.hh"

WCManager::WCManager(){

}
WCManager::~WCManager(){
    
}

void WCManager::fillHit(const int& cellid, const float& x, const float& y, const float& E, const float& ADC){
    int layer = cellid/100000;
    if(layer < 0 || layer >= NLayers) return;
    SumEX[layer] += x*E;
    SumEY[layer] += y*E;
    SumE[layer] += E;
    umap_id_e[cellid] = E;
    umap_id_adc[cellid] = ADC;
}


float WCManager::getWCX(const int& layer){
    if(layer < 0 || layer >= NLayers){
        std::cout<<"Wrong layer input"<<std::endl;
        throw(layer);
    }
    if(SumE[layer] == 0){
        std::cout<<"SumE is zero"<<std::endl;
        throw(layer);  
    }
    return SumEX[layer]/SumE[layer];
}

float WCManager::getWCY(const int& layer){
    if(layer < 0 || layer >= NLayers){
        std::cout<<"Wrong layer input"<<std::endl;
        throw(layer);
    }
    if(SumE[layer] == 0){
        std::cout<<"SumE is zero"<<std::endl;
        throw(layer);  
    }
    return SumEY[layer]/SumE[layer];
}

void WCManager::clear(){
    SumEX={};
    SumEY={};
    SumE={};
    umap_id_e.clear();
}
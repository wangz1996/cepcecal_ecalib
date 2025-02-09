#include "WCManager.hh"

WCManager::WCManager(){

}
WCManager::~WCManager(){
    
}

void WCManager::fillHit(const int& layer, const float& x, const float& y, const float& E){
    if(layer < 0 || layer >= NLayers) return;
    SumEX[layer] += x*E;
    SumEY[layer] += y*E;
    SumE[layer] += E;
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
}
#ifndef WCManager_HH
#define WCManager_HH

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <array>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include "TH1D.h"

class WCManager{

public:
  WCManager();
  ~WCManager();

  void fillHit(const int& cellid, const float& x, const float& y, const float& E,const float& ADC);
  float getWCX(const int& layer);
  float getWCY(const int& layer);
  float getE(const int& cellid){return umap_id_e[cellid];}
  float getADC(const int& cellid){return umap_id_adc[cellid];}
  bool layerNoEmpty(const int& layer){return SumE[layer]>0;}
  void clear();

private:
  //Weighted center Sum_E*x/Sum_E
  std::unordered_map<int,float> umap_id_e;
  std::unordered_map<int,float> umap_id_adc;
  static constexpr int NLayers = 30;
  std::array<float,NLayers> SumEX;
  std::array<float,NLayers> SumEY;
  std::array<float,NLayers> SumE;
};

#endif
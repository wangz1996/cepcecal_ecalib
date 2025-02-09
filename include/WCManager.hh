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

  void fillHit(const int& layer, const float& x, const float& y, const float& E);
  float getWCX(const int& layer);
  float getWCY(const int& layer);
  bool layerNoEmpty(const int& layer){return SumE[layer]>0;}
  void clear();

private:
  //Weighted center Sum_E*x/Sum_E
  static constexpr int NLayers = 30;
  std::array<float,NLayers> SumEX;
  std::array<float,NLayers> SumEY;
  std::array<float,NLayers> SumE;
};

#endif
#ifndef SFMANAGER_HH
#define SFMANAGER_HH

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TMinuit.h"
#include "TSystem.h"

class SFManager{
public:
  SFManager();
  ~SFManager();

  std::vector<double> getOriginalValues(TH1D* hist);
  // static void fcnForMinuit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  float getSF(TH1D* h1,TH1D* h2);
  void recordSF(const int& id, const float& sf);
  void saveSF();

  static TH1D *g_h2;
  static std::vector<double> g_originalValues;
  
private:
  int ierflg = 0;
  static constexpr double vstart = 1.0;
  static constexpr double step   = 0.01;
  static constexpr double minVal = 0.1;
  static constexpr double maxVal = 3.0;
  std::unordered_map<int,float> umap_id_sf;
  float sf;
  int cellid;
  
};



#endif
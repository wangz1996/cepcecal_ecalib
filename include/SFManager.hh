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

class SFManager{
public:
  SFManager();
  ~SFManager();

  std::vector<double> getOriginalValues(TH1D* hist);
  void fcnForMinuit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  float getSF(TH1D* h1,TH1D* h2);

private:
  int ierflg = 0;
  static double vstart = 1.0;   
  static double step   = 0.01;  
  static minVal = 0.1;
  static maxVal = 5.0;
  std::vector<double> *g_originalValues=nullptr;
  std::vector<double> *values=nullptr;
}

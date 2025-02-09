#ifndef TREEMANAGER_HH
#define TREEMANAGER_HH

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"

#include "WCManager.hh"

class TreeManager{

public:
  TreeManager();
  ~TreeManager();

  void readTree(const std::string &filename, const std::string &treename);
  void eventLoop();

private:
  std::vector<int>* CellID=nullptr;
  std::vector<double>* Hit_X=nullptr;
  std::vector<double>* Hit_Y=nullptr;
  std::vector<double>* Hit_Z=nullptr;
  std::vector<double>* Hit_E=nullptr;
  std::vector<double>* Hit_ADC=nullptr;
  static constexpr int NLayers = 30;

  TTree *t;
  TFile *f;

  WCManager *wcm;

  static constexpr std::array<int,30> posArr{35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,14,13,12,11,10,9,8,7,6,5,4,3};
};

#endif
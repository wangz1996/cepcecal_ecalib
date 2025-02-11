#ifndef TREEMANAGER_HH
#define TREEMANAGER_HH

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <filesystem>
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include <ROOT/RDataFrame.hxx>
#include "WCManager.hh"
#include "SFManager.hh"

class TreeManager{

public:
  TreeManager();
  ~TreeManager();

  void readTree(const std::string &filename, const std::string &treename);
  void eventLoop();
  void setDoScale(const std::string& _sf_path){do_scale=true;sf_path=_sf_path;}

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
  SFManager *sfm;

  static constexpr std::array<int,30> posArr{35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,14,13,12,11,10,9,8,7,6,5,4,3};
  std::unordered_map<int,TH1D*> umap_id_hist;
  std::unordered_map<int,TH1D*> umap_id_hscaled;
  std::unordered_map<int,float> umap_id_sf;

  bool do_scale;
  std::string sf_path;

  const std::set<int> oddSet{7,8,9,10,11,12,13,14,18,19};
  const std::set<int> evenSet{7,8,9,10,11,12,13,14};

  std::unordered_map<int, int> stdcell{
      {4, 12},
      {5, 10},
      {8, 11}, // 14
      {9, 18},
      {6, 9},
      {7, 13},
      {10, 11},
      {11, 11}, // 8
      {12, 10},
      {13, 12},
      {14, 9},
      {15, 10},
      {16, 11},
      {17, 12}};
};

#endif
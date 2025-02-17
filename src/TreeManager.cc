#include "TreeManager.hh"

TreeManager::TreeManager()
{
    wcm = new WCManager();
    sfm = new SFManager();
    // fout = new TFile("hist.root","RECREATE");
}
TreeManager::~TreeManager()
{
    // fout->Close();
}

void TreeManager::readTree(const std::string &filename, const std::string &treename)
{
    f = TFile::Open(TString(filename), "READ");
    t = (TTree *)f->Get(TString(treename));
    // Set Branch Address
    t->SetBranchAddress("CellID", &CellID);
    t->SetBranchAddress("Hit_X", &Hit_X);
    t->SetBranchAddress("Hit_Y", &Hit_Y);
    t->SetBranchAddress("Hit_Z", &Hit_Z);
    t->SetBranchAddress("Hit_Energy", &Hit_E);
    t->SetBranchAddress("Hit_ADC", &Hit_ADC);

    if (do_scale)
    {
        ROOT::RDataFrame df("sf", sf_path);
        df.Foreach([this](const int &id, const float &sf)
                   { this->umap_id_sf[id] = sf; }, {"cellid", "sf"});
    }
}

void TreeManager::eventLoop()
{
    for (size_t ientry = 0; ientry < t->GetEntries(); ientry++)
    {
        t->GetEntry(ientry);
        // clear the wcm
        wcm->clear();
        for (size_t i = 0; i < Hit_E->size(); i++)
        {
            // int layer = CellID->at(i)/100000;
            wcm->fillHit(CellID->at(i), Hit_X->at(i), Hit_Y->at(i), Hit_E->at(i), Hit_ADC->at(i));
        }
        // Find the cell which wc falls into
        for (size_t iLayer = 0; iLayer < NLayers; iLayer++)
        {
            int cellid = 0;
            if (!wcm->layerNoEmpty(iLayer))
                continue;
            if (iLayer % 2 == 1)
            { // odd layers
                int x_i_odd = round((wcm->getWCX(iLayer) + 108.65) / 5.3);
                int y_i_odd = round((wcm->getWCY(iLayer) + 90.8) / 45.4);
                if (y_i_odd != 2)
                {
                    continue;
                }
                int chnid = posArr[x_i_odd];
                if (oddSet.count(chnid) == 0)
                {
                    continue;
                }
                cellid = iLayer * 100000 + 2 * 10000 + chnid;
                if (do_sum and (chnid == 11 or chnid == 12))
                {
                    sum_ADC = 0;
                    sum_Energy = 0;
                    sum_SADC = 0;
                    sum_SEnergy = 0;
                    for (auto j : oddSet)
                    {
                        center_ID = iLayer * 100000 + 2 * 10000 + j;
                        if (umap_id_sf.count(center_ID) == 0)
                        {
                            umap_id_sf[center_ID] = 1.;
                        }
                        double value_toSADC_add = wcm->getADC(center_ID) / 1000.0 * umap_id_sf[center_ID];
                        double value_toSE_add = wcm->getE(center_ID)  * umap_id_sf[center_ID];
                      
                        double value_toADC_add = wcm->getADC(center_ID) / 1000.0;
                        double value_toE_add = wcm->getE(center_ID) ;
                     
                        sum_ADC += value_toADC_add;
                        sum_Energy += value_toE_add;
                        sum_SADC += value_toSADC_add;
                        sum_SEnergy += value_toSE_add;
                    }
                }
            }
            else
            { // even layers
                int x_i_even = round((wcm->getWCX(iLayer) + 90.8) / 45.4);
                int y_i_even = round((wcm->getWCY(iLayer) + 108.65) / 5.3);
                if (x_i_even != 2)
                {
                    continue;
                }
                int chnid = posArr[y_i_even];
                if (evenSet.count(chnid) == 0)
                {
                    continue;
                }
                cellid = iLayer * 100000 + 2 * 10000 + chnid;
                if (do_sum and (chnid == 9 or chnid == 10))
                {
                    sum_ADC = 0;
                    sum_Energy = 0;
                    sum_SADC = 0;
                    sum_SEnergy = 0;
                    for (auto j : evenSet)
                    {
                        center_ID = iLayer * 100000 + 2 * 10000 + j;
                        if (umap_id_sf.count(center_ID) == 0)
                        {
                            umap_id_sf[center_ID] = 1.;
                        }
                        
                        double value_toSADC_add = wcm->getADC(center_ID) / 1000.0 * umap_id_sf[center_ID];
                        double value_toSE_add = wcm->getE(center_ID)  * umap_id_sf[center_ID];
                      
                        double value_toADC_add = wcm->getADC(center_ID) / 1000.0;
                        double value_toE_add = wcm->getE(center_ID) ;
                     
                        sum_ADC += value_toADC_add;
                        sum_Energy += value_toE_add;
                        sum_SADC += value_toSADC_add;
                        sum_SEnergy += value_toSE_add;
                    }
                }
            }
            // Now we know the cellid of the WC
            // We can fill the histogram now
            if (umap_id_hist.count(cellid) == 0)
            {
                umap_id_hist[cellid] = new TH1D(TString::Format("h%d", cellid), "", 47, 3, 50);
            }
            if (do_scale && umap_id_hscaled.count(cellid) == 0)
            {
                umap_id_hscaled[cellid] = new TH1D(TString::Format("hs%d", cellid), "", 47, 3, 50);
            }
            umap_id_hist[cellid]->Fill(wcm->getADC(cellid) / 1000.);
            if (do_scale)
            {
                umap_id_sf[cellid] = umap_id_sf.count(cellid) == 0 ? 1. : umap_id_sf[cellid];
                umap_id_hscaled[cellid]->Fill(wcm->getADC(cellid) / 1000. * umap_id_sf[cellid]);
            }
            if (do_sum && (umap_sum_SADC.count(iLayer) == 0 or umap_sum_ADC.count(iLayer) == 0 or umap_sum_Energy.count(iLayer) == 0 or umap_sum_SEnergy.count(iLayer) == 0))
            {
                umap_sum_SADC[iLayer] = new TH1D(TString::Format("hsumSADC%d", iLayer), "", 47, 3, 70);
                umap_sum_ADC[iLayer] = new TH1D(TString::Format("hsumADC%d", iLayer), "", 47, 3, 70);   
                umap_sum_Energy[iLayer] = new TH1D(TString::Format("hsumEnergy%d", iLayer), "", 47, 3, 180);
                umap_sum_SEnergy[iLayer] = new TH1D(TString::Format("hsumSEnergy%d", iLayer), "", 47, 3, 180);
            }
            if (do_sum)
            {
                umap_sum_ADC[iLayer]->Fill(sum_ADC);
                umap_sum_Energy[iLayer]->Fill(sum_Energy);
                if (do_scale)
                {
                    umap_sum_SADC[iLayer]->Fill(sum_SADC);
                    umap_sum_SEnergy[iLayer]->Fill(sum_SEnergy);
                }
            }
        }
    }

    // Calculate scale factors
    for (auto &[id, hist] : umap_id_hist)
    {
        auto layer = id / 100000;
        if (stdcell.count(layer) == 0)
            continue;
        auto stdcellid = layer * 100000 + 20000 + stdcell[layer];
        auto sf = sfm->getSF(hist, umap_id_hist[stdcellid]);
        if (id == 620010)
            std::cout << sf << std::endl;
        sfm->recordSF(id, sf);
    }
    sfm->saveSF();
    // Write the histograms to a root file
    std::cout << "Starting to save histograms..." << std::endl;
    auto fout = new TFile("hist.root", "RECREATE");
    fout->cd();
    for (size_t i = 0; i < NLayers; i++)
    {
        fout->mkdir(TString::Format("layer%d", i));
    }
    for (auto &[id, hist] : umap_id_hist)
    {
        int layer = id / 100000;
        fout->cd(TString::Format("layer%d", layer));
        hist->Write();
    }

    if (do_scale)
    {
        for (auto &[id, hist] : umap_id_hscaled)
        {
            int layer = id / 100000;
            fout->cd(TString::Format("layer%d", layer));
            hist->Write();
        }
    }
    if (do_sum )
    {
        for (auto &[id, hist] : umap_sum_SADC)
        {
            int layer = id ;
            fout->cd(TString::Format("layer%d", layer));
            hist->Write();
        }
        for (auto &[id, hist] : umap_sum_ADC)
        {
            int layer = id ;
            fout->cd(TString::Format("layer%d", layer));
            hist->Write();
        }
        for (auto &[id, hist] : umap_sum_SEnergy)
        {
            int layer = id ;
            fout->cd(TString::Format("layer%d", layer));
            hist->Write();
        }
        for (auto &[id, hist] : umap_sum_Energy)
        {
            int layer = id ;
            fout->cd(TString::Format("layer%d", layer));
            hist->Write();
        }
    }
    fout->Close();
}
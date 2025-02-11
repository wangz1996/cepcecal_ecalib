#include "SFManager.hh"
TH1D *SFManager::g_h2 = nullptr;
std::vector<double> SFManager::g_originalValues;
SFManager::SFManager(){
    gSystem->Load("libMinuit");
}
SFManager::~SFManager(){
    
}
std::vector<double> SFManager::getOriginalValues(TH1D* hist){
    std::vector<double> values;
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        double binContent = hist->GetBinContent(bin); 
        double binCenter = hist->GetBinCenter(bin);   
        for (int i = 0; i < static_cast<int>(binContent); ++i) {
            values.push_back(binCenter);
        }
    }
    return values;
}
void fcnForMinuit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    double scaleFactor = par[0];
    TH1D h_temp("h_temp", "Scaled h1 distribution", 
        SFManager::g_h2->GetNbinsX(), 
        SFManager::g_h2->GetXaxis()->GetXmin(), 
        SFManager::g_h2->GetXaxis()->GetXmax());
    for (double val : SFManager::g_originalValues) {
            h_temp.Fill(val * scaleFactor);
        }
    double integral = h_temp.Integral();
    if (integral > 1e-12) { 
            h_temp.Scale(1.0 / integral);
    } else {
         f = 1e12;
         return;
        }
    double chi2 = 0.0;
    for (int bin = 1; bin <= SFManager::g_h2->GetNbinsX(); ++bin) {
            double h2_val = SFManager::g_h2->GetBinContent(bin);
            double h1_val = h_temp.GetBinContent(bin);
        if (h2_val > 0) {
                double diff = (h1_val - h2_val);
                chi2 += ((diff * diff) / h2_val);
        }
    }
    f = chi2;
}
float SFManager::getSF(TH1D* h1,TH1D* h2){
	TH1D *g_h1= (TH1D*)h1->Clone();
	SFManager::g_h2= (TH1D*)h2->Clone();
    SFManager::g_originalValues = getOriginalValues(g_h1);
    SFManager::g_h2->Scale(1.0/SFManager::g_h2->Integral());
    TMinuit *gMinuit = new TMinuit(1);
    gMinuit->SetFCN(fcnForMinuit);
	gMinuit->SetPrintLevel(-1);
    int ierflg = 0;
    gMinuit->mnparm(0, "scaleFactor", vstart, step, minVal, maxVal, ierflg);
	gMinuit->SetPrintLevel(-1);
    gMinuit->Migrad();

    double best_sf, best_sf_err;
    gMinuit->GetParameter(0, best_sf, best_sf_err);
	return best_sf;
}

void SFManager::recordSF(const int& id,const float& sf){
    umap_id_sf[id]=sf;
}
void SFManager::saveSF(){
    // Save the scale factor to a file
    auto fout=new TFile("sf_40gev.root","RECREATE");
    auto tout=new TTree("sf","sf");
    tout->Branch("cellid",&cellid,"cellid/I");
    tout->Branch("sf",&sf,"sf/F");
    for(auto& [key,val]:umap_id_sf){
        cellid=key;
        sf=val;
        tout->Fill();
    }
    fout->cd();
    tout->Write();
    fout->Close();
}
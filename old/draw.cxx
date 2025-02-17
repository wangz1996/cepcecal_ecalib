#include "standard.hh"
TH1D *g_h2;
std::vector<double> g_originalValues;
const std::set<int> smooth_layer{6,7,10,13,14,15,17};
std::unordered_map<int,float> readsf(){
	std::unordered_map<int,float> umap_sf;
	ROOT::RDataFrame df("sf","sf_backup.root");
	df.Foreach([&umap_sf](int id,float sf){
		umap_sf.insert(std::pair<int,float>(id,sf));
	},{"cellId","sf"});
	return umap_sf;
}
float getMaximumBinCenter(TH1D *h,float thred=2){
	std::vector<std::array<float,2>> vec;
	for(int i=0;i<h->GetNbinsX();i++){
		vec.emplace_back(std::array<float,2>{float(h->GetBinCenter(i+1)),float(h->GetBinContent(i+1))});
	}
	std::sort(vec.begin(),vec.end(),[](std::array<float,2> a,std::array<float,2> b){return a[1]>b[1];});
	for(auto i:vec){
		if(i[0]>thred){
			return i[0];
		}
	}
	return 0.;
}
void shower_profile(){
	TFile *fdata=TFile::Open("hist_all.root","READ");
	TFile *fmc=TFile::Open("MC/hist_all.root","READ");
	TH1D *hdata=(TH1D*)fdata->Get("layer5/h520011");
	hdata->Scale(1./hdata->Integral());
	TH1D *hmc=(TH1D*)fmc->Get("layer5/h520011");
	hmc->Scale(1./hmc->Integral());

	TCanvas *c1=new TCanvas("c1","test",800,600);
	c1->cd();
	hdata->SetLineColor(kRed);
	hdata->Draw();
	hmc->SetLineColor(kBlue);
	hmc->Draw("same");
	c1->SaveAs("profile_5.png");

}
void fcnForMinuit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    // 只有 1 个自由参数：scaleFactor
    double scaleFactor = par[0];

    // 创建一个临时 hist，用来填充 '缩放后' 的数据分布
    // 假设 g_h2 有 N 个 bin，范围和 bin 宽度都跟 h1/h2 一样
    // 你也可以直接用和 h1 相同的 binning
    TH1D h_temp("h_temp", "Scaled h1 distribution", 
                g_h2->GetNbinsX(), 
                g_h2->GetXaxis()->GetXmin(), 
                g_h2->GetXaxis()->GetXmax());

    // 将 g_originalValues 中的每个值乘以 scaleFactor 后填入临时 hist
    for (double val : g_originalValues) {
        h_temp.Fill(val * scaleFactor);
    }

    // 归一化
    double integral = h_temp.Integral();
    if (integral > 1e-12) { // 防止除零
        h_temp.Scale(1.0 / integral);
    } else {
        // 如果意外出现 integral=0，通常是非常态情形；我们可以给一个很大的 chi2 以免误导
        f = 1e12;
        return;
    }
    
    // 计算 chi2
    double chi2 = 0.0;
    for (int bin = 1; bin <= g_h2->GetNbinsX(); ++bin) {
        double h2_val = g_h2->GetBinContent(bin);
        double h1_val = h_temp.GetBinContent(bin);
        if (h2_val > 0) {
            double diff = (h1_val - h2_val);
            chi2 += ((diff * diff) / h2_val);
        }
    }
    //use difference of maximum x value as chi2
    // chi2 = abs(h_temp.GetBinCenter(h_temp.GetMaximumBin()) - g_h2->GetBinCenter(g_h2->GetMaximumBin()));
    // std::cout<<par[0]<<" "<<h_temp.GetBinCenter(h_temp.GetMaximumBin())<<" "<<g_h2->GetBinCenter(g_h2->GetMaximumBin())<<" "<<chi2<<std::endl;
    // 最终把 chi2 写到 f 中，Minuit 就会拿这个作为“目标函数”
    f = chi2;
}
std::vector<double> getOriginalValues(TH1D* hist) {
    std::vector<double> values;

    // 遍历所有 bin
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        double binContent = hist->GetBinContent(bin); // 获取 bin 的填充值
        double binCenter = hist->GetBinCenter(bin);   // 获取 bin 的中心值（原始值）

        // 如果 bin 有填充的数值，按次数添加到 vector
        for (int i = 0; i < static_cast<int>(binContent); ++i) {
            values.push_back(binCenter);
        }
    }

    return values;
}
float getSF(TH1D* h1,TH1D* h2){
	TH1D *g_h1= (TH1D*)h1->Clone();
	g_h2= (TH1D*)h2->Clone();
    g_originalValues = getOriginalValues(g_h1);
    g_h2->Scale(1.0/g_h2->Integral());
    TMinuit *gMinuit = new TMinuit(1);
    gMinuit->SetFCN(fcnForMinuit);
	gMinuit->SetPrintLevel(-1);
    int ierflg = 0;
    double vstart = 1.0;   // 初始值
    double step   = 0.01;  // 步长
    double minVal = 0.1;
    double maxVal = 5.0;
    gMinuit->mnparm(0, "scaleFactor", vstart, step, minVal, maxVal, ierflg);
	gMinuit->SetPrintLevel(-1);
    gMinuit->Migrad();

    double best_sf, best_sf_err;
    gMinuit->GetParameter(0, best_sf, best_sf_err);
	return best_sf;
}
std::string makeKey(const float& x,const float& y,const float& z){
	return std::to_string(int(round(x+1e-8)))+"_"+std::to_string(int(round(y+1e-8)))+"_"+std::to_string(int(round(z+1e-8)));
}
void shift(const std::string& fpath,const std::string& tname,bool applysf){
	std::unordered_map<int,float> sfs;
	if(applysf)sfs = readsf();
	TFile *fsf = new TFile("sf_40gev.root","RECREATE");
	TTree *t=new TTree("sf","sf");
	TString fout_name = applysf ? "hist_all_energy_scaled.root" : "hist_all_energy.root";
	TFile *fout=new TFile(fout_name,"RECREATE");
	int _cellid=0;
	float _sf=0.;
	t->Branch("cellId",&_cellid,"cellId/I");
	t->Branch("sf",&_sf,"sf/F");
	ROOT::RDataFrame df(tname,fpath);
	const float beam_x=-2.;
	const float beam_y=9.3;
	const float beam_xsig=10.;
	const float beam_ysig=10.;
	//Save map from round x y z to cellid
	std::unordered_map<std::string,int> umap_xyz_id;
	df.Foreach([&umap_xyz_id](const std::vector<int>& id,const std::vector<double>& Hit_X,const std::vector<double>& Hit_Y,const std::vector<double>& Hit_Z){
		for(size_t i=0;i<Hit_X.size();i++){
			umap_xyz_id[makeKey(Hit_X[i],Hit_Y[i],Hit_Z[i])] = id[i];
		}
	},{"CellID","Hit_X","Hit_Y","Hit_Z"});
	std::unordered_map<int,TH1D*> umap_id_th1d;
	//create a umap to store hist of each layer with center cell hit
	std::unordered_map<int,TH1D*> umap_layer_th1d;
	//Channels need to be summed for odd layers
	std::array<int,10> odd_chns={7,8,9,10,11,12,13,14,18,19};
	//Channels need to be summed for even layers
	std::array<int,8> even_chns={7,8,9,10,11,12,13,14};
	const int N=30;
	df.Foreach([&umap_id_th1d,&umap_xyz_id,&sfs,applysf,&umap_layer_th1d,odd_chns,even_chns](const std::vector<int>& id,const std::vector<double>& Hit_X,const std::vector<double>& Hit_Y,const std::vector<double>& Hit_Z,const std::vector<double>& Hit_E,const std::vector<double>& Hit_ADC){
		std::array<float,30> weight_per_layer={};
		std::array<float,30> wsx_per_layer={};
		std::array<float,N> sx_per_layer={};
		std::array<float,30> wsy_per_layer={};
		std::array<float,N> sy_per_layer={};
		std::array<float,N> z_per_layer={};
		std::unordered_map<int,float> umap_id_ADC;
		//Calculate weight center
		for(size_t i=0;i<Hit_X.size();i++){
			int layer = id[i]/100000;
			int chip = (id[i]%100000)/10000;
			int chn = id[i]%100;
			wsx_per_layer[layer] += Hit_X[i]*Hit_E[i];
			sx_per_layer[layer] += Hit_E[i];
			wsy_per_layer[layer] += Hit_Y[i]*Hit_E[i];
			sy_per_layer[layer] += Hit_E[i];
			z_per_layer[layer] = Hit_Z[i];
			// umap_id_ADC[id[i]] += Hit_ADC[i]/1000.;
			umap_id_ADC[id[i]] += Hit_E[i];
			// if(layer==15){
			// 	std::cout<<Hit_ADC[i]/1000.<<std::endl;
			// }
		}
		std::array<float,N> wcx_per_layer={};
		std::array<float,N> wcy_per_layer={};
		for(int i=0;i<N;i++){
			wcx_per_layer[i] = wsx_per_layer[i]/sx_per_layer[i];
			wcy_per_layer[i] = wsy_per_layer[i]/sy_per_layer[i];
			float wcx,wcy; //weight center x and y
			if(i%2==0){//even layer 
				wcy = round((wcy_per_layer[i] + 108.65)/5.3)*5.3 - 108.65;
				wcx = round((wcx_per_layer[i] + 90.8)/45.4)*45.4 - 90.8;
			}
			else{//odd layer
				wcx = round((wcx_per_layer[i] + 108.65)/5.3)*5.3 - 108.65;
				wcy = round((wcy_per_layer[i] + 90.8)/45.4)*45.4 - 90.8;
			}
			int i_wcx,i_wcy,i_z;
			i_wcx = int(round(wcx));
			i_wcy = int(round(wcy));
			i_z = int(round(z_per_layer[i]+1e-8));
			if(abs(i_wcx)>1000 || abs(i_wcy)>1000){
				continue;
			}
			//Fill hist if the center falls into some channel
			//continue if it doesn't fall into middle of that cell
			bool fall=false;
			if(i%2==0){ //even layer
				if(abs(wcy_per_layer[i]-wcy)<(1.5) && abs(wcx_per_layer[i]-wcx)<(25.))fall=true;
				//judge if it falls into the center of the layer
				if(abs(wcy_per_layer[i]-9.3)<2.5 && abs(wcx_per_layer[i]-0.)<25.){
					float sum_adc=0.;
					for(auto j:even_chns){
					    int cellid = i*100000 + 20000 + j;
						if(sfs.count(cellid)==0){
							sfs[cellid] = 1.;
						}
					    if(umap_id_ADC.count(cellid)!=0)sum_adc += ( applysf ? umap_id_ADC.at(cellid)*sfs.at(cellid) : umap_id_ADC.at(cellid) );
					}
					if(umap_layer_th1d[i]==nullptr){
						umap_layer_th1d[i] = new TH1D(TString("hsum_")+TString(std::to_string(i)),"test",100,0,230);
					}
					umap_layer_th1d[i]->Fill(sum_adc);
				}
			}
			else{ //odd layer
				if(abs(wcx_per_layer[i]-wcx)<(1.5) && abs(wcy_per_layer[i]-wcy)<(25.))fall=true;
				//judge if it falls into the center of the layer
				if(abs(wcx_per_layer[i]-(-2.))<2.5 && abs(wcy_per_layer[i]-0.)<25.){
					float sum_adc=0.;
					for(auto j:odd_chns){
					    int cellid = i*100000 + 20000 + j;
						if(sfs.count(cellid)==0){
							sfs[cellid] = 1.;
						}
					    if(umap_id_ADC.count(cellid)!=0)sum_adc += ( applysf ? umap_id_ADC.at(cellid)*sfs.at(cellid) : umap_id_ADC.at(cellid) );
					}
					if(umap_layer_th1d[i]==nullptr){
						umap_layer_th1d[i] = new TH1D(TString("hsum_")+TString(std::to_string(i)),"test",100,0,230);
					}
					umap_layer_th1d[i]->Fill(sum_adc);
				}
			}
			if(!fall){
				// std::cout<<"ERROR"<<std::endl;
				continue;
			}
			int cellid = umap_xyz_id[makeKey(i_wcx,i_wcy,i_z)];
			// std::cout<<cellid<<" "<<i_wcx<<" "<<i_wcy<<" "<<i_z<<std::endl;
			if(umap_id_th1d[cellid]==nullptr){
				umap_id_th1d[cellid] = new TH1D(TString("h")+TString(std::to_string(cellid)),"test",50,0,50);
			}
			if(applysf && sfs.count(cellid)!=0)umap_id_th1d[cellid]->Fill(umap_id_ADC[cellid] * sfs.at(cellid));
			else {umap_id_th1d[cellid]->Fill(umap_id_ADC[cellid]);}
		}
	},{"CellID","Hit_X","Hit_Y","Hit_Z","Hit_Energy","Hit_ADC"});

	fout->cd();
	for(int i=0;i<N;i++){
		fout->mkdir(TString("layer")+TString(std::to_string(i)));
	}
	std::map<int,TH1D*> map_id_th1d(umap_id_th1d.begin(), umap_id_th1d.end());
	
	std::unordered_map<int,float> umap_stdmax;
	std::unordered_map<int,TH1D*> umap_layer_std;
	for(auto i:map_id_th1d){
	    int id = i.first;
	    int layer= id/100000;
		int chip = (id%100000)/10000;
		int chn = id%100;
		//if(smooth_layer.count(layer)==1 && applysf==false)i.second->Smooth(1,"G");
		if(chip!=2)continue;
		if(chn<7 || chn > 20)continue;
		if(chn>=15 && chn<=17) continue;
		if(layer%2==0 && chn>=15)continue;
		if(chn == stdcell[layer]){
			umap_stdmax[layer] = getMaximumBinCenter(i.second);
			umap_layer_std[layer] = i.second;
		}
	}

	for(auto i:map_id_th1d){
		_cellid=i.first;
		int layer= _cellid/100000;
		int chip = (_cellid%100000)/10000;
		int chn = _cellid%100;
		if(chip!=2)continue;
		if(chn<7 || chn >= 20)continue;
		if(chn>=15 && chn<=17) continue;
		if(layer%2==0 && chn>=15)continue;
		_sf = 1.;
		if(umap_stdmax[layer]>1.){
			float _sf1 = umap_stdmax[layer]/getMaximumBinCenter(i.second);
			_sf = getSF(i.second,umap_layer_std[layer]);
			// if(layer==11 && chn!=19){_sf=1;}
			std::cout<<_cellid<<" "<<_sf1<<" "<<_sf<<std::endl;

			t->Fill();
		}
	}
	for(auto i:umap_layer_th1d){
		int layer = i.first;
		fout->cd(TString("layer")+TString(std::to_string(layer)));
		i.second->Write();
	}
	for(auto i:map_id_th1d){
		int id = i.first;
		int layer= id/100000;
		int chip = (id%100000)/10000;
		int chn = id%100;
		fout->cd(TString("layer")+TString(std::to_string(layer)));
		if(chip!=2)continue;
		i.second->Write();
	}
	fsf->cd();
	t->Write();
	fout->Close();
}
void draw(){
	//for(int i=5;i<=12;i+=1){
	//	dweight(i,"/mnt2/USTC/jxwang/CEPC_ScECAL/CEPCforZhen/Calib_IncludeADC/40GeV_hl_electron.root","Calib_Hit","Hit_ADC");
	//}
	// shower_profile();
	shift("/mnt2/USTC/jxwang/CEPC_ScECAL/CEPCforZhen/Calib_IncludeADC/40GeV_hl_electron.root","Calib_Hit",false);
	shift("/mnt2/USTC/jxwang/CEPC_ScECAL/CEPCforZhen/Calib_IncludeADC/40GeV_hl_electron.root","Calib_Hit",true);
	// shift("/mnt2/USTC/jxwang/CEPC_ScECAL/CEPCforZhen/MCTruth/electron_40GeV_100000.root","MC_Truth",false);
}

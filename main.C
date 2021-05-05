#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string.h>
#include<algorithm>
#include <sstream>
#include <cstring>

#include <TString.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TTree.h>



using namespace std;

class WaveAnalysis{
public:

    string title1 = ""; // ref
    string title2 = ""; // dut
    string fout_title = "";

    int STEP; //402 or 802
    double WAVE1_POLARITY = -1; // positive = 1, negitive = -1
    double WAVE2_POLARITY = 1; // positive = 1, negitive = -1
    
    string SENSOR_TYPE = "LGAD"; 

    double DELTA_TIME = 0.025; // [ns]

    double THRESHOLD = 15; //[mv]
    //double THRESHOLD = 1; //[mv]
    //double Q_THREHOLD = 0.001; //[fC]
    double Q_THREHOLD = 0.5; //[fC]
    double PEAK_THRESHOLD = 15; //[mV]

    double N = 3000;
    double BASELINE_START = 0; //[ns]
    double BASELINE_STOP = 4; //[ns]

    double AMPLIFIER_TIMS = 100;
    int SIGNAL_FLAG = 1;

    TFile *FOUT;
    TDirectory *WAVED;
    TTree *TREE;

    WaveAnalysis(string t1, string t2);

    ~WaveAnalysis();

    void Charge_distribution(double gate_width);
    void CFD();
    void FullAnalysis();

private:

    int get_steps(string title);
    double *get_peak_parameters(double *waveform_point);
    double get_baseline(double *waveform_point, double win_start, double win_stop);
    double get_charge(TString file,double baseline,int step,int gate_width);
    void get_waveform_data(string file,double *waveform,double wave_polarity);
    double get_CFD(double *waveform,double baseline,double CFD,double peakV,double peakP);

};


WaveAnalysis::WaveAnalysis(string t1,string t2){
    
    title1 = t1;
    title2 = t2;

    string fout_title = title2;
    replace(fout_title.begin(),fout_title.end(),'/','-');
    FOUT = new TFile(Form("./output/WaveAnalysis_%s.root",fout_title.data()),"recreate");
    WAVED = FOUT->mkdir("waveform");


}

WaveAnalysis::~WaveAnalysis(){}


void WaveAnalysis::Charge_distribution(double gate_width){

    SIGNAL_FLAG = 1;

    double *peak;
    double peakV;
    double peakP;
    double charge,rise_time, baseline;

    TREE = new TTree("tree","tree");
    TREE->Branch("charge",&charge,"charge/D");
    TREE->Branch("peak",&peakV,"peak/D");
    TREE->Branch("rise_time",&rise_time,"rise_time/D");
    TREE->Branch("baseline",&baseline,"noise/D");

    STEP = get_steps(title2);

    string file;

    double waveform_point[STEP];



    for(int i =0; i<N; i++){
        cout<<"**********     Entry number : "<<i+1<<"     **********";
        char Istr[5]={0};
        sprintf(Istr,"%05d",i);
    
        file= title2+Istr;
        get_waveform_data(file,waveform_point, WAVE2_POLARITY);
        peak = get_peak_parameters(waveform_point);
        peakV = peak[0];  
        peakP = peak[1];
        
        baseline = get_baseline(waveform_point,BASELINE_START,BASELINE_STOP);
        cout<<"Basekine: "<<baseline<<endl;
        int gateL = peakP-gate_width/DELTA_TIME/2;
        int gateR = peakP+gate_width/DELTA_TIME/2;
        
        double Q = 0;
        int count_point = 0;
        for(int i_temp=peakP;i_temp<gateR+1;i_temp++){
            Q +=waveform_point[i_temp];
            count_point++;
            if(waveform_point[i_temp]<baseline){
                cout<<"gateR:   "<<i_temp*DELTA_TIME<<"  ns"<<endl;
                break;
            }
        }

        for (int j_temp=peakP-1; j_temp>gateL;j_temp--){
            Q +=waveform_point[j_temp];
            count_point++;
            if(waveform_point[j_temp]<baseline){
                cout<<"gateL:   "<<j_temp*DELTA_TIME<<"  ns"<<endl;

                rise_time = peakP*DELTA_TIME-j_temp*DELTA_TIME;

                cout<<"start time:  "<<j_temp*DELTA_TIME<<" ns"<<endl;
                cout<<"peak time:  "<<peakP*DELTA_TIME<<" ns"<<endl;
                cout<<"rise time:  "<<rise_time<<" ns"<<endl;

                break;
            }
        }

        Q-=count_point*baseline;
        Q*=20*DELTA_TIME/AMPLIFIER_TIMS; //   Q = wave * 1/1000 / 50 ohm * 1E-9 *1E15 = 20 * wave

        charge = Q;

        if (peakV>THRESHOLD && Q>Q_THREHOLD){
            
            cout<<"Entry number: "<<i<<endl;
            cout<<"Charge:  "<<Q<<"  fC"<<endl;
            cout<<"Peak:    "<<peakV<< "  mV"<<endl;
            cout<<"Rise Time: "<<rise_time<<"  ns"<<endl;
            cout<<"Noise: "<<baseline<<" mV"<<endl;

            TREE->Fill();

        }
        
    }


    FOUT->cd();

    TREE->Write();

}



void WaveAnalysis::CFD(){

    SIGNAL_FLAG = 0;

    string file1,file2;
    double CFD;
    int FillCounter = 0;

    STEP = get_steps(title2);

    double waveform_point1[STEP],waveform_point2[STEP];
    double *peak1,*peak2,peakV1,peakV2,peakP1,peakP2;
    double baseline1,baseline2;
    double CFDt1,CFDt2,CFDt;

    double xmin = -6; //[ps]
    double xmax = -4; //[ps]
    int nbin = 100;

    TH1F *h1=new TH1F("Cap10_0.1","CFD-0.1",nbin,xmin,xmax);
    TH1F *h2=new TH1F("Cap10_0.2","CFD-0.2",nbin,xmin,xmax);
    TH1F *h3=new TH1F("Cap10_0.3","CFD-0.3",nbin,xmin,xmax);
    TH1F *h4=new TH1F("Cap10_0.4","CFD-0.4",nbin,xmin,xmax);
    TH1F *h5=new TH1F("Cap10_0.5","CFD-0.5",nbin,xmin,xmax);
    TH1F *h6=new TH1F("Cap10_0.6","CFD-0.6",nbin,xmin,xmax);
    TH1F *h7=new TH1F("Cap10_0.7","CFD-0.7",nbin,xmin,xmax);
    TH1F *h8=new TH1F("Cap10_0.8","CFD-0.8",nbin,xmin,xmax);
    TH1F *h9=new TH1F("Cap10_0.9","CFD-0.9",nbin,xmin,xmax);

    for(int i=0; i<N; i++){
        //if (i == 2737){continue;}
        cout<<"**********     Entry number : "<<i+1<<"     **********";

        char Istr[5]={0};
        sprintf(Istr,"%05d",i);
        file1 = title1+Istr;
        file2 = title2+Istr;
        
        get_waveform_data(file1,waveform_point1, WAVE1_POLARITY);
        get_waveform_data(file2,waveform_point2, WAVE2_POLARITY);

        peak1 = get_peak_parameters(waveform_point1);
        peakV1 = peak1[0];
        peakP1 = peak1[1];

        peak2 = get_peak_parameters(waveform_point2);
        peakV2 = peak2[0];
        peakP2 = peak2[1];
        
        baseline1 = get_baseline(waveform_point1,BASELINE_START,BASELINE_STOP);
        baseline2 = get_baseline(waveform_point2,BASELINE_START,BASELINE_STOP);

        if (peakV1>THRESHOLD && peakV1<700 && peakV2>THRESHOLD && peakV2<700){
            FillCounter++;
            for (CFD=0.1; CFD<0.91; CFD += 0.1){
                CFDt1 = get_CFD(waveform_point1,baseline1,CFD,peakV1,peakP1);
                CFDt2 = get_CFD(waveform_point2,baseline2,CFD,peakV2,peakP2);
                CFDt = CFDt1-CFDt2;

                cout<<"CFDt :"<<CFDt<<endl;

                int intCFD = CFD*10+0.1;
                
                switch (intCFD)
                {
                    case 1 :
                        h1->Fill(CFDt);;
                        break;
                    case 2 :
                        h2->Fill(CFDt);;
                        break;
                    case 3 :
                        h3->Fill(CFDt);;
                        break;
                    case 4 :
                        h4->Fill(CFDt);;
                        break;
                    case 5 :
                        h5->Fill(CFDt);;
                        break;
                    case 6 :
                        h6->Fill(CFDt);;
                        break;
                    case 7 :
                        h7->Fill(CFDt);;
                        break;
                    case 8 :
                        h8->Fill(CFDt);;
                        break;
                    case 9 :
                        h9->Fill(CFDt);;
                        break;
                    default:
                        cout << " **** CFD Error! ****"<<endl;
                }
            }
        }
    }

    TCanvas *canvas_time = new TCanvas("canvas_time","Time_Resolution",200,10,2000,500);
    canvas_time->Divide(5,2);
    canvas_time->SetGrid();

    TF1 *g = new TF1("g","gaus",0,5);
    double par[27];

    canvas_time->cd(1);
    //h1->Draw();
    h1->Fit(g);
    g->GetParameters(&par[0]);
    
    h5->SetLineWidth(2);
    h5->SetLineColor(1);
    h5->GetXaxis()->SetTitle("T_{DUT}-T_{Ref} [ns]");
    h5->GetYaxis()->SetTitle("Count");
    h5->GetYaxis()->SetTitleOffset(1.3);


    canvas_time->cd(2); h2->Draw();h2->Fit(g); g->GetParameters(&par[3]);
    canvas_time->cd(3); h3->Draw();h3->Fit(g); g->GetParameters(&par[6]);
    canvas_time->cd(4); h4->Draw();h4->Fit(g); g->GetParameters(&par[9]);
    canvas_time->cd(5); h5->Draw();h5->Fit(g); g->GetParameters(&par[12]);
    canvas_time->cd(6); h6->Draw();h6->Fit(g); g->GetParameters(&par[15]);
    canvas_time->cd(7); h7->Draw();h7->Fit(g); g->GetParameters(&par[18]);
    canvas_time->cd(8); h8->Draw();h8->Fit(g); g->GetParameters(&par[21]);
    canvas_time->cd(9); h9->Draw();h9->Fit(g); g->GetParameters(&par[24]);
    gStyle->SetOptFit(1011);
    
    double CFDFraction[9]={0.1,.2,.3,.4,.5,.6,.7,.8,.9},parF[9]={1000*par[2],1000*par[5],1000*par[8],1000*par[11],1000*par[14],1000*par[17],1000*par[20],1000*par[23],1000*par[26]};
    TGraph *gF=new TGraph(9,CFDFraction,parF);
    canvas_time->cd(10);
    gF->Draw("AC*");
    gF->SetTitle("CFD Fraction ");
    gF->GetXaxis()->SetTitle("CFD Fraction");
    gF->GetYaxis()->SetTitle("Resolution [ps]");

    FOUT->cd();
    canvas_time->Write("time_resolution");

}

void WaveAnalysis::FullAnalysis(){
    Charge_distribution(4);
    CFD();
}


int WaveAnalysis::get_steps(string title){

    string tmp_csvs = title+"00000.txt";

    const char *tmp_csvchar = tmp_csvs.c_str();

    TGraph *tmp_g = new TGraph(tmp_csvchar,"%lf,%lf");

    int steps = tmp_g->GetN();

    return steps; 

}



void WaveAnalysis::get_waveform_data(string file,double *waveform, double wave_polarity){

    cout<<"\nReading :  "<<file+".txt"<<endl;
    
    string tmp_csvs = file+".txt";
    const char *tmp_csvchar = tmp_csvs.c_str();

    TGraph *tmp_g = new TGraph(tmp_csvchar,"%lf,%lf");

    int num_point = tmp_g->GetN();
    STEP = num_point;
    cout<<"\nNum Points: "<<num_point<<endl;

    double *tmp_waveformX = tmp_g->GetX();
    double *tmp_waveformY = tmp_g->GetY();

    double tmp_time[num_point];


    for (int i=0; i<num_point; i++){
        tmp_time[i] = i*DELTA_TIME;
        waveform[i] = tmp_waveformY[i]*(1000)*wave_polarity;
    }

    if(SIGNAL_FLAG == 1){
        TGraph *wave_g = new TGraph(num_point,tmp_time,waveform);
        wave_g->SetNameTitle(tmp_csvchar,tmp_csvchar);
        wave_g->GetXaxis()->SetTitle("Time [ns]");
        wave_g->GetYaxis()->SetTitle("Amplitude [mV]");
        WAVED->cd();
        wave_g->Write(tmp_csvchar);
    }


}



//void WaveAnalysis::get_waveform_data(string file,double *waveform){
//    
//    cout<<"\nReading :  "<<file<<endl;
//
//    ifstream in;
//    cout<<file<<endl;
//    in.open(file+".csv");
//    int i;
//    double t[STEP];
//
//    //line 1 to line 5 are comments
//    char tmp[255];
//    for( int lineNo=0; lineNo<5; lineNo++)
//    {
//      in.getline(tmp, 255);
//      //cout<<"tmp  -----  "<<tmp<<endl;
//    }// */
//
//    string line;
//    for(i=0;i<415;i++){
//        //cout<<"########## DEBUG i "<<i<<endl;
//        getline(in,line);//读取每行数据
//        string number;
//        istringstream readstr(line); //string数据流化
//        //将一行数据按'，'分割
//        for(int j = 0;j < 2;j++){ //可根据数据的实际情况取循环获取
//            getline(readstr,number,','); //循环读取数据
//            if (j == 1) waveform[i] = (atof(number.c_str())); // str to float
//           // cout<<"########## DEBUG j "<<j<<endl;
//        }
//
//       // cout<<"########## DEBUG"<<endl;
//        waveform[i] *= -1000;//
//        //cout<<"-------------waveform  "<<i<<"   "<<waveform[i]<<endl;
//        t[i]=i*DELTA_TIME;
//    }
//
//       // TGraph *g=new TGraph(STEP,t,waveform);
//       // TCanvas *c = new TCanvas("Waveform","Waveform",200,10,1000,500);
//       // c->SetGrid();
//       // g->Draw("AC*");
//       // g->SetTitle("Waveform of ");
//       // g->GetXaxis()->SetTitle("t/ns");
//       // g->GetYaxis()->SetTitle("V/mV");
//        //}
//
//        in.close();    
//}

double *WaveAnalysis::get_peak_parameters(double *waveform_point){

    int i = 0;
    double temp;
    static double peak[2];
    peak[0] = 0;
    peak[1] = 0;

    for(peak[0] = 0,i=0;i<STEP;i++){
        temp = waveform_point[i];
        //cout<<min<<endl;
        if(temp > peak[0]){
           peak[0] = temp;
           peak[1]=i;
        }
    }

    return(peak);
    
}

double WaveAnalysis::get_baseline(double *waveform_point, double win_start, double win_stop){
    
    int n_start,n_stop,width=0;
    n_start = (int)(win_start/DELTA_TIME);
    n_stop = (int)(win_stop/DELTA_TIME);

    double sum=0.0,baseline=0.0;

    for(int i=n_start; i<n_stop; i++){
        sum += waveform_point[i];
        width++;
    }

    baseline = sum/width;
    
    return(baseline);

}

double WaveAnalysis::get_CFD(double *waveform,double baseline,double CFD,double peakV,double peakP){

    double CFDt;
    double v_cfd = CFD*(peakV-baseline)+baseline;
    int i_cfd=0, fit_point_num = 6;
    double fitT1,fitT2,fitStatus = 0;

    vector <double> value;
    vector <double> t;

    for (int i=peakP; i>0; i--){
        if (waveform[i] < v_cfd){ i_cfd = i; break; }
    }

    for(int j = i_cfd-fit_point_num/2;j<i_cfd+fit_point_num/2+1;j++){
        value.push_back(waveform[j]);
        t.push_back(j*DELTA_TIME);
    }

    fitT1 = (i_cfd-fit_point_num/2)*DELTA_TIME;
    fitT2 = (i_cfd+fit_point_num/2+1)*DELTA_TIME;
    cout<<"Fit range[ns]:  "<< fitT1<<" - "<<fitT2<<endl;

    TGraph *g = new TGraph(fit_point_num,&t[0],&value[0]);
    TF1 *f = new TF1("f","pol3",fitT1,fitT2);
    g->Fit(f,"R+");
    if(f->GetChisquare()>1e-4) fitStatus=1;
    CFDt = f->GetX(v_cfd, fitT1, fitT2, 1e-14);

    delete g;
    delete f;
    return(CFDt);
}

int main()
{

    gStyle->SetOptStat("e");
    //gStyle->SetOptStat(0000);
    gStyle->SetOptFit(111);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.18);
    gStyle->SetStatH(0.22);

    // No Canvas Border                  
    gStyle->SetCanvasBorderMode(0);      
    gStyle->SetCanvasBorderSize(0);      
    // White BG                          
    gStyle->SetCanvasColor(10);          
    // Format for axes                   
    gStyle->SetLabelFont(42,"xyz");      
    gStyle->SetLabelSize(0.04,"xyz");    
    gStyle->SetLabelOffset(0.01,"xyz");  
    gStyle->SetNdivisions(510,"xyz");    
    gStyle->SetTitleFont(42,"xyz");      
    gStyle->SetTitleColor(1,"xyz");      
    gStyle->SetTitleSize(0.05,"xyz");    
    gStyle->SetTitleOffset(1.5,"xyz");  
    // No pad borders                    
    gStyle->SetPadBorderMode(0);         
    gStyle->SetPadBorderSize(0);         
    // White BG                          
    gStyle->SetPadColor(10);             
    // Margins for labels etc.           
    gStyle->SetPadLeftMargin(0.15);      
    gStyle->SetPadBottomMargin(0.15);    
    //gStyle->SetPadRightMargin(0.05);     
    //gStyle->SetPadTopMargin(0.05);       
    // No error bars in x direction      
    gStyle->SetErrorX(0);                
    // Format legend                     
    gStyle->SetLegendBorderSize(0); 




//W7-II-L1-15-70
    
    //string W7_350V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_350V/C2BiasV";         
    //string w7_350V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_350V/C3BiasV";         
    //WaveAnalysis W7_350V(W7_350V_t1,w7_350V_t2);                                                                                                                                                                                                 
    //W7_350V.FullAnalysis();  

    //string W7_330V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_330V/C2BiasV";         
    //string w7_330V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_330V/C3BiasV";         
    //WaveAnalysis W7_330V(W7_330V_t1,w7_330V_t2);                                                                                                                                                                                                 
    //W7_330V.FullAnalysis();  


    //string W7_310V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_310V/C2BiasV";         
    //string w7_310V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_310V/C3BiasV";         
    //WaveAnalysis W7_310V(W7_310V_t1,w7_310V_t2);                                                                                                                                                                                                 
    //W7_310V.FullAnalysis();  

    //string W7_290V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_290V/C2BiasV";         
    //string w7_290V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_290V/C3BiasV";         
    //WaveAnalysis W7_290V(W7_290V_t1,w7_290V_t2);                                                                                                                                                                                                 
    //W7_290V.FullAnalysis();  

    //string W7_270V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_270V/C2BiasV";         
    //string w7_270V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_270V/C3BiasV";         
    //WaveAnalysis W7_270V(W7_270V_t1,w7_270V_t2);                                                                                                                                                                                                 
    //W7_270V.FullAnalysis();  


    //string W7_250V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_250V/C2BiasV";         
    //string w7_250V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_250V/C3BiasV";         
    //WaveAnalysis W7_250V(W7_250V_t1,w7_250V_t2);                                                                                                                                                                                                 
    //W7_250V.FullAnalysis();  

    //string W7_230V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_230V/C2BiasV";         
    //string w7_230V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_230V/C3BiasV";         
    //WaveAnalysis W7_230V(W7_230V_t1,w7_230V_t2);                                                                                                                                                                                                 
    //W7_230V.FullAnalysis();  

    //string W7_210V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_210V/C2BiasV";         
    //string w7_210V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_210V/C3BiasV";         
    //WaveAnalysis W7_210V(W7_210V_t1,w7_210V_t2);                                                                                                                                                                                                 
    //W7_210V.FullAnalysis();  

    //string W7_200V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_200V/C2BiasV";         
    //string w7_200V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_200V/C3BiasV";         
    //WaveAnalysis W7_200V(W7_200V_t1,w7_200V_t2);                                                                                                                                                                                                 
    //W7_200V.FullAnalysis();  

    //string W7_190V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_190V/C2BiasV";         
    //string w7_190V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_190V/C3BiasV";         
    //WaveAnalysis W7_190V(W7_190V_t1,w7_190V_t2);                                                                                                                                                                                                 
    //W7_190V.FullAnalysis();  

    //string W7_170V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_170V/C2BiasV";         
    //string w7_170V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_170V/C3BiasV";         
    //WaveAnalysis W7_170V(W7_170V_t1,w7_170V_t2);                                                                                                                                                                                                 
    //W7_170V.FullAnalysis();  

    //string W7_150V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_150V/C2BiasV";         
    //string w7_150V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_150V/C3BiasV";         
    //WaveAnalysis W7_150V(W7_150V_t1,w7_150V_t2);                                                                                                                                                                                                 
    //W7_150V.FullAnalysis();  


//W7-III-L1-15-70

   // string W7_90V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L1-15_70_beta/220V_90V/C2BiasV";         
   // string w7_90V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L1-15_70_beta/220V_90V/C3BiasV";         
   // WaveAnalysis W7_90V(W7_90V_t1,w7_90V_t2);                                                                                                                                                                                                 
   // W7_90V.FullAnalysis();  

   // string W7_80V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L1-15_70_beta/220V_80V/C2BiasV";         
   // string w7_80V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L1-15_70_beta/220V_80V/C3BiasV";         
   // WaveAnalysis W7_80V(W7_80V_t1,w7_80V_t2);                                                                                                                                                                                                 
   // W7_80V.FullAnalysis();  

   // string W7_70V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L1-15_70_beta/220V_70V/C2BiasV";         
   // string w7_70V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L1-15_70_beta/220V_70V/C3BiasV";         
   // WaveAnalysis W7_70V(W7_70V_t1,w7_70V_t2);                                                                                                                                                                                                 
   // W7_70V.FullAnalysis();  


//W8-IV-L1-15-70

  //string w8_180V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15-70-beta/200V_180V/C2BiasV";         
  //string w8_180V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15-70-beta/200V_180V/C3BiasV";         
  //WaveAnalysis W8_180V(w8_180V_t1,w8_180V_t2);                                                                                                                                                                                                 
  //W8_180V.FullAnalysis();  


  // string w8_170V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_170V/C2BiasV";         
  // string w8_170V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_170V/C3BiasV";         
  // WaveAnalysis W8_170V(w8_170V_t1,w8_170V_t2);                                                                                                                                                                                                 
  // W8_170V.FullAnalysis();  


  // string w8_160V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_160V/C2BiasV";         
  // string w8_160V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_160V/C3BiasV";         
  // WaveAnalysis W8_160V(w8_160V_t1,w8_160V_t2);                                                                                                                                                                                                 
  // W8_160V.FullAnalysis();  


  // string w8_150V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_150V/C2BiasV";         
  // string w8_150V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_150V/C3BiasV";         
  // WaveAnalysis W8_150V(w8_150V_t1,w8_150V_t2);                                                                                                                                                                                                 
  // W8_150V.FullAnalysis();  


  // string w8_140V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_140V/C2BiasV";         
  // string w8_140V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_140V/C3BiasV";         
  // WaveAnalysis W8_140V(w8_140V_t1,w8_140V_t2);                                                                                                                                                                                                 
  // W8_140V.FullAnalysis();  


  // string w8_130V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_130V/C2BiasV";         
  // string w8_130V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_130V/C3BiasV";         
  // WaveAnalysis W8_130V(w8_130V_t1,w8_130V_t2);                                                                                                                                                                                                 
  // W8_130V.FullAnalysis();  


  // string w8_120V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_120V/C2BiasV";         
  // string w8_120V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_120V/C3BiasV";         
  // WaveAnalysis W8_120V(w8_120V_t1,w8_120V_t2);                                                                                                                                                                                                 
  // W8_120V.FullAnalysis();  


  // string w8_110V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_110V/C2BiasV";         
  // string w8_110V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_110V/C3BiasV";         
  // WaveAnalysis W8_110V(w8_110V_t1,w8_110V_t2);                                                                                                                                                                                                 
  // W8_110V.FullAnalysis();  


  //string w8_100V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_100V/C2BiasV";         
  //string w8_100V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_100V/C3BiasV";         
  //WaveAnalysis W8_100V(w8_100V_t1,w8_100V_t2);                                                                                                                                                                                                 
  //W8_100V.FullAnalysis();  


  // string w8_90V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_90V/C2BiasV";         
  // string w8_90V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_90V/C3BiasV";         
  // WaveAnalysis W8_90V(w8_90V_t1,w8_90V_t2);                                                                                                                                                                                                 
  // W8_90V.FullAnalysis();  


  // string w8_80V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_80V/C2BiasV";         
  // string w8_80V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_80V/C3BiasV";         
  // WaveAnalysis W8_80V(w8_80V_t1,w8_80V_t2);                                                                                                                                                                                                 
  // W8_80V.FullAnalysis();  


  // string w8_70V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_70V/C2BiasV";         
  // string w8_70V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_70V/C3BiasV";         
  // WaveAnalysis W8_70V(w8_70V_t1,w8_70V_t2);                                                                                                                                                                                                 
  // W8_70V.FullAnalysis();  


  // string w8_60V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_60V/C2BiasV";         
  // string w8_60V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_60V/C3BiasV";         
  // WaveAnalysis W8_60V(w8_60V_t1,w8_60V_t2);                                                                                                                                                                                                 
  // W8_60V.FullAnalysis();  




//---------------------------------------------------------------------------------------------------------------
//Large Array
//---------------------------------------------------------------------------------------------------------------


//W7-III-L4-15-100 

 // string w7_140V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L4-15_100_beta/220V_140V/C2BiasV";         
 // string w7_140V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L4-15_100_beta/220V_140V/C3BiasV";         
 // WaveAnalysis W7_140V(w7_140V_t1,w7_140V_t2);                                                                                                                                                                                                 
 // W7_140V.FullAnalysis();  


 // string w7_120V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L4-15_100_beta/220V_120V/C2BiasV";         
 // string w7_120V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L4-15_100_beta/220V_120V/C3BiasV";         
 // WaveAnalysis W7_120V(w7_120V_t1,w7_120V_t2);                                                                                                                                                                                                 
 // W7_120V.FullAnalysis();  

 // string w7_100V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L4-15_100_beta/220V_100V/C2BiasV";         
 // string w7_100V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W7-III-F4-L4-15_100_beta/220V_100V/C3BiasV";         
 // WaveAnalysis W7_100V(w7_100V_t1,w7_100V_t2);                                                                                                                                                                                                 
 // W7_100V.FullAnalysis();  

//W8-IV-L4-15-100

    //string w8_180V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_180V/C2BiasV";         
    //string w8_180V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_180V/C3BiasV";         
    //WaveAnalysis W8_180V(w8_180V_t1,w8_180V_t2);                                                                                                                                                                                                 
    //W8_180V.FullAnalysis();  

    //string w8_160V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_160V/C2BiasV";         
    //string w8_160V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_160V/C3BiasV";         
    //WaveAnalysis W8_160V(w8_160V_t1,w8_160V_t2);                                                                                                                                                                                                 
    //W8_160V.FullAnalysis();  

    //string w8_140V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_140V/C2BiasV";         
    //string w8_140V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_140V/C3BiasV";         
    //WaveAnalysis W8_140V(w8_140V_t1,w8_140V_t2);                                                                                                                                                                                                 
    //W8_140V.FullAnalysis();  

    //string w8_120V_t1 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_120V/C2BiasV";         
    //string w8_120V_t2 = "/home/admin/STDB/BetaTest/IMEV1beta/B58_1_A3_W8-IV-E4-L4-15_100_beta/220V_120V/C3BiasV";         
    //WaveAnalysis W8_120V(w8_120V_t1,w8_120V_t2);                                                                                                                                                                                                 
    //W8_120V.FullAnalysis();  


//-----------------------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------------------

    //string PIN_400V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/400V/C3BiasV";
    //string PIN_400V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/400V/C3BiasV";
    //WaveAnalysis PIN_400V(PIN_400V_t1,PIN_400V_t2);
    //PIN_400V.Charge_distribution(4);

    //string PIN_350V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/350V/C3BiasV";
    //string PIN_350V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/350V/C3BiasV";
    //WaveAnalysis PIN_350V(PIN_350V_t1,PIN_350V_t2);
    //PIN_350V.Charge_distribution(4);

    //string PIN_300V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/300V/C3BiasV";
    //string PIN_300V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/300V/C3BiasV";
    //WaveAnalysis PIN_300V(PIN_300V_t1,PIN_300V_t2);
    //PIN_300V.Charge_distribution(4);

    //string PIN_250V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/250V/C3BiasV";
    //string PIN_250V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/250V/C3BiasV";
    //WaveAnalysis PIN_250V(PIN_250V_t1,PIN_250V_t2);
    //PIN_250V.Charge_distribution(4);

    //string PIN_200V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/200V/C3BiasV";
    //string PIN_200V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/200V/C3BiasV";
    //WaveAnalysis PIN_200V(PIN_200V_t1,PIN_200V_t2);
    //PIN_200V.Charge_distribution(4);

    //string PIN_150V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/150V/C3BiasV";
    //string PIN_150V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/150V/C3BiasV";
    //WaveAnalysis PIN_150V(PIN_150V_t1,PIN_150V_t2);
    //PIN_150V.Charge_distribution(4);

    //string PIN_100V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/100V/C3BiasV";
    //string PIN_100V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-100-beta/100V/C3BiasV";
    //WaveAnalysis PIN_100V(PIN_100V_t1,PIN_100V_t2);
    //PIN_100V.Charge_distribution(4);

    //string PIN_50V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-50-beta/50V/C3BiasV";
    //string PIN_50V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-PIN-15-50-beta/50V/C3BiasV";
    //WaveAnalysis PIN_50V(PIN_50V_t1,PIN_50V_t2);
    //PIN_50V.Charge_distribution(4);

    //string w8_180V_t1 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-L1-15-70-beta/w8-IV-E4-L1-15_70_beta_180V/C3BiasV";
    //string w8_180V_t2 = "/home/admin/STDB/BetaTest/data/IMEV1beta/One_Chanel/w8-IV-E4-L1-15-70-beta/w8-IV-E4-L1-15_70_beta_180V/C3BiasV";
    //WaveAnalysis W8_180V(w8_180V_t1,w8_180V_t2);
    //W8_180V.Charge_distribution(4);



//-----------------------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------------------

    string SiC_PIN_500V_t1 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/20210215/NJU_SiC_PIN_500V/C2--_2.5GT_Ref_Clock_--";
    string SiC_PIN_500V_t2 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/20210215/NJU_SiC_PIN_500V/C3--_2.5GT_Ref_Clock_--";
    WaveAnalysis SiC_PIN_500V(SiC_PIN_500V_t1,SiC_PIN_500V_t2);
    SiC_PIN_500V.FullAnalysis();
    

    //string SiC_Noise_500V_t1 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/2021/Noise/C2--_2.5GT_Ref_Clock_--";
    //string SiC_Noise_500V_t2 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/2021/Noise/C3--_2.5GT_Ref_Clock_--";
    //WaveAnalysis SiC_Noise_500V(SiC_Noise_500V_t1,SiC_Noise_500V_t2);
    //SiC_Noise_500V.Charge_distribution(4);

    //string NDL_200V_t1 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/2021/NJU_SiC_PIN_500V/C2--_2.5GT_Ref_Clock_--";
    //string NDL_200V_t2 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/2021/NJU_SiC_PIN_500V/C2--_2.5GT_Ref_Clock_--";
    //WaveAnalysis NDL_200V(NDL_200V_t1,NDL_200V_t2);
    //NDL_200V.Charge_distribution(4);

    //string NDL_Noise_200V_t1 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/2021/Noise/C2--_2.5GT_Ref_Clock_--";
    //string NDL_Noise_200V_t2 = "/home/admin/STDB/BetaTest/data/4H-SiC-beta/2021/Noise/C2--_2.5GT_Ref_Clock_--";
    //WaveAnalysis NDL_Noise_200V(NDL_Noise_200V_t1,NDL_Noise_200V_t2);
    //NDL_Noise_200V.Charge_distribution(4);



}


#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TLine.h>
#include <TPaveText.h>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooPlot.h>
#include <RooLandau.h>
#include <RooGaussian.h>
#include <RooFFTConvPdf.h>

using namespace RooFit;
using namespace std;



void plot_ana(){


    //gStyle->SetOptStat("e");
    gStyle->SetOptStat(1111);
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



    // LGAD
    //double charge_xmin = 0.;
    //double charge_xmax = 80.;
    //double charge_bin = 40;

    //double peak_xmin = 0.;
    //double peak_xmax = 200.;
    //double peak_bin = 50;
    //
    //double rise_edge_xmin = 0.;
    //double rise_edge_xmax = 500.;
    //double rise_edge_bin = 50;
    //
    //double rise_time_xmin = 0.;
    //double rise_time_xmax = 3.;
    //double rise_time_bin = 60;

    //double baseline_xmin = -7.;
    //double baseline_xmax = 7.;
    //double baseline_bin = 70;

    // PIN
    double charge_xmin = 0.;
    double charge_xmax = 10.;
    double charge_bin = 40;

    double peak_xmin = 0.;
    double peak_xmax = 60.;
    double peak_bin = 60;
    
    double rise_edge_xmin = 0.;
    double rise_edge_xmax = 150.;
    double rise_edge_bin = 50;

    double rise_time_xmin = 0.; 
    double rise_time_xmax = 3.; 
    double rise_time_bin = 60;  
    
    double baseline_xmin = -7.;
    double baseline_xmax = 7.;
    double baseline_bin = 70;

    TH1F *h_charge = new TH1F("","",charge_bin,charge_xmin,charge_xmax);
    TH1F *h_baseline = new TH1F("","",baseline_bin,baseline_xmin,baseline_xmax);
    TH1F *h_rise_time = new TH1F("","",rise_time_bin,rise_time_xmin,rise_time_xmax);
    TH1F *h_peak = new TH1F("","",peak_bin,peak_xmin,peak_xmax);
    TH1F *h_rise_edge = new TH1F("","",rise_edge_bin,rise_edge_xmin,rise_edge_xmax);

    TH2F *h2_charge_peak = new TH2F("","",charge_bin,charge_xmin,charge_xmax,peak_bin,peak_xmin,peak_xmax);
    TH2F *h2_charge_rise_time = new TH2F("","",charge_bin,charge_xmin,charge_xmax,rise_time_bin,rise_time_xmin,rise_time_xmax);

    h2_charge_peak->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_charge_peak->GetYaxis()->SetLimits(peak_xmin,peak_xmax);

    h2_charge_rise_time->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_charge_rise_time->GetYaxis()->SetLimits(rise_time_xmin,rise_time_xmax);


    h_charge->SetLineColor(kRed);
    h_charge->SetLineStyle(kDashed);
    h_charge->SetMarkerColor(kRed);
    h_charge->SetLineWidth(2);

    h_baseline->SetLineColor(kBlue);
    h_baseline->SetLineStyle(kDashed);
    h_baseline->SetMarkerColor(kBlue);
    h_baseline->SetLineWidth(2);
    
    //h_rise_time->SetLineColor(kMagenta);
    //h_rise_time->SetLineStyle(kDashed);
    h_rise_time->SetMarkerColor(kMagenta);
    h_rise_time->SetLineWidth(2);
    
    h_peak->SetLineWidth(2);
    h_peak->SetLineColor(kBlack);

    //h_rise_edge->SetLineColor(417);
    //h_rise_edge->SetLineStyle(kDashed);
    h_rise_edge->SetMarkerColor(417);
    h_rise_edge->SetLineWidth(2);

    h2_charge_rise_time->SetMarkerColor(kOrange);
    h2_charge_peak->SetMarkerColor(kCyan);


    h_charge->GetXaxis()->SetTitle("Collected Charges [fC]");
    h_charge->GetXaxis()->CenterTitle();
    h_charge->GetXaxis()->SetTitleOffset(1.5);
    h_charge->GetYaxis()->SetTitle("Entries");
    h_charge->GetYaxis()->CenterTitle();
    h_charge->GetYaxis()->SetTitleOffset(1.5);

    h_baseline->GetXaxis()->SetTitle("Baseline [mv]");
    h_baseline->GetXaxis()->CenterTitle();
    h_baseline->GetXaxis()->SetTitleOffset(1.5);
    h_baseline->GetYaxis()->SetTitle("Entries");
    h_baseline->GetYaxis()->CenterTitle();
    h_baseline->GetYaxis()->SetTitleOffset(1.5);

    h_rise_time->GetXaxis()->SetTitle("t_{Rise} [ns]");
    h_rise_time->GetXaxis()->CenterTitle();
    h_rise_time->GetXaxis()->SetTitleOffset(1.5);
    h_rise_time->GetYaxis()->SetTitle("Entries");
    h_rise_time->GetYaxis()->CenterTitle();
    h_rise_time->GetYaxis()->SetTitleOffset(1.5);

    h_peak->GetXaxis()->SetTitle("Ampl peak [mV]");
    h_peak->GetXaxis()->CenterTitle();
    h_peak->GetXaxis()->SetTitleOffset(1.5);
    h_peak->GetYaxis()->SetTitle("Entries");
    h_peak->GetYaxis()->CenterTitle();
    h_peak->GetYaxis()->SetTitleOffset(1.5);
    
    h_rise_edge->GetXaxis()->SetTitle("dV/dt [mv/ns]");
    h_rise_edge->GetXaxis()->CenterTitle();
    h_rise_edge->GetXaxis()->SetTitleOffset(1.5);
    h_rise_edge->GetYaxis()->SetTitle("Entries");
    h_rise_edge->GetYaxis()->CenterTitle();
    h_rise_edge->GetYaxis()->SetTitleOffset(1.5);
    
    h2_charge_rise_time->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_charge_rise_time->GetXaxis()->CenterTitle();
    h2_charge_rise_time->GetXaxis()->SetTitleOffset(1.5);
    h2_charge_rise_time->GetYaxis()->SetTitle("Rise Time [ns]");
    h2_charge_rise_time->GetYaxis()->CenterTitle();
    h2_charge_rise_time->GetYaxis()->SetTitleOffset(1.5);

    h2_charge_peak->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_charge_peak->GetXaxis()->CenterTitle();
    h2_charge_peak->GetXaxis()->SetTitleOffset(1.5);
    h2_charge_peak->GetYaxis()->SetTitle("Peak Amplitude [mV]");
    h2_charge_peak->GetYaxis()->CenterTitle();
    h2_charge_peak->GetYaxis()->SetTitleOffset(1.5);


    //TFile *f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-IMEV1beta-One_Chanel-w8-IV-E4-L1-15-70-beta-w8-IV-E4-L1-15_70_beta_180V-C3BiasV.root","READ");
    //TFile *f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-IMEV1beta-One_Chanel-w8-IV-E4-PIN-15-100-beta-250V-C3BiasV.root","READ");
    //TFile *f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-4H-SiC-beta-4H_SiC_PIN_1-500V-C3BiasV.root","READ");
    
    TFile *f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-4H-SiC-beta-2021-NJU_SiC_PIN_500V-C3--_2.5GT_Ref_Clock_--.root","READ");
    //TFile *f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-4H-SiC-beta-2021-Noise-C3--_2.5GT_Ref_Clock_--.root","READ");
    
    //TFile *f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-4H-SiC-beta-2021-NJU_SiC_PIN_500V-C2--_2.5GT_Ref_Clock_--.root","READ");
    //TFile *f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-4H-SiC-beta-2021-Noise-C2--_2.5GT_Ref_Clock_--.root","READ");

    TTree *tree = (TTree *)f_in->Get("tree");

    double charge, peak, rise_time, baseline;

    tree->SetBranchAddress("charge",&charge);
    tree->SetBranchAddress("peak",&peak);
    tree->SetBranchAddress("rise_time",&rise_time);
    tree->SetBranchAddress("baseline",&baseline);

    for(int i=0; i<tree->GetEntries(); i++){

        tree->GetEntry(i);

        h_charge->Fill(charge);
        h_baseline->Fill(baseline);
        h_rise_time->Fill(rise_time);
        h_peak->Fill(peak);
        h_rise_edge->Fill(peak/rise_time);

        h2_charge_peak->Fill(charge,peak);
        h2_charge_rise_time->Fill(charge,rise_time);


    }


    // get mpv for rise_time
    double mpv_rise_time = h_rise_time->GetMaximumBin()*(rise_time_xmax-rise_time_xmin)/rise_time_bin; //[ns]
    TPaveText *rise_time_pt = new TPaveText(1.5,0.7*h_rise_time->GetMaximum(),2.5,0.9*h_rise_time->GetMaximum());
    rise_time_pt->AddText(Form("Rising Time MPV: %.2f [ns]", mpv_rise_time));

    // get mpv for rise_edge
    double mpv_rise_edge = h_rise_edge->GetMaximumBin()*(rise_edge_xmax-rise_edge_xmin)/rise_edge_bin;
    TPaveText *rise_edge_pt = new TPaveText(300,0.7*h_rise_edge->GetMaximum(),450,0.9*h_rise_edge->GetMaximum());
    rise_edge_pt->AddText(Form("Rising Edge MPV: %.2f [mv/ns]", mpv_rise_edge));





//----------------------------------------------------------------------------------------------------------------------


    // fit charge by RooFit
    RooRealVar v_charge("v_charge", "v_charge", charge_xmin, charge_xmax);
    
    RooDataHist data_charge("data_charge","data_charge",v_charge,h_charge);
 
    RooRealVar ml_charge("mean landau", "mean landau", 3, 0, 6);
    RooRealVar sl_charge("sigma landau", "sigma landau", 1., 0.1, 10);
    RooLandau landau_charge("lx", "lx", v_charge, ml_charge, sl_charge);
 
    RooRealVar mg_charge("mean cgauss", "mean gauss", 0);
    RooRealVar sg_charge("sigma cgauss", "sigma gauss", 5., 0.1, 10);
    RooGaussian gauss_charge("gauss_charge", "gauss_charge", v_charge, mg_charge, sg_charge);
 
    RooFFTConvPdf lxg("lxg", "landau (X) gauss", v_charge, landau_charge, gauss_charge);
 
    lxg.fitTo(data_charge);
 
    RooPlot *frame_charge = v_charge.frame(Title("landau (x) gauss convolution"));
    data_charge.plotOn(frame_charge);
    lxg.plotOn(frame_charge);
    //landau_charge.plotOn(frame_charge, LineStyle(kDashed));

    landau_charge.paramOn(frame_charge);

    frame_charge->GetXaxis()->SetTitle("Collected Charges [fC]");
    frame_charge->GetXaxis()->CenterTitle();
    frame_charge->GetXaxis()->SetTitleOffset(1.5);
    frame_charge->GetYaxis()->SetTitle("Entries");
    frame_charge->GetYaxis()->CenterTitle();
    frame_charge->GetYaxis()->SetTitleOffset(1.5);









    //fit baseline
    RooRealVar v_baseline("v_baseline", "v_baseline", baseline_xmin, baseline_xmax);

    RooDataHist data_baseline("data_baselien","data_baseline",v_baseline,h_baseline);

    RooRealVar mg_baseline("mean gauss", "mean gauss", -4.,-7.,7.);
    RooRealVar sg_baseline("sigma gauss", "sigma gauss", 5., 0.1, 10);
    RooGaussian gauss_baseline("gauss_baseline", "gauss_baseline", v_baseline, mg_baseline, sg_baseline);

    gauss_baseline.fitTo(data_baseline);

    RooPlot *frame_baseline = v_baseline.frame();
    data_baseline.plotOn(frame_baseline);
    gauss_baseline.plotOn(frame_baseline);

    frame_baseline->GetXaxis()->SetTitle("Baseline [mv]");
    frame_baseline->GetXaxis()->CenterTitle();
    frame_baseline->GetXaxis()->SetTitleOffset(1.5);
    frame_baseline->GetYaxis()->SetTitle("Entries");
    frame_baseline->GetYaxis()->CenterTitle();
    frame_baseline->GetYaxis()->SetTitleOffset(1.5);

    gauss_baseline.paramOn(frame_baseline);
    







    TCanvas *c = new TCanvas("c","c",2000,1000);
    c->Divide(4,2);
    c->cd(1);
    c->GetPad(1)->SetGrid();
    //h_charge->Draw();
    frame_charge->Draw();

    c->cd(2);
    c->GetPad(2)->SetGrid();
    //h_baseline->Draw();
    frame_baseline->Draw();

    c->cd(3);
    c->GetPad(3)->SetGrid();
    h_rise_time->Draw();
    //rise_time_pt->Draw();

    c->cd(4);
    c->GetPad(4)->SetGrid();
    h_peak->Draw();
    
    c->cd(5);
    c->GetPad(5)->SetGrid();
    h_rise_edge->Draw();
    //rise_edge_pt->Draw();
    
    c->cd(6);
    c->GetPad(6)->SetGrid();
    h2_charge_rise_time->Draw();

    c->cd(7);
    c->GetPad(7)->SetGrid();
    h2_charge_peak->Draw();












    //TCanvas *c1 = new TCanvas("c1","c1",500,500);
    //c1->SetGrid();
    //c1->cd();
    //h_charge->Draw();
    //charge_landau->Draw("same");

    //TCanvas *c2 = new TCanvas("c2","c2",500,500);
    //c2->SetGrid();
    //c2->cd();
    //h_baseline->Draw();
    //baseline_gaus->Draw("same");

    //TCanvas *c3 = new TCanvas("c3","c3",500,500);
    //c3->SetGrid();
    //c3->cd();
    //h_rise_time->Draw();
    //c3->SaveAs("./output/NDL_Trise.png");

    //TCanvas *c4 = new TCanvas("c4","c4",500,500);
    //c4->SetGrid();
    //c4->cd();
    //h_rise_edge->Draw();
    //c4->SaveAs("./output/NDL_dVdt.png");

    //
    //TCanvas *c5 = new TCanvas("c5","c5",500,500);
    //c5->SetGrid();
    //c5->cd();
    //h2_charge_rise_time->Draw();


    //TCanvas *c6 = new TCanvas("c6","c6",500,500);
    //c6->SetGrid();
    //c6->cd();
    //h2_charge_peak->Draw();


}

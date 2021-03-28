#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TLine.h>
#include <TBox.h>

void plot_ana_compare(){


    gStyle->SetOptStat(0000);

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
    double charge_xmin = 0.;
    double charge_xmax = 80.;
    double charge_bin = 40;

    double peak_xmin = 0.;
    double peak_xmax = 700.;
    double peak_bin = 100;
    
    double rise_time_xmin = 0.;
    double rise_time_xmax = 3.;
    double rise_time_bin = 30;

    double baseline_xmin = -7.;
    double baseline_xmax = 7.;
    double baseline_bin = 70;


    TH1F *h_lgad_charge = new TH1F("","",charge_bin,charge_xmin,charge_xmax);
    TH1F *h_lgad_baseline = new TH1F("","",baseline_bin,baseline_xmin,baseline_xmax);
    TH1F *h_lgad_rise_edge = new TH1F("","",50,0,500);
    TH2F *h2_lgad_charge_peak = new TH2F("","",charge_bin,charge_xmin,charge_xmax,peak_bin,peak_xmin,peak_xmax);
    TH2F *h2_lgad_charge_rise_time = new TH2F("","",charge_bin,charge_xmin,charge_xmax,rise_time_bin,rise_time_xmin,rise_time_xmax);
    TH2F *h2_lgad_rise_time_peak = new TH2F("","",rise_time_bin,rise_time_xmin,rise_time_xmax,peak_bin,peak_xmin,peak_xmax);

    h2_lgad_charge_peak->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_lgad_charge_peak->GetYaxis()->SetLimits(peak_xmin,peak_xmax);

    h2_lgad_charge_rise_time->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_lgad_charge_rise_time->GetYaxis()->SetLimits(rise_time_xmin,rise_time_xmax);

    h_lgad_charge->SetLineColor(kRed);
    h_lgad_charge->SetLineStyle(kDashed);
    h_lgad_charge->SetMarkerColor(kRed);
    h_lgad_charge->SetLineWidth(2);

    h_lgad_baseline->SetLineColor(kRed);
    h_lgad_baseline->SetLineStyle(kDashed);
    h_lgad_baseline->SetMarkerColor(kRed);
    h_lgad_baseline->SetLineWidth(2);
    
    h_lgad_rise_edge->SetLineColor(kRed);
    h_lgad_rise_edge->SetLineStyle(kDashed);
    h_lgad_rise_edge->SetMarkerColor(kRed);
    h_lgad_rise_edge->SetLineWidth(2);
    
    h2_lgad_charge_rise_time->SetMarkerColor(kRed);
    h2_lgad_charge_peak->SetMarkerColor(kRed);
    h2_lgad_rise_time_peak->SetMarkerColor(kRed);


    h_lgad_charge->GetXaxis()->SetTitle("Collected Charges [fC]");
    h_lgad_charge->GetXaxis()->CenterTitle();
    h_lgad_charge->GetXaxis()->SetTitleOffset(1.5);
    h_lgad_charge->GetYaxis()->SetTitle("Normalized Entries");
    h_lgad_charge->GetYaxis()->CenterTitle();
    h_lgad_charge->GetYaxis()->SetTitleOffset(1.5);

    h_lgad_baseline->GetXaxis()->SetTitle("Baseline [mv]");
    h_lgad_baseline->GetXaxis()->CenterTitle();
    h_lgad_baseline->GetXaxis()->SetTitleOffset(1.5);
    h_lgad_baseline->GetYaxis()->SetTitle("Normalized Entries");
    h_lgad_baseline->GetYaxis()->CenterTitle();
    h_lgad_baseline->GetYaxis()->SetTitleOffset(1.5);

    h_lgad_rise_edge->GetXaxis()->SetTitle("Rising Edge [mv/ns]");
    h_lgad_rise_edge->GetXaxis()->CenterTitle();
    h_lgad_rise_edge->GetXaxis()->SetTitleOffset(1.5);
    h_lgad_rise_edge->GetYaxis()->SetTitle("Normalized Entries");
    h_lgad_rise_edge->GetYaxis()->CenterTitle();
    h_lgad_rise_edge->GetYaxis()->SetTitleOffset(1.5);

    h2_lgad_charge_rise_time->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_lgad_charge_rise_time->GetXaxis()->CenterTitle();
    h2_lgad_charge_rise_time->GetXaxis()->SetTitleOffset(1.5);
    h2_lgad_charge_rise_time->GetYaxis()->SetTitle("Rise Time [ns]");
    h2_lgad_charge_rise_time->GetYaxis()->CenterTitle();
    h2_lgad_charge_rise_time->GetYaxis()->SetTitleOffset(1.5);

    h2_lgad_charge_peak->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_lgad_charge_peak->GetXaxis()->CenterTitle();
    h2_lgad_charge_peak->GetXaxis()->SetTitleOffset(1.5);
    h2_lgad_charge_peak->GetYaxis()->SetTitle("Peak Amplitude [mV]");
    h2_lgad_charge_peak->GetYaxis()->CenterTitle();
    h2_lgad_charge_peak->GetYaxis()->SetTitleOffset(1.5);

    h2_lgad_rise_time_peak->GetXaxis()->SetTitle("Rise Time [ns]");
    h2_lgad_rise_time_peak->GetXaxis()->CenterTitle();
    h2_lgad_rise_time_peak->GetXaxis()->SetTitleOffset(1.5);
    h2_lgad_rise_time_peak->GetYaxis()->SetTitle("Peak Amplitude [mV]");
    h2_lgad_rise_time_peak->GetYaxis()->CenterTitle();
    h2_lgad_rise_time_peak->GetYaxis()->SetTitleOffset(1.5);

    TFile *lgad_f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-IMEV1beta-One_Chanel-w8-IV-E4-L1-15-70-beta-w8-IV-E4-L1-15_70_beta_180V-C3BiasV.root","READ");
    TTree *lgad_tree = (TTree *)lgad_f_in->Get("tree");

    double lgad_charge, lgad_peak, lgad_rise_time, lgad_baseline;

    lgad_tree->SetBranchAddress("charge",&lgad_charge);
    lgad_tree->SetBranchAddress("peak",&lgad_peak);
    lgad_tree->SetBranchAddress("rise_time",&lgad_rise_time);
    lgad_tree->SetBranchAddress("baseline",&lgad_baseline);

    for(int i=0; i<lgad_tree->GetEntries(); i++){

        lgad_tree->GetEntry(i);

        h_lgad_charge->Fill(lgad_charge);
        h_lgad_baseline->Fill(lgad_baseline);
        h_lgad_rise_edge->Fill(lgad_peak/lgad_rise_time);

        h2_lgad_charge_peak->Fill(lgad_charge,lgad_peak);
        h2_lgad_charge_rise_time->Fill(lgad_charge,lgad_rise_time);
        h2_lgad_rise_time_peak->Fill(lgad_rise_time,lgad_peak);

    }


    h_lgad_charge->Scale(1/h_lgad_charge->Integral());
    h_lgad_baseline->Scale(1/h_lgad_baseline->Integral());
    h_lgad_rise_edge->Scale(1/h_lgad_rise_edge->Integral());
    

    TH1F *h_pin_charge = new TH1F("","",40,0,20);
    TH1F *h_pin_baseline = new TH1F("","",70,-7.,7);
    TH1F *h_pin_rise_edge = new TH1F("","",250,0,500);
    TH2F *h2_pin_charge_peak = new TH2F("","",40,0,20,100,0,300);
    TH2F *h2_pin_charge_rise_time = new TH2F("","",40,0,20,30,0,3);
    TH2F *h2_pin_rise_time_peak = new TH2F("","",30,0,3,100,0,300);

    h2_pin_charge_peak->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_pin_charge_peak->GetYaxis()->SetLimits(peak_xmin,peak_xmax);

    h2_pin_charge_rise_time->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_pin_charge_rise_time->GetYaxis()->SetLimits(rise_time_xmin,rise_time_xmax);


    h_pin_charge->SetLineColor(417);
    h_pin_charge->SetLineStyle(kDashed);
    h_pin_charge->SetMarkerColor(kRed);
    h_pin_charge->SetLineWidth(2);

    h_pin_baseline->SetLineColor(417);
    h_pin_baseline->SetLineStyle(kDashed);
    h_pin_baseline->SetMarkerColor(417);
    h_pin_baseline->SetLineWidth(2);
    
    h_pin_rise_edge->SetLineColor(417);
    h_pin_rise_edge->SetLineStyle(kDashed);
    h_pin_rise_edge->SetMarkerColor(417);
    h_pin_rise_edge->SetLineWidth(2);
    
    h2_pin_charge_rise_time->SetMarkerColor(417);
    h2_pin_charge_peak->SetMarkerColor(417);
    h2_pin_rise_time_peak->SetMarkerColor(417);


    h_pin_charge->GetXaxis()->SetTitle("Collected Charges [fC]");
    h_pin_charge->GetXaxis()->CenterTitle();
    h_pin_charge->GetXaxis()->SetTitleOffset(1.5);
    h_pin_charge->GetYaxis()->SetTitle("Normalized Entries");
    h_pin_charge->GetYaxis()->CenterTitle();
    h_pin_charge->GetYaxis()->SetTitleOffset(1.5);

    h_pin_baseline->GetXaxis()->SetTitle("Baseline [mv]");
    h_pin_baseline->GetXaxis()->CenterTitle();
    h_pin_baseline->GetXaxis()->SetTitleOffset(1.5);
    h_pin_baseline->GetYaxis()->SetTitle("Normalized Entries");
    h_pin_baseline->GetYaxis()->CenterTitle();
    h_pin_baseline->GetYaxis()->SetTitleOffset(1.5);

    h2_pin_charge_rise_time->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_pin_charge_rise_time->GetXaxis()->CenterTitle();
    h2_pin_charge_rise_time->GetXaxis()->SetTitleOffset(1.5);
    h2_pin_charge_rise_time->GetYaxis()->SetTitle("Rise Time [ns]");
    h2_pin_charge_rise_time->GetYaxis()->CenterTitle();
    h2_pin_charge_rise_time->GetYaxis()->SetTitleOffset(1.5);

    h2_pin_charge_peak->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_pin_charge_peak->GetXaxis()->CenterTitle();
    h2_pin_charge_peak->GetXaxis()->SetTitleOffset(1.5);
    h2_pin_charge_peak->GetYaxis()->SetTitle("Peak Amplitude [mV]");
    h2_pin_charge_peak->GetYaxis()->CenterTitle();
    h2_pin_charge_peak->GetYaxis()->SetTitleOffset(1.5);


    TFile *pin_f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-IMEV1beta-One_Chanel-w8-IV-E4-PIN-15-100-beta-250V-C3BiasV.root","READ");
    TTree *pin_tree = (TTree *)pin_f_in->Get("tree");

    double pin_charge, pin_peak, pin_rise_time, pin_baseline;

    pin_tree->SetBranchAddress("charge",&pin_charge);
    pin_tree->SetBranchAddress("peak",&pin_peak);
    pin_tree->SetBranchAddress("rise_time",&pin_rise_time);
    pin_tree->SetBranchAddress("baseline",&pin_baseline);

    for(int i=0; i<pin_tree->GetEntries(); i++){

        pin_tree->GetEntry(i);

        h_pin_charge->Fill(pin_charge);
        h_pin_baseline->Fill(pin_baseline);
        h_pin_rise_edge->Fill(pin_peak/pin_rise_time);

        h2_pin_charge_peak->Fill(pin_charge,pin_peak);
        h2_pin_charge_rise_time->Fill(pin_charge,pin_rise_time);
        h2_pin_rise_time_peak->Fill(pin_rise_time,pin_peak);

    }


    h_pin_charge->Scale(1/h_pin_charge->Integral());
    h_pin_baseline->Scale(1/h_pin_baseline->Integral());
    h_pin_rise_edge->Scale(1/h_pin_rise_edge->Integral());




    TH1F *h_sic_charge = new TH1F("","",40,0,20);
    TH1F *h_sic_baseline = new TH1F("","",70,-7.,7);
    TH1F *h_sic_rise_edge = new TH1F("","",250,0,500);
    TH2F *h2_sic_charge_peak = new TH2F("","",40,0,20,100,0,300);
    TH2F *h2_sic_charge_rise_time = new TH2F("","",40,0,20,30,0,3);
    TH2F *h2_sic_rise_time_peak = new TH2F("","",30,0,3,100,0,300);

    h2_sic_charge_peak->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_sic_charge_peak->GetYaxis()->SetLimits(peak_xmin,peak_xmax);

    h2_sic_charge_rise_time->GetXaxis()->SetLimits(charge_xmin,charge_xmax);
    h2_sic_charge_rise_time->GetYaxis()->SetLimits(rise_time_xmin,rise_time_xmax);


    h_sic_charge->SetLineColor(kCyan);
    h_sic_charge->SetLineStyle(kDashed);
    h_sic_charge->SetMarkerColor(kCyan);
    h_sic_charge->SetLineWidth(2);

    h_sic_baseline->SetLineColor(kCyan);
    h_sic_baseline->SetLineStyle(kDashed);
    h_sic_baseline->SetMarkerColor(kCyan);
    h_sic_baseline->SetLineWidth(2);
    
    h_sic_rise_edge->SetLineColor(kCyan);
    h_sic_rise_edge->SetLineStyle(kDashed);
    h_sic_rise_edge->SetMarkerColor(kCyan);
    h_sic_rise_edge->SetLineWidth(2);
    
    h2_sic_charge_rise_time->SetMarkerColor(kCyan);
    h2_sic_charge_peak->SetMarkerColor(kCyan);
    h2_sic_rise_time_peak->SetMarkerColor(kCyan);


    h_sic_charge->GetXaxis()->SetTitle("Collected Charges [fC]");
    h_sic_charge->GetXaxis()->CenterTitle();
    h_sic_charge->GetXaxis()->SetTitleOffset(1.5);
    h_sic_charge->GetYaxis()->SetTitle("Normalized Entries");
    h_sic_charge->GetYaxis()->CenterTitle();
    h_sic_charge->GetYaxis()->SetTitleOffset(1.5);

    h_sic_baseline->GetXaxis()->SetTitle("Baseline [mv]");
    h_sic_baseline->GetXaxis()->CenterTitle();
    h_sic_baseline->GetXaxis()->SetTitleOffset(1.5);
    h_sic_baseline->GetYaxis()->SetTitle("Normalized Entries");
    h_sic_baseline->GetYaxis()->CenterTitle();
    h_sic_baseline->GetYaxis()->SetTitleOffset(1.5);

    h2_sic_charge_rise_time->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_sic_charge_rise_time->GetXaxis()->CenterTitle();
    h2_sic_charge_rise_time->GetXaxis()->SetTitleOffset(1.5);
    h2_sic_charge_rise_time->GetYaxis()->SetTitle("Rise Time [ns]");
    h2_sic_charge_rise_time->GetYaxis()->CenterTitle();
    h2_sic_charge_rise_time->GetYaxis()->SetTitleOffset(1.5);

    h2_sic_charge_peak->GetXaxis()->SetTitle("Collected Charges [fC]");
    h2_sic_charge_peak->GetXaxis()->CenterTitle();
    h2_sic_charge_peak->GetXaxis()->SetTitleOffset(1.5);
    h2_sic_charge_peak->GetYaxis()->SetTitle("Peak Amplitude [mV]");
    h2_sic_charge_peak->GetYaxis()->CenterTitle();
    h2_sic_charge_peak->GetYaxis()->SetTitleOffset(1.5);


    TFile *sic_f_in = new TFile("./output/WaveAnalysis_-home-admin-STDB-BetaTest-data-4H-SiC-beta-4H_SiC_PIN_1-500V-C3BiasV.root","READ");
    TTree *sic_tree = (TTree *)sic_f_in->Get("tree");

    double sic_charge, sic_peak, sic_rise_time, sic_baseline;

    sic_tree->SetBranchAddress("charge",&sic_charge);
    sic_tree->SetBranchAddress("peak",&sic_peak);
    sic_tree->SetBranchAddress("rise_time",&sic_rise_time);
    sic_tree->SetBranchAddress("baseline",&sic_baseline);

    for(int i=0; i<sic_tree->GetEntries(); i++){

        sic_tree->GetEntry(i);

        h_sic_charge->Fill(sic_charge);
        h_sic_baseline->Fill(sic_baseline);
        h_sic_rise_edge->Fill(sic_peak/sic_rise_time);

        h2_sic_charge_peak->Fill(sic_charge,sic_peak);
        h2_sic_charge_rise_time->Fill(sic_charge,sic_rise_time);
        h2_sic_rise_time_peak->Fill(sic_rise_time,sic_peak);

    }


    h_sic_charge->Scale(1/h_sic_charge->Integral());
    h_sic_baseline->Scale(1/h_sic_baseline->Integral());
    h_sic_rise_edge->Scale(1/h_sic_rise_edge->Integral());


    h_lgad_charge->GetYaxis()->SetRangeUser(0.,0.3);


    TF1 *landau_conv_gaus = new TF1("landau_lgad","gaus*landau",0,80);


    TF1 *gaus_lgad = new TF1("gaus_lgad","gaus",-7,7);
    TF1 *gaus_pin = new TF1("gaus_pin","gaus",-7,7);
    TF1 *gaus_sic = new TF1("gaus_sic","gaus",-7,7);

    h_lgad_baseline->Fit("gaus_lgad","R0+");
    h_pin_baseline->Fit("gaus_pin","R0+");
    h_sic_baseline->Fit("gaus_sic","R0+");

    gaus_lgad->SetLineColor(kRed);
    gaus_pin->SetLineColor(417);
    gaus_sic->SetLineColor(kCyan);



    TLegend *l1 = new TLegend(0.55,0.7,0.9,0.85);
    l1->AddEntry(h_lgad_charge,"U=100V Si LGAD 50um ");
    l1->AddEntry(h_pin_charge,"U=250V Si PIN  50um");
    l1->AddEntry(h_sic_charge,"U=500V 4H-SiC PIN  100um");
    
    TLegend *l2 = new TLegend(0.6,0.7,0.85,0.85);
    l2->AddEntry(h_lgad_charge,"Si LGAD 50um ");
    l2->AddEntry(h_pin_charge,"Si PIN  50um");
    l2->AddEntry(h_sic_charge,"4H-SiC PIN  50um");



    TCanvas *c1 = new TCanvas("c1","c1",500,500);
    c1->SetGrid();
    c1->cd();
    h_lgad_charge->Draw("hist");
    h_pin_charge->Draw("same hist");
    h_sic_charge->Draw("same hist");
    //landau_lgad->Draw("same");
    l1->Draw("same");

    TCanvas *c2 = new TCanvas("c2","c2",500,500);
    c2->SetGrid();
    c2->cd();
    h_sic_baseline->Draw();
    h_lgad_baseline->Draw("same hist");
    h_pin_baseline->Draw("same hist");
    
    gaus_lgad->Draw("same");
    gaus_pin->Draw("same");
    gaus_sic->Draw("same");

    l1->Draw("same");


    TCanvas *c3 = new TCanvas("c3","c3",500,500);
    c3->SetGrid();
    c3->cd();
    h2_lgad_charge_rise_time->Draw();
    h2_pin_charge_rise_time->Draw("same");
    h2_sic_charge_rise_time->Draw("same");
    l2->Draw("same");


    TCanvas *c4 = new TCanvas("c4","c4",500,500);
    c4->SetGrid();
    c4->cd();
    h2_lgad_charge_peak->Draw();
    h2_pin_charge_peak->Draw("same");
    h2_sic_charge_peak->Draw("same");
    l2->Draw("same");

    TCanvas *c5 = new TCanvas("c5","c5",500,500);
    c5->SetGrid();
    c5->cd();
    h2_lgad_rise_time_peak->Draw();
    h2_pin_rise_time_peak->Draw("same");
    h2_sic_rise_time_peak->Draw("same");
    l2->Draw("same");

    TCanvas *c6 = new TCanvas("c6","c6",500,500);
    c6->SetGrid();
    c6->cd();
    h_lgad_rise_edge->Draw("hist");
    h_pin_rise_edge->Draw("same hist");
    h_sic_rise_edge->Draw("same hist");
    l1->Draw("same");
    

    c1->SaveAs("./fig/compare_charge.pdf");
    c2->SaveAs("./fig/compare_baseline.pdf");
    c3->SaveAs("./fig/compare_charge_rise_time.png");
    c4->SaveAs("./fig/compare_charge_peak.png");
    c5->SaveAs("./fig/compare_rise_time_peak.png");
    c6->SaveAs("./fig/compare_rise_edge.pdf");


}

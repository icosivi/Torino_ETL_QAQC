//-------c++----------------//
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <numeric>
#include <functional>
#include <sys/stat.h>
#include <dirent.h>

//------ROOT----------------//
#include <TTree.h>
#include <TTreeReader.h>
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TColor.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>


void categories(){

  gStyle->SetOptStat(0);

  int palette[3];
  palette[0] = 3;
  palette[1] = 5;
  palette[2] = 2;
  
  gStyle->SetTitleOffset( 1.3, "z" );
  gStyle->SetLabelOffset( 0., "z" );
  //gStyle->SetTitleSize(0.06,"z"); 
  gStyle->SetLabelSize(0.035,"z");

  gStyle->SetPaintTextFormat("1.1f");

   
  ///////////////////////////////////////////////////////////////////// 
  // THIS SCRIPT ASSUMES raw voltages in [V] and raw currents in [A] //
  // Parameters below are in [V] and [uA]                            //
  /////////////////////////////////////////////////////////////////////


  float V_current_level = 230. ; // Leakage current measured at this V
  float V_current_monitor = 200. ; //Voltage up to which the current level is monitored
  float VBD_expected = 230. ; //minimum VBD to be considered for GOOD or MEDIUM categories, if VBD<VBD_expected sensor is BAD
  float V_min_kfactor = 100. ; //k-factor not used to calculate VBD if VBD<V_min_kfactor
  float I_thr = 10.; //sensor discarded if I > I_thr [uA] in the voltage operation range [ 0-V_current_monitor ]
  float I_compliance = 1000; // VBD calculation begins when I < I_compliance  [uA]
  float I_compliance_minimum = 50; // [uA] VBD calculation performed only if compliance was set above this threshold
  float k_thr = 8.; // 8  k value to define VBD using k-factor method
  float current_conversion_value = 1E6; // conversion from [A] (raw data) to [uA] (used in the final plots)
  int start_bd_calculation = 5; //BD calculation start from this sampled bias point: avoid considering the very first voltages of the bias sweep 


  float low_iv_range = 0.1;  //low and high ranges for I@100V plot [uA]
  float high_iv_range = 500;
  float low_vbd_range = 0; //low and high ranges for VBD plot [V]
  float high_vbd_range = 300;

  bool bcurrent = false;
  bool bvoltage = true;
  bool bcategory = false;
  bool bnoisy = false;
  bool save = false;

  bool verbose = false;

  if(bcategory) gStyle->SetPalette(3,palette); // custom palette used for Categories
  
  if(bvoltage) gStyle->SetPalette(kBlackBody); // palette for VBD

  if(bcurrent || bnoisy){ 
  gStyle->SetPalette(kBlackBody); // palette for I@V_current_level
  TColor::InvertPalette();
  }


  bool invert_polarity=true;

  TFile *file_qa;
  TTree *tree_qa;
  
  file_qa = TFile::Open("../../root_files/UFSD4_16x16_IVtree.root");
  tree_qa = dynamic_cast<TTree*>(file_qa->Get("Tree"));
  TTreeReader reader_qa("Tree", file_qa);
  

  // This script works with root files having the following branches. IBACK is the total current measured by the backplane
  TTreeReaderArray<float> I_qa(reader_qa, "IBACK");
  TTreeReaderArray<float> V_qa(reader_qa, "V");
  TTreeReaderValue<int> wafer_qa(reader_qa, "wafer");
  TTreeReaderValue<int> col_qa(reader_qa, "column");
  TTreeReaderValue<int> row_qa(reader_qa, "row");
  
  
  int counter_qa = 0;
  float I_qa_100V = 0.;
  float VBD_qa_u = 0.;
  float VBD_qa_d = 0.;
  float k_qa_u[2] = {0.,0.};
  float k_qa_d[2] = {0.,0.};

  int current_levels_counter[6]={0,0,0,0,0,0};

  float vcount_vbd_qa[18][5][5];
  float vcount_i_qa[18][5][5];
  float vcount_bump_qa[18][5][5];
  float vvbd_qa[18][5][5];
  float vi_qa[18][5][5];
  float vbump_qa[18][5][5];
  //float vnoisy_qa[18][5][5];
  TH2F *hI_qa_100V[18];
  TH2F *hV_qa[18];
  TH2F *hbump_qa[18];
  TH2F *hnoisy_qa[18];
  TH2F *hcount_qa[18];
  float min_I[18];

  for(int i=0; i<18; i++){
    for(int j=0; j<5; j++){
      for(int k=0; k<5; k++){
        vbump_qa[i][j][k]=0;
        vcount_vbd_qa[i][j][k]=0;
        vcount_i_qa[i][j][k]=0;
        vcount_bump_qa[i][j][k]=0;
        vvbd_qa[i][j][k]=0;
        vi_qa[i][j][k]=0;
        //vnoisy_qa[i][j][k]=0;
      }
    }

    hI_qa_100V[i]=new TH2F( Form("I_qa_%iV_W%i",int(V_current_level),i+1), Form("I@%iV on-wafer W%i",int(V_current_level),i+1), 5,2,7,5,2,7);
    hV_qa[i]=new TH2F( Form("V_qa_W%i",i+1), Form("VBD on-wafer W%i",i+1), 5,2,7,5,2,7);
    min_I[i] = 1000;

    hbump_qa[i]=new TH2F( Form("bump_qa_W%i",i+1), Form("Categories on-wafer W%i",i+1), 5,2,7,5,2,7);
    hnoisy_qa[i]=new TH2F( Form("noisy_qa_W%i",i+1), Form("Bad pads on-wafer W%i",i+1), 5,2,7,5,2,7);

  }


  int dumb_counter=0; // Total number of sensors/measurements in the root file
  int non_empty_counter=0; // Sensors with both I and V arrays having size !=0 (the sensor was in fact measured)
  int non_zero_counter=0; // Sensors with both I and V arrays having values !=0 (the sensor was measured and can be biased)

  int total_sensors_counter_qa=0; //Fraction of sensors whose category can be properly defined
  int good_sensors_counter_qa=0; //Fraction of GOOD sensors
  int bad_sensors_counter_qa=0; //Fraction of BAD sensors
  int medium_sensors_counter_qa=0; //Fraction of MEDIUM sensors

  bool iv_quality; // Boolean reflecting the quality of the IV curve: sensor is BAD whenever iv_quality is set to False

  while( reader_qa.Next() ){

    dumb_counter++ ;

    if(int(I_qa.GetSize())>=start_bd_calculation && int(V_qa.GetSize())>=start_bd_calculation){
  
       non_empty_counter++ ; 

      float non_zero_current = 0;
      for(int i=1;  i<int(I_qa.GetSize()); i++) non_zero_current += I_qa.At(i) ;

      float non_zero_voltage = 0;
      for(int i=1;  i<int(V_qa.GetSize()); i++) non_zero_voltage += V_qa.At(i) ;


      if(non_zero_voltage!=0 && non_zero_current!=0){

        non_zero_counter++ ;
    
        iv_quality = true;
        I_qa_100V = -1000;
        VBD_qa_u = -1000.;
        VBD_qa_d = -1000.;
        
        int type_counter = 0;
        int type_counter_I_qa = 0;
        int type_counter_V_qa_u = 0;
        int type_counter_V_qa_d = 0;

        k_qa_u[0] = 0.;
        k_qa_u[1] = 0.;
        k_qa_d[0] = 0.;
        k_qa_d[1] = 0.;
      
          
          if( I_qa.At(I_qa.GetSize()-1) ){ 

            for(int i=0; i<int(I_qa.GetSize()); i++) I_qa.At(i) = -I_qa.At(i);
            if(verbose) cout<<"Sensor from wafer "<<*wafer_qa<<" row "<<*row_qa<<" column "<<*col_qa<<" with negative current detected. Changing current polarity."<<endl;

          }
      
          
          //////////// I@V_current_level calculation ///////////////////////
          for(int i=0; i<int(I_qa.GetSize()); i++){
      
          	if( i!=0 && V_qa.At(i-1)<V_current_level && V_qa.At(i)>=V_current_level ){
      
              if( V_qa.At(i)==V_current_level ) I_qa_100V = I_qa.At(i);
              else I_qa_100V = (((I_qa.At(i)-I_qa.At(i-1))/(V_qa.At(i)-V_qa.At(i-1)))*(V_current_level-V_qa.At(i)) + I_qa.At(i));
              break;
          	}
      
          }
      
          for(int i=1;  i<int(I_qa.GetSize()); i++){
      
            if( V_qa.At(i)<=V_current_monitor && I_qa.At(i) > I_thr*(1./current_conversion_value) ){
      
              iv_quality=false;
              break;
    
            }
          }
    
          I_qa_100V = current_conversion_value*I_qa_100V;
          if( I_qa_100V<0. ) iv_quality=false;
          
          //////////////// VBD calculation //////////////////////////////
      
          if( I_qa.At( I_qa.GetSize()-1 )>I_compliance_minimum*(1./current_conversion_value) ){
      
            
            ///////// VBD "up" (Calculation start from the end of the Voltage array downwards) /////////
            for(int i=(int(I_qa.GetSize())-2); i>=start_bd_calculation; i--){
        
              if(V_qa.At( I_qa.GetSize()-1 )<V_min_kfactor){
                  
                  VBD_qa_u = V_qa.At( I_qa.GetSize()-1 ) ;
                  break ;
              }
              
              if( I_qa.At(i)<I_compliance*(1./current_conversion_value) && V_qa.At(i)>=V_min_kfactor ){
        
                k_qa_u[0] = ( (I_qa.At(i)-I_qa.At(i-1))/(V_qa.At(i)-V_qa.At(i-1)) )*(V_qa.At(i)/I_qa.At(i)) ;
                k_qa_u[1] = ( (I_qa.At(i+1)-I_qa.At(i))/(V_qa.At(i+1)-V_qa.At(i)) )*(V_qa.At(i)/I_qa.At(i)) ;
                if(*wafer_qa==1 && *row_qa==4 && *col_qa==3) cout<<k_qa_u[0]<<" "<<k_qa_u[1]<<" "<<V_qa.At(i)<<endl;
        
                if( k_qa_u[0]<k_thr && k_qa_u[1]>=k_thr ){
                  
                  VBD_qa_u = V_qa.At(i) ;
                  break ;
                }
              }
            }
        
      
            ///////// VBD "down" (Calculation start from the beginning of the Voltage array upwards)/////////
            for(int i=start_bd_calculation; i<=(int(I_qa.GetSize())-2); i++){
        
              if(V_qa.At( I_qa.GetSize()-1 )<V_min_kfactor){
                  
                  VBD_qa_d = V_qa.At( I_qa.GetSize()-1 ) ;
                  break ;
              }
        
              if( I_qa.At(i)<I_compliance*(1./current_conversion_value) && V_qa.At(i)>=V_min_kfactor ){
                  
                k_qa_d[0] = ( (I_qa.At(i)-I_qa.At(i-1))/(V_qa.At(i)-V_qa.At(i-1)) )*(V_qa.At(i)/I_qa.At(i)) ;
                k_qa_d[1] = ( (I_qa.At(i+1)-I_qa.At(i))/(V_qa.At(i+1)-V_qa.At(i)) )*(V_qa.At(i)/I_qa.At(i)) ;
        
                if( k_qa_d[0]<k_thr && k_qa_d[1]>=k_thr ){
                  
                  VBD_qa_d = V_qa.At(i) ;
                  break ;
                }
              }
            }
          }
      
      
      
          if(VBD_qa_u!=-1000 && VBD_qa_d!=-1000){ 
            
            vvbd_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += VBD_qa_u;
            vcount_vbd_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 1;
      
          }
      
          if(I_qa_100V>=0){ 
      
            vi_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += I_qa_100V;
            vcount_i_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 1;
      
          }
      
         
          if( VBD_qa_u>VBD_expected && VBD_qa_d!=-1000 &&  VBD_qa_d<=VBD_expected ){ 

            if(verbose) cout<<"!!! SENSOR WITH CURRENT JUMPS ALERT  !!!: "<<"from wafer "<<*wafer_qa<<" row "<<*row_qa<<" column "<<*col_qa<<endl;
            //vnoisy_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 1;

          }
          
          
          if( iv_quality && VBD_qa_u>VBD_expected && VBD_qa_d!=-1000 ){

            if( VBD_qa_d>VBD_expected ){
      
              vbump_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 1; 
              vcount_bump_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 1; 
      
            }else{ 
      
              vbump_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 2;
              vcount_bump_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 1;
      
            }
      
          }else if( !iv_quality || (VBD_qa_d!=-1000 && VBD_qa_u!=-1000 && VBD_qa_u<=VBD_expected ) ){
      
            vbump_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 100;
            vcount_bump_qa[*wafer_qa-1][*col_qa-2][*row_qa-2] += 1;
      
          }else if(verbose) cout<<"Sensor from wafer "<<*wafer_qa<<" row "<<*row_qa<<" column "<<*col_qa<<" has Current within acceptance, but VBD could not be calculated. PLEASE CHECK."<<endl;
      }else if(verbose) cout<<"Sensor from wafer "<<*wafer_qa<<" row "<<*row_qa<<" column "<<*col_qa<<" has Voltage and/or Current always equal to zero. PLEASE CHECK."<<endl;
    }else if(verbose) cout<<"Sensor from wafer "<<*wafer_qa<<" row "<<*row_qa<<" column "<<*col_qa<<Form(": Voltage and/or Current vectors have less than %i elements. PLEASE CHECK.",start_bd_calculation)<<endl;
  }



  for(int i=0; i<18; i++){
    for(int j=0; j<5; j++){
      for(int k=0; k<5; k++){

        if( i!=9 && i!=10 ){
  
          if(vcount_vbd_qa[i][j][k]!=0) hV_qa[i]->Fill(j+2,k+2, vvbd_qa[i][j][k]/vcount_vbd_qa[i][j][k] );
          if(vcount_i_qa[i][j][k]!=0) hI_qa_100V[i]->Fill(j+2,k+2, vi_qa[i][j][k]/vcount_i_qa[i][j][k] );

          //if(vnoisy_qa[i][j][k]!=0) hnoisy_qa[i]->Fill(j+2,k+2,1);
          //else hnoisy_qa[i]->Fill(j+2,k+2,0);
  
          if(vcount_bump_qa[i][j][k]!=0) total_sensors_counter_qa++;

          if(vcount_bump_qa[i][j][k]==0) hnoisy_qa[i]->Fill(j+2,k+2,-1000);
          
          if( vbump_qa[i][j][k]/vcount_bump_qa[i][j][k]==1 ){
          
            hbump_qa[i]->Fill(j+2,k+2,1);
            hnoisy_qa[i]->Fill(j+2,k+2,0);
            good_sensors_counter_qa++;
  
          }else if( vbump_qa[i][j][k]/vcount_bump_qa[i][j][k]>2 ){
  
            hbump_qa[i]->Fill(j+2,k+2,3);
            hnoisy_qa[i]->Fill(j+2,k+2,0);
            bad_sensors_counter_qa++;
  
  
          }else if( vbump_qa[i][j][k]/vcount_bump_qa[i][j][k]>1 && vbump_qa[i][j][k]/vcount_bump_qa[i][j][k]<=2 ){
  
            hbump_qa[i]->Fill(j+2,k+2,2);
            hnoisy_qa[i]->Fill(j+2,k+2,1);
            medium_sensors_counter_qa++;
  
          }

          if( vi_qa[i][j][k]/vcount_i_qa[i][j][k]<10. && vi_qa[i][j][k]/vcount_i_qa[i][j][k]>=0. ) current_levels_counter[0]++ ;
          else if( vi_qa[i][j][k]/vcount_i_qa[i][j][k]>=10. && vi_qa[i][j][k]/vcount_i_qa[i][j][k]<50. ) current_levels_counter[1]++ ;
          else if( vi_qa[i][j][k]/vcount_i_qa[i][j][k]>=50. && vi_qa[i][j][k]/vcount_i_qa[i][j][k]<100. ) current_levels_counter[2]++ ;
          else if( vi_qa[i][j][k]/vcount_i_qa[i][j][k]>=100. && vi_qa[i][j][k]/vcount_i_qa[i][j][k]<500. ) current_levels_counter[3]++ ;
          else if( vi_qa[i][j][k]/vcount_i_qa[i][j][k]>=500. ) current_levels_counter[4]++ ;
          else if( vcount_i_qa[i][j][k]==0. && vcount_vbd_qa[i][j][k]!=0. )  current_levels_counter[5]++ ;

        }
      }
    }
  }


  if(verbose){
  
    cout<<"Fraction of GOOD sensors (on-wafer): "<<double(good_sensors_counter_qa)/double(total_sensors_counter_qa)<<endl;
    cout<<"Fraction of MEDIUM sensors (on-wafer): "<<double(medium_sensors_counter_qa)/double(total_sensors_counter_qa)<<endl;
    cout<<"Fraction of BAD sensors (on-wafer): "<<double(bad_sensors_counter_qa)/double(total_sensors_counter_qa)<<endl;
    cout<<"\n";
  
  }

  /*cout<<"Fraction with I<10 uA: "<<double(current_levels_counter[0])/double(total_sensors_counter_qa)<<endl;
  cout<<"Fraction with I in 10-50 uA range: "<<double(current_levels_counter[1])/double(total_sensors_counter_qa)<<endl;
  cout<<"Fraction with I in 50-100 uA range: "<<double(current_levels_counter[2])/double(total_sensors_counter_qa)<<endl;
  cout<<"Fraction with I in 100-500 uA range: "<<double(current_levels_counter[3])/double(total_sensors_counter_qa)<<endl;
  cout<<"Fraction with I>500 uA: "<<double(current_levels_counter[4])/double(total_sensors_counter_qa)<<endl;
  cout<<"Fraction not reaching 100V: "<<double(current_levels_counter[5])/double(total_sensors_counter_qa)<<endl;
  cout<<"\n";*/
 
      
  TCanvas *cI_100V[18];
  TCanvas *cV[18];
  TCanvas *cbump[18];
  TCanvas *cnoisy[18];


  for(int i=0; i<18; i++){
   
    //if( i!=9 && i!=10 ){
    //if( i==4 || i==5 || i==6 || i==7 || i==8 ){
    if( i==0 ){

      //hI_qa_100V[i]->GetZaxis()->SetRangeUser( 0.1, hI_qa_100V[i]->GetMaximum() ); //Alternative colored axis range
      hI_qa_100V[i]->GetZaxis()->SetRangeUser( low_iv_range, high_iv_range );
      hI_qa_100V[i]->GetZaxis()->SetTitle("[uA]");
      hI_qa_100V[i]->GetXaxis()->SetTitle("column"); 
      hI_qa_100V[i]->GetYaxis()->SetTitle("row"); 
      hI_qa_100V[i]->SetMarkerSize(3.);

  
      //hV_qa[i]->GetZaxis()->SetRangeUser( hV_qa[i]->GetMinimum(), hV_qa[i]->GetMaximum() ); //Alternative colored axis range
      hV_qa[i]->GetZaxis()->SetRangeUser( low_vbd_range, high_vbd_range ); 
      hV_qa[i]->GetZaxis()->SetTitle("[V]");  
      hV_qa[i]->GetXaxis()->SetTitle("column"); 
      hV_qa[i]->GetYaxis()->SetTitle("row"); 
      hV_qa[i]->SetMarkerSize(3.);


      hbump_qa[i]->GetZaxis()->SetRangeUser( 1, 3 );
      hbump_qa[i]->GetZaxis()->SetTitle("Category");
      hbump_qa[i]->GetXaxis()->SetTitle("column");
      hbump_qa[i]->GetYaxis()->SetTitle("row");


      hnoisy_qa[i]->GetZaxis()->SetRangeUser( -0.1, 1 );
      hnoisy_qa[i]->GetZaxis()->SetTitle("Presence of bad pad(s)");
      hnoisy_qa[i]->GetXaxis()->SetTitle("column");
      hnoisy_qa[i]->GetYaxis()->SetTitle("row");


      if(bcurrent){

        cI_100V[i]=new TCanvas(Form("c_I_%iV_W%i",int(V_current_level),i+1), Form("c I@%iV_W%i",int(V_current_level),i+1), 1000,1000);
        cI_100V[i]->SetRightMargin(0.15);
        cI_100V[i]->cd();
        hI_qa_100V[i]->Draw("textcolz");
        gPad->SetGrid(1,1);
        gPad->SetLogz(1);
        gPad->Update();
        cI_100V[i]->Update();
        hI_qa_100V[i]->GetXaxis()->SetNdivisions(4);
        hI_qa_100V[i]->GetYaxis()->SetNdivisions(3);
  
        if(save) cI_100V[i]->SaveAs( Form("pics/I_%iV_W%i_FINAL.png",int(V_current_level),i+1) );

      }


      if(bvoltage){

        gStyle->SetPaintTextFormat("1.1f");
        cV[i]=new TCanvas(Form("c_V_W%i",i+1), Form("c VBD_W%i",i+1), 1000,1000);
        cV[i]->SetRightMargin(0.15);
        cV[i]->cd();
        hV_qa[i]->Draw("textcolz");
        gPad->SetGrid(1,1);
        gPad->Update();
        cV[i]->Update();
        hV_qa[i]->GetXaxis()->SetNdivisions(4);
        hV_qa[i]->GetYaxis()->SetNdivisions(3);
  
        if(save) cV[i]->SaveAs( Form("pics/VBD_W%i_FINAL.png",i+1) );

      }


      if(bcategory){
      
        cbump[i]=new TCanvas(Form("c_bump_W%i",i+1), Form("c bump_W%i",i+1), 1000,1000);
        cbump[i]->SetRightMargin(0.15);
        cbump[i]->cd();
        hbump_qa[i]->Draw("colz");
        gPad->SetGrid(1,1);
        gPad->Update();
        cbump[i]->Update();
        hbump_qa[i]->GetXaxis()->SetNdivisions(4);
        hbump_qa[i]->GetYaxis()->SetNdivisions(3);
       
        if(save) cbump[i]->SaveAs( Form("pics/categories_W%i_FINAL.png",i+1) );   

      }


      if(bnoisy){
      
        cnoisy[i]=new TCanvas(Form("c_noisy_W%i",i+1), Form("c noisy_W%i",i+1), 1000,1000);
        cnoisy[i]->SetRightMargin(0.15);
        cnoisy[i]->cd();
        hnoisy_qa[i]->Draw("colz");
        gPad->SetGrid(1,1);
        gPad->Update();
        cnoisy[i]->Update();
        hnoisy_qa[i]->GetXaxis()->SetNdivisions(4);
        hnoisy_qa[i]->GetYaxis()->SetNdivisions(3);
       
        if(save) cnoisy[i]->SaveAs( Form("pics/noisy_W%i_FINAL.png",i+1) );

      }

    }

  }

}
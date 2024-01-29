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
  //gStyle->SetPalette(3,palette); // custom palette used for Categories
  
  gStyle->SetPalette(kBlackBody); // palette for I@V_current_level && VBD
  TColor::InvertPalette(); // Uncomment when plotting I@V_current_level

  gStyle->SetTitleOffset( 1.3, "z" );
  gStyle->SetLabelOffset( 0., "z" );
  //gStyle->SetTitleSize(0.06,"z"); 
  gStyle->SetLabelSize(0.035,"z");

  gStyle->SetPaintTextFormat("1.1f");

   
  ///////////////////////////////////////////////////////////////////// 
  // THIS SCRIPT ASSUMES raw voltages in [V] and raw currents in [A] //
  // Parameters below are in [V] and [uA]                            //
  /////////////////////////////////////////////////////////////////////


  float V_current_level = 100. ; //k-factor not used to calculate VBD if VBD<V_current_level + Leakage current measured at this V
  float V_current_monitor = 200. ; //Voltage up to which the current level is monitored
  float VBD_expected = 200. ; //minimum VBD to be considered for GOOD or MEDIUM categories, if VBD<VBD_expected sensor is BAD
  float V_difference = 20. ; //max difference in V between VBD_up and VBD_down for a sensor to be considered GOOD
  float I_thr = 10.; //sensor discarded if I > I_thr [uA] in the voltage operation range [ 0-V_current_monitor ]
  float I_compliance = 1000; // VBD calculation begins when I < I_compliance  [uA]
  float I_compliance_minimum = 50; // [uA] VBD calculation performed only if compliance was set above this threshold
  float k_thr = 8.; //k value to define VBD using k-factor method
  float current_conversion_value = 1E6; // conversion from [A] (raw data) to [uA] (used in the final plots)
  int start_bd_calculation = 5; //BD calculation start from this sampled bias point: avoid considering the very first voltages of the bias sweep 


  float low_iv_range = 0.1;  //low and high ranges for I@100V plot [uA]
  float high_iv_range = 500;
  float low_vbd_range = 0; //low and high ranges for VBD plot [V]
  float high_vbd_range = 300;

  bool bcurrent = true;
  bool bvoltage = false;
  bool bcategory = false;
  bool save = true;

  bool invert_polarity=true;

  TFile *file_fbk;
  TTree *tree_fbk;
  
  file_fbk = TFile::Open("root_files/UFSD4_16x16_IVtree.root");
  tree_fbk = dynamic_cast<TTree*>(file_fbk->Get("Tree"));
  TTreeReader reader_fbk("Tree", file_fbk);
  

  // This script works with root files having the following branches. IBACK is the total current measured by the backplane
  TTreeReaderArray<float> I_fbk(reader_fbk, "IBACK");
  TTreeReaderArray<float> V_fbk(reader_fbk, "V");
  TTreeReaderValue<int> wafer_fbk(reader_fbk, "wafer");
  TTreeReaderValue<int> col_fbk(reader_fbk, "column");
  TTreeReaderValue<int> row_fbk(reader_fbk, "row");
  
  
  int counter_fbk = 0;
  float I_fbk_100V = 0.;
  float VBD_fbk_u = 0.;
  float VBD_fbk_d = 0.;
  float k_fbk_u[2] = {0.,0.};
  float k_fbk_d[2] = {0.,0.};

  int current_levels_counter[6]={0,0,0,0,0,0};

  float vcount_vbd_fbk[18][5][5];
  float vcount_i_fbk[18][5][5];
  float vcount_bump_fbk[18][5][5];
  float vvbd_fbk[18][5][5];
  float vi_fbk[18][5][5];
  float vbump_fbk[18][5][5];
  TH2F *hI_fbk_100V[18];
  TH2F *hV_fbk[18];
  TH2F *hbump_fbk[18];
  TH2F *hcount_fbk[18];
  float min_I[18];

  for(int i=0; i<18; i++){
    for(int j=0; j<5; j++){
      for(int k=0; k<5; k++){
        vbump_fbk[i][j][k]=0;
        vcount_vbd_fbk[i][j][k]=0;
        vcount_i_fbk[i][j][k]=0;
        vcount_bump_fbk[i][j][k]=0;
        vvbd_fbk[i][j][k]=0;
        vi_fbk[i][j][k]=0;
      }
    }

    hI_fbk_100V[i]=new TH2F( Form("I_fbk_%iV_W%i",int(V_current_level),i+1), Form("I@%iV on-wafer W%i",int(V_current_level),i+1), 5,2,7,5,2,7);
    hV_fbk[i]=new TH2F( Form("V_fbk_W%i",i+1), Form("VBD on-wafer W%i",i+1), 5,2,7,5,2,7);
    min_I[i] = 1000;

    hbump_fbk[i]=new TH2F( Form("bump_fbk_W%i",i+1), Form("Categories on-wafer W%i",i+1), 5,2,7,5,2,7);

  }


  int dumb_counter=0; // Total number of sensors/measurements in the root file
  int non_empty_counter=0; // Sensors with both I and V arrays having size !=0 (the sensor was in fact measured)
  int non_zero_counter=0; // Sensors with both I and V arrays having values !=0 (the sensor was measured and can be biased)

  int total_sensors_counter_fbk=0; //Fraction of sensors whose category can be properly defined
  int good_sensors_counter_fbk=0; //Fraction of GOOD sensors
  int bad_sensors_counter_fbk=0; //Fraction of BAD sensors
  int medium_sensors_counter_fbk=0; //Fraction of MEDIUM sensors

  bool iv_quality; // Boolean reflecting the quality of the IV curve: sensor is BAD whenever iv_quality is set to False

  while( reader_fbk.Next() ){

    dumb_counter++ ;

    if(I_fbk.GetSize()>=start_bd_calculation && V_fbk.GetSize()>=start_bd_calculation){
  
       non_empty_counter++ ; 

      float non_zero_current = 0;
      for(int i=1;  i<I_fbk.GetSize(); i++) non_zero_current += I_fbk.At(i) ;

      float non_zero_voltage = 0;
      for(int i=1;  i<V_fbk.GetSize(); i++) non_zero_voltage += V_fbk.At(i) ;


      if(non_zero_voltage!=0 && non_zero_current!=0){

        non_zero_counter++ ;
    
        iv_quality = true;
        I_fbk_100V = -1000;
        VBD_fbk_u = -1000.;
        VBD_fbk_d = -1000.;
        
        int type_counter = 0;
        int type_counter_I_fbk = 0;
        int type_counter_V_fbk_u = 0;
        int type_counter_V_fbk_d = 0;

        k_fbk_u[0] = 0.;
        k_fbk_u[1] = 0.;
        k_fbk_d[0] = 0.;
        k_fbk_d[1] = 0.;
      
          
          if(invert_polarity) for(int i=0; i<I_fbk.GetSize(); i++) I_fbk.At(i) = -I_fbk.At(i);
      
          
          //////////// I@V_current_level calculation ///////////////////////
          for(int i=0; i<I_fbk.GetSize(); i++){
      
          	if( i!=0 && V_fbk.At(i-1)<V_current_level && V_fbk.At(i)>=V_current_level ){
      
              if( V_fbk.At(i)==V_current_level ) I_fbk_100V = I_fbk.At(i);
              else I_fbk_100V = (((I_fbk.At(i)-I_fbk.At(i-1))/(V_fbk.At(i)-V_fbk.At(i-1)))*(V_current_level-V_fbk.At(i)) + I_fbk.At(i));
              break;
          	}
      
          }
      
          for(int i=1;  i<I_fbk.GetSize(); i++){
      
            if( V_fbk.At(i)<=V_current_monitor && I_fbk.At(i) > I_thr*(1./current_conversion_value) ){
      
              iv_quality=false;
              break;
    
            }
          }
    
          I_fbk_100V = current_conversion_value*I_fbk_100V;
          if( I_fbk_100V<0. ) iv_quality=false;
          
          //////////////// VBD calculation //////////////////////////////
      
          if( I_fbk.At( I_fbk.GetSize()-1 )>I_compliance_minimum*(1./current_conversion_value) ){
      
            
            ///////// VBD "up" (Calculation start from the end of the Voltage array downwards) /////////
            for(int i=(I_fbk.GetSize()-2); i>=start_bd_calculation; i--){
        
              if(V_fbk.At( I_fbk.GetSize()-1 )<V_current_level){
                  
                  VBD_fbk_u = V_fbk.At( I_fbk.GetSize()-1 ) ;
                  break ;
              }
              
              if( I_fbk.At(i)<I_compliance*(1./current_conversion_value) && V_fbk.At(i)>=V_current_level ){
        
                k_fbk_u[0] = ( (I_fbk.At(i)-I_fbk.At(i-1))/(V_fbk.At(i)-V_fbk.At(i-1)) )*(V_fbk.At(i)/I_fbk.At(i)) ;
                k_fbk_u[1] = ( (I_fbk.At(i+1)-I_fbk.At(i))/(V_fbk.At(i+1)-V_fbk.At(i)) )*(V_fbk.At(i)/I_fbk.At(i)) ;
        
                if( k_fbk_u[0]<k_thr && k_fbk_u[1]>=k_thr ){
                  
                  VBD_fbk_u = V_fbk.At(i) ;
                  break ;
                }
              }
            }
        
      
            ///////// VBD "down" (Calculation start from the beginning of the Voltage array upwards)/////////
            for(int i=start_bd_calculation; i<=(I_fbk.GetSize()-2); i++){
        
              if(V_fbk.At( I_fbk.GetSize()-1 )<V_current_level){
                  
                  VBD_fbk_d = V_fbk.At( I_fbk.GetSize()-1 ) ;
                  break ;
              }
        
              if( I_fbk.At(i)<I_compliance*(1./current_conversion_value) && V_fbk.At(i)>=V_current_level ){
                  
                k_fbk_d[0] = ( (I_fbk.At(i)-I_fbk.At(i-1))/(V_fbk.At(i)-V_fbk.At(i-1)) )*(V_fbk.At(i)/I_fbk.At(i)) ;
                k_fbk_d[1] = ( (I_fbk.At(i+1)-I_fbk.At(i))/(V_fbk.At(i+1)-V_fbk.At(i)) )*(V_fbk.At(i)/I_fbk.At(i)) ;
        
                if( k_fbk_d[0]<k_thr && k_fbk_d[1]>=k_thr ){
                  
                  VBD_fbk_d = V_fbk.At(i) ;
                  break ;
                }
              }
            }
          }
      
      
      
          if(VBD_fbk_u!=-1000 && VBD_fbk_d!=-1000){ 
            
            vvbd_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += VBD_fbk_u;
            vcount_vbd_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 1;
      
          }
      
          if(I_fbk_100V>=0){ 
      
            vi_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += I_fbk_100V;
            vcount_i_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 1;
      
          }
      
         
          if( VBD_fbk_u>VBD_expected && VBD_fbk_d!=-1000 &&  VBD_fbk_d<=VBD_expected ) cout<<"!!! SENSOR WITH CURRENT JUMPS ALERT  !!!: "<<"from wafer "<<*wafer_fbk<<" row "<<*row_fbk<<" column "<<*col_fbk<<endl;
          
          
          if( iv_quality && VBD_fbk_u>VBD_expected && VBD_fbk_d!=-1000 ){

            if( VBD_fbk_d>VBD_expected ){
      
              vbump_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 1; 
              vcount_bump_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 1; 
      
            }else{ 
      
              vbump_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 2;
              vcount_bump_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 1;
      
            }
      
          }else if( !iv_quality || (VBD_fbk_d!=-1000 && VBD_fbk_u!=-1000 && VBD_fbk_u<=VBD_expected ) ){
      
            vbump_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 100;
            vcount_bump_fbk[*wafer_fbk-1][*col_fbk-2][*row_fbk-2] += 1;
      
          }else cout<<"Sensor from wafer "<<*wafer_fbk<<" row "<<*row_fbk<<" column "<<*col_fbk<<" has Current within acceptance, but VBD could not be calculated. PLEASE CHECK."<<endl;
      }else cout<<"Sensor from wafer "<<*wafer_fbk<<" row "<<*row_fbk<<" column "<<*col_fbk<<" has Voltage and/or Current always equal to zero. PLEASE CHECK."<<endl;
    }else cout<<"Sensor from wafer "<<*wafer_fbk<<" row "<<*row_fbk<<" column "<<*col_fbk<<Form(": Voltage and/or Current vectors have less than %i elements. PLEASE CHECK.",start_bd_calculation)<<endl;
  }



  for(int i=0; i<18; i++){
    for(int j=0; j<5; j++){
      for(int k=0; k<5; k++){

        if( i!=9 && i!=10 ){
  
          if(vcount_vbd_fbk[i][j][k]!=0) hV_fbk[i]->Fill(j+2,k+2, vvbd_fbk[i][j][k]/vcount_vbd_fbk[i][j][k] );
          if(vcount_i_fbk[i][j][k]!=0) hI_fbk_100V[i]->Fill(j+2,k+2, vi_fbk[i][j][k]/vcount_i_fbk[i][j][k] );
  
          if(vcount_bump_fbk[i][j][k]!=0) total_sensors_counter_fbk++;
          
          if( vbump_fbk[i][j][k]/vcount_bump_fbk[i][j][k]==1 ){
          
            hbump_fbk[i]->Fill(j+2,k+2,1);
            good_sensors_counter_fbk++;
  
          }else if( vbump_fbk[i][j][k]/vcount_bump_fbk[i][j][k]>2 ){
  
            hbump_fbk[i]->Fill(j+2,k+2,3);
            bad_sensors_counter_fbk++;
  
  
          }else if( vbump_fbk[i][j][k]/vcount_bump_fbk[i][j][k]>1 && vbump_fbk[i][j][k]/vcount_bump_fbk[i][j][k]<=2 ){
  
            hbump_fbk[i]->Fill(j+2,k+2,2);
            medium_sensors_counter_fbk++;
  
          }

          if( vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]<10. && vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]>=0. ) current_levels_counter[0]++ ;
          else if( vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]>=10. && vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]<50. ) current_levels_counter[1]++ ;
          else if( vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]>=50. && vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]<100. ) current_levels_counter[2]++ ;
          else if( vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]>=100. && vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]<500. ) current_levels_counter[3]++ ;
          else if( vi_fbk[i][j][k]/vcount_i_fbk[i][j][k]>=500. ) current_levels_counter[4]++ ;
          else if( vcount_i_fbk[i][j][k]==0. && vcount_vbd_fbk[i][j][k]!=0. )  current_levels_counter[5]++ ;

        }
      }
    }
  }



  cout<<"Fraction of GOOD sensors (on-wafer): "<<double(good_sensors_counter_fbk)/double(total_sensors_counter_fbk)<<endl;
  cout<<"Fraction of MEDIUM sensors (on-wafer): "<<double(medium_sensors_counter_fbk)/double(total_sensors_counter_fbk)<<endl;
  cout<<"Fraction of BAD sensors (on-wafer): "<<double(bad_sensors_counter_fbk)/double(total_sensors_counter_fbk)<<endl;
  cout<<"\n";

  /*cout<<"Fraction with I<10 uA: "<<double(current_levels_counter[0])/double(total_sensors_counter_fbk)<<endl;
  cout<<"Fraction with I in 10-50 uA range: "<<double(current_levels_counter[1])/double(total_sensors_counter_fbk)<<endl;
  cout<<"Fraction with I in 50-100 uA range: "<<double(current_levels_counter[2])/double(total_sensors_counter_fbk)<<endl;
  cout<<"Fraction with I in 100-500 uA range: "<<double(current_levels_counter[3])/double(total_sensors_counter_fbk)<<endl;
  cout<<"Fraction with I>500 uA: "<<double(current_levels_counter[4])/double(total_sensors_counter_fbk)<<endl;
  cout<<"Fraction not reaching 100V: "<<double(current_levels_counter[5])/double(total_sensors_counter_fbk)<<endl;
  cout<<"\n";*/
 
      
  TCanvas *cI_100V[18];
  TCanvas *cV[18];
  TCanvas *cbump[18];

  for(int i=0; i<18; i++){
   
    if( i!=9 && i!=10 ){

      //hI_fbk_100V[i]->GetZaxis()->SetRangeUser( 0.1, hI_fbk_100V[i]->GetMaximum() ); //Alternative colored axis range
      hI_fbk_100V[i]->GetZaxis()->SetRangeUser( low_iv_range, high_iv_range );
      hI_fbk_100V[i]->GetZaxis()->SetTitle("[uA]");
      hI_fbk_100V[i]->GetXaxis()->SetTitle("column"); 
      hI_fbk_100V[i]->GetYaxis()->SetTitle("row"); 
      hI_fbk_100V[i]->SetMarkerSize(3.);

  
      //hV_fbk[i]->GetZaxis()->SetRangeUser( hV_fbk[i]->GetMinimum(), hV_fbk[i]->GetMaximum() ); //Alternative colored axis range
      hV_fbk[i]->GetZaxis()->SetRangeUser( low_vbd_range, high_vbd_range ); 
      hV_fbk[i]->GetZaxis()->SetTitle("[V]");  
      hV_fbk[i]->GetXaxis()->SetTitle("column"); 
      hV_fbk[i]->GetYaxis()->SetTitle("row"); 
      hV_fbk[i]->SetMarkerSize(3.);


      hbump_fbk[i]->GetZaxis()->SetRangeUser( 1, 3 );
      hbump_fbk[i]->GetZaxis()->SetTitle("Category");
      hbump_fbk[i]->GetXaxis()->SetTitle("column");
      hbump_fbk[i]->GetYaxis()->SetTitle("row");


      if(bcurrent){

        cI_100V[i]=new TCanvas(Form("c_I_%iV_W%i",int(V_current_level),i+1), Form("c I@%iV_W%i",int(V_current_level),i+1), 1000,1000);
        cI_100V[i]->SetRightMargin(0.15);
        cI_100V[i]->cd();
        hI_fbk_100V[i]->Draw("textcolz");
        gPad->SetGrid(1,1);
        gPad->SetLogz(1);
        gPad->Update();
        cI_100V[i]->Update();
        hI_fbk_100V[i]->GetXaxis()->SetNdivisions(4);
        hI_fbk_100V[i]->GetYaxis()->SetNdivisions(3);
  
        if(save) cI_100V[i]->SaveAs( Form("pics/I_%iV_W%i_FINAL.png",int(V_current_level),i+1) );

      }


      if(bvoltage){

        gStyle->SetPaintTextFormat("1.1f");
        cV[i]=new TCanvas(Form("c_V_W%i",i+1), Form("c VBD_W%i",i+1), 1000,1000);
        cV[i]->SetRightMargin(0.15);
        cV[i]->cd();
        hV_fbk[i]->Draw("textcolz");
        gPad->SetGrid(1,1);
        gPad->Update();
        cV[i]->Update();
        hV_fbk[i]->GetXaxis()->SetNdivisions(4);
        hV_fbk[i]->GetYaxis()->SetNdivisions(3);
  
        if(save) cV[i]->SaveAs( Form("pics/VBD_W%i_FINAL.png",i+1) );

      }


      if(bcategory){
      
        cbump[i]=new TCanvas(Form("c_bump_W%i",i+1), Form("c bump_W%i",i+1), 1000,1000);
        cbump[i]->SetRightMargin(0.15);
        cbump[i]->cd();
        hbump_fbk[i]->Draw("colz");
        gPad->SetGrid(1,1);
        gPad->Update();
        cbump[i]->Update();
        hbump_fbk[i]->GetXaxis()->SetNdivisions(4);
        hbump_fbk[i]->GetYaxis()->SetNdivisions(3);
       
        if(save) cbump[i]->SaveAs( Form("pics/categories_W%i_FINAL.png",i+1) );

        

      }

    }

  }

}
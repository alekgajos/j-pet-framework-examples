/**
 *  @copyright Copyright 2019 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file TimeWindowCreator.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include "TimeWindowCreator.h"
#include <EventIII.h>
#include <tuple>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>

using namespace jpet_options_tools;
using namespace std;

TimeWindowCreator::TimeWindowCreator(const char* name): JPetUserTask(name) {}

TimeWindowCreator::~TimeWindowCreator() {}

double fixTime(double dt){
  int offset = round(dt / 2.7027);
  return dt - offset * 2.7027;
}

bool TimeWindowCreator::init()
{
  INFO("TimeSlot Creation Started");
  fOutputEvents = new JPetTimeWindow("JPetSigCh");
  
  // Reading values from the user options if available
  // Min allowed signal time
  if (isOptionSet(fParams.getOptions(), kMinTimeParamKey)) {
    fMinTime = getOptionAsFloat(fParams.getOptions(), kMinTimeParamKey);
  } else {
    WARNING(
      Form("No value of the %s parameter provided by the user. Using default value of %lf.",
        kMinTimeParamKey.c_str(), fMinTime
      )
    );
  }
  // Max allowed signal time
  if (isOptionSet(fParams.getOptions(), kMaxTimeParamKey)) {
    fMaxTime = getOptionAsFloat(fParams.getOptions(), kMaxTimeParamKey);
  } else {
    WARNING(
      Form("No value of the %s parameter provided by the user. Using default value of %lf.",
        kMaxTimeParamKey.c_str(), fMaxTime
      )
    );
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)) {
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }

  /************************************************************************/
  /* Read calibration of inter-PM time offsets                            */
  /************************************************************************/
  ifstream f("/home/alek/FTAB/new_firmware/configs/time_calibration.txt", ios::in);

  int side, scin, pm;
  double offset;
  while(!f.eof()){
    f >> side >> scin >> pm >> offset;
    fTimeOffsets[{side-1,scin}][pm-1] = offset;
  }

  f.close();
  
  // find earliest PM in each matrix
  for(auto& side_scin: fTimeOffsets){
    double min = 1000.;
    int earliest = 5;
    for(int i=0; i<4; ++i){
      if( side_scin.second[i] < min){
        earliest = i;
        min = side_scin.second[i];
      }
    }
    fEarliestPM[side_scin.first] = earliest;
  }

  /************************************************************************/
  /* My histos                                                            */
  /************************************************************************/
  getStatistics().createHistogram(
                                  new TH1F("sig_tresholds", "Number of leading edge thresholds in signal",
                                           5, -0.5, 4.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("sig_tots", "Number of good tots in signal",
                                           5, -0.5, 4.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("sig_tot", "Signal tot",
                                           200, 0., 200.)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("sigs_per_matrix", "Number of signals per matrix",
                                           5, -0.5, 4.5)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("inter_pm_tdiffs", "Time differences between PM-s",
                                           200, -20., 20.)
                                  );

  getStatistics().createHistogram(
                                  new TH1F("hit_tot", "Hit tot - Scin 13",
                                           200, 0., 400.)
                                  );
  
  getStatistics().createHistogram(
                                  new TH1F("h_pmt_no", "Observed PM", 4, 0.5, 4.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_strip_no", "Observed Scin", 13, 0.5, 13.5)
                                  );
  getStatistics().createHistogram(
                                  new TH1F("h_thr_tdiff", "Time difference between A and B thresholds;#Delta t_{AB} [ns]", 1000, -10., 10.)
                                  );

  getStatistics().createHistogram(new TH1F("Inter-module Dt", "Inter-module Dt;#Delta t [ns]",
                                           100000, -50000, 50000));

  getStatistics().createHistogram(new TH1F("dt_mean", "Scin 13 Dt using mean time;#Delta t [ns]",
                                           1000, -20, 20));

  getStatistics().createHistogram(new TH1F("dt_earliest", "Scin 13 Dt using earliest time;#Delta t [ns]",
                                           1000, -20, 20));
  
  
  // Control histograms
  if (fSaveControlHistos) { initialiseHistograms(); }
  return true;
}

bool TimeWindowCreator::exec()
{
  if (auto event = dynamic_cast<EventIII* const> (fEvent)) {
    int kTDCChannels = event->GetTotalNTDCChannels();
    if (fSaveControlHistos){
      getStatistics().getHisto1D("sig_ch_per_time_slot")->Fill(kTDCChannels);
    }
    // Loop over all TDC channels in file
    auto tdcChannels = event->GetTDCChannelsArray();

    vector<JPetSigCh> mySigChs;

    for (int i = 0; i < kTDCChannels; ++i) {
      auto tdcChannel = dynamic_cast<TDCChannel* const> (tdcChannels->At(i));
      auto channelNumber = tdcChannel->GetChannel();
      
      //      Skip trigger signals - every 105th
      if (channelNumber % 105 == 104) continue;

      // Check if channel exists in database from loaded local file
      if (getParamBank().getChannels().count(channelNumber) == 0) {
        if (fSaveControlHistos){
          getStatistics().getHisto1D("channel_missing")->Fill(channelNumber);
        }
        WARNING(
          Form("DAQ Channel %d appears in data but does not exist in the detector setup.", channelNumber)
        );
        continue;
      }

      // Get channel for corresponding number
      auto& channel = getParamBank().getChannel(channelNumber);

      // Loop over all entries on leading edge in current TDCChannel and create SigCh
      for (int j = 0; j < tdcChannel->GetLeadHitsNum(); j++) {
        auto leadTime = tdcChannel->GetLeadTime(j);
        if (leadTime > fMaxTime || leadTime < fMinTime ) {
	  continue; }
		  
        auto leadSigCh = generateSigCh(
				       leadTime, channel, JPetSigCh::Leading
        );

	mySigChs.push_back(leadSigCh);
        if (fSaveControlHistos){
          fillChannelHistos(channel, JPetSigCh::Leading);
        }
      }

      // Loop over all entries on trailing edge in current TDCChannel and create SigCh
      for (int j = 0; j < tdcChannel->GetTrailHitsNum(); j++) {
        auto trailTime = tdcChannel->GetTrailTime(j);
        if (trailTime > fMaxTime || trailTime < fMinTime ) { continue; }

        auto trailSigCh = generateSigCh(
          trailTime, channel, JPetSigCh::Trailing
        );

	mySigChs.push_back(trailSigCh);
        if (fSaveControlHistos){
          fillChannelHistos(channel, JPetSigCh::Trailing);
        }
      }
      
    }


    /********************************************************************/
    /* My study                                                         */
    /********************************************************************/
    //
    using times = vector<double>;
    using timesByEdge = array<times, 2>;
    using timesByThreshold = array<timesByEdge, 3>;
    using timesByPM = array<timesByThreshold,5>;
    using timesByScin = array<timesByPM,28>;
    using timesBySide = array<timesByScin,2>;

    timesBySide all_sigchs;
    	
    for(auto& sc: mySigChs){
      
      int scin = sc.getChannel().getPM().getScin().getID() - 200;
      JPetPM::Side side = sc.getChannel().getPM().getSide();

      int pm = sc.getChannel().getPM().getMatrixPosition();

      int thr = sc.getChannel().getThresholdNumber();
      JPetSigCh::EdgeType type = sc.getType();
      double time = sc.getTime() / 1000.;

      getStatistics().getHisto1D("h_pmt_no")->Fill(pm);
      getStatistics().getHisto1D("h_strip_no")->Fill(scin);

      int side_no = 0; // 0 - Side A, 1 - Side B
      if(side == JPetPM::SideA){
        side_no = 0;
      }else{
        side_no = 1;
      }

      int edge = 0; // 0 - Leading, 1 - Trailing
      if(type==JPetSigCh::Leading){
        edge = 0;
      }else{
        edge = 1;
      }

      all_sigchs[side_no][scin][pm][thr][edge].push_back( time );  
    }
    
    /**************************************************************************/
    /* Actual clustering                                                      */
    /**************************************************************************/
    using signalsByScin = array<vector<Signal>,14>;
    using signalsBySide = array<signalsByScin,2>;    

    signalsBySide Signals;
    
    for(int side=0;side<=1;++side){
      for(int scin=1;scin<=13;scin+=12){

        for(int pm = 1; pm <= 4; ++pm){
         
	  times& thr_a_leads = all_sigchs[side][scin][pm][1][0]; 
	  times& thr_b_leads = all_sigchs[side][scin][pm][2][0]; 

	  times& thr_a_trails = all_sigchs[side][scin][pm][1][1]; 
	  times& thr_b_trails = all_sigchs[side][scin][pm][2][1]; 

          /**********************************************************************/
          /* Correct hit times for inter-PM offsets                             */
          /**********************************************************************/
          applyTimeCalibration(thr_a_leads, side, scin, pm, 1);
          applyTimeCalibration(thr_b_leads, side, scin, pm, 2);
          applyTimeCalibration(thr_a_trails, side, scin, pm, 1);
          applyTimeCalibration(thr_b_trails, side, scin, pm, 2);
        }          

        // find leading edge candidates          
        for(int pm1=1; pm1<=4; ++pm1){

          if( all_sigchs[side][scin][pm1][1][0].size() != 1 ){
            continue;
          }

          double t1 = all_sigchs[side][scin][pm1][1][0].front();
          
          for(int pm2=pm1+1; pm2<=4; ++pm2){

            if( all_sigchs[side][scin][pm2][1][0].size() != 1 ){
              continue;
            }

            double t2 = all_sigchs[side][scin][pm2][1][0].front();

            if (scin == 13) {
              getStatistics().getHisto1D("inter_pm_tdiffs")->Fill(t1 - t2);
            }

            if( fabs(t1 -t2) < 5.0 ){

              // we have a new signal candidate or add to an exisiting one
              if( Signals[side][scin].empty() ){
                Signal sig;
                sig.leads[0][pm1-1] = t1;
                sig.leads[0][pm2-1] = t2;
                sig.n_thresholds = 2;
                Signals[side][scin].push_back(sig);
              }else{
                for(Signal& sig: Signals[side][scin]){

                  for(int p=0;p<4;++p){
                    if( sig.leads[0][p] >= 0. ){

                      if( p+1 != pm1 || p+1 != pm2 ){

                        if( fabs( sig.leads[0][p] - t1 ) < 5.0 &&
                            fabs( sig.leads[0][p] - t1 ) < 5.0 ){

                          // extend the exisiting signal
                          sig.leads[0][pm1-1] = t1;
                          sig.leads[0][pm2-1] = t2;
                                                    
                        }
                        
                      }
                      
                    }
                  }
                }
                
              }
            }

              
          }
            
        }

        // completely assemble the found signals
        for(Signal& sig: Signals[side][scin]){

          // count fired thresholds and calc average time
          sig.t_mean = 0.;
          sig.n_thresholds = 0;
          for(int pm=0;pm<4;++pm){
            if(sig.leads[0][pm] >= 0.){
              sig.n_thresholds++;
              sig.t_mean += sig.leads[0][pm];
            }
          }
          sig.t_mean /= sig.n_thresholds;

          // fill time on earliest element
          sig.t_earliest = sig.leads[0][fEarliestPM[{side, scin}]];

          
          // find TOT
          double tot = 0.;
          int n_tots = 0;
          for(int pm=0;pm<4;++pm){
            if(sig.leads[0][pm] >= 0.){
              if(all_sigchs[side][scin][pm+1][1][1].size() == 1){
                tot += all_sigchs[side][scin][pm+1][1][1].front() - sig.leads[0][pm];
                n_tots++;
              }
            }
          }
          sig.tot = tot / n_tots;
          sig.n_tots = n_tots;

          // fill histos
          getStatistics().getHisto1D("sig_tresholds")->Fill(sig.n_thresholds);
          getStatistics().getHisto1D("sig_tots")->Fill(sig.n_tots);
          getStatistics().getHisto1D("sig_tot")->Fill(sig.tot);
        }

        getStatistics().getHisto1D("sigs_per_matrix")->Fill(Signals[side][scin].size());
        
      }      
    }
  
    // require a tagging signal on Scin 1 side 0 (side 1 does not work!)
    if(Signals[0][1].empty()){
      // return true;
    }

    if(Signals[0][13].size() == 1 && Signals[1][13].size() == 1){

      Signal sigA = Signals[0][13].front();
      Signal sigB = Signals[1][13].front();
      
      double t_hit = 0.5 * (sigA.t_mean + sigB.t_mean);

      // double inter_module_dt = t_hit - Signals[0][1].front().t_mean;
      // getStatistics().getHisto1D("Inter-module Dt")->Fill(inter_module_dt);

      // if(fabs(inter_module_dt) < 10.0) { // a good coincidence!

        double dt_mean = sigA.t_mean - sigB.t_mean;
        double dt_earliest = sigA.t_earliest - sigB.t_earliest;

        getStatistics().getHisto1D("dt_mean")->Fill(dt_mean);
        getStatistics().getHisto1D("dt_earliest")->Fill(dt_earliest);
        getStatistics().getHisto1D("hit_tot")->Fill(sigA.tot + sigB.tot);
      // }
      
    }

    
    /********************************************************************/
    /* End of my study                                                  */
    /********************************************************************/
    

    fCurrEventNumber++;
  } else { return false; }
  return true;
}

/**
* Sets up Signal Channel fields
*/
void TimeWindowCreator::fillChannelHistos(
  const JPetChannel& channel, JPetSigCh::EdgeType edge
) {
  getStatistics().getHisto1D("channel_occ")->Fill(channel.getID());
  getStatistics().getHisto1D("channel_thrnum")->Fill(channel.getThresholdNumber());
  getStatistics().getHisto1D("matrix_occ")->Fill(channel.getPM().getMatrixPosition());
  getStatistics().getHisto1D("pm_occ")->Fill(channel.getPM().getID());
  getStatistics().getHisto1D("scin_occ")->Fill(channel.getPM().getScin().getID());
  getStatistics().getHisto1D("slot_occ")->Fill(channel.getPM().getScin().getID());

  if(edge == JPetSigCh::Leading){
    getStatistics().getHisto1D("channel_occ_leads")->Fill(channel.getID());
    getStatistics().getHisto1D("channel_thrnum_leads")->Fill(channel.getThresholdNumber());
    getStatistics().getHisto1D("matrix_occ_leads")->Fill(channel.getPM().getMatrixPosition());
    getStatistics().getHisto1D("pm_occ_leads")->Fill(channel.getPM().getID());
    getStatistics().getHisto1D("scin_occ_leads")->Fill(channel.getPM().getScin().getID());
    getStatistics().getHisto1D("slot_occ_leads")->Fill(channel.getPM().getScin().getID());
  } else if(edge == JPetSigCh::Trailing){
    getStatistics().getHisto1D("channel_occ_trails")->Fill(channel.getID());
    getStatistics().getHisto1D("channel_thrnum_trails")->Fill(channel.getThresholdNumber());
    getStatistics().getHisto1D("matrix_occ_trails")->Fill(channel.getPM().getMatrixPosition());
    getStatistics().getHisto1D("pm_occ_trails")->Fill(channel.getPM().getID());
    getStatistics().getHisto1D("scin_occ_trails")->Fill(channel.getPM().getScin().getID());
    getStatistics().getHisto1D("slot_occ_trails")->Fill(channel.getPM().getScin().getID());
  }

  if(channel.getPM().getSide() == JPetPM::SideA){
    getStatistics().getHisto1D("pm_occ_sides")->Fill(1);
  }else if(channel.getPM().getSide() == JPetPM::SideB){
    getStatistics().getHisto1D("pm_occ_sides")->Fill(2);
  }
}

/**
* Sets up Signal Channel fields
*/
JPetSigCh TimeWindowCreator::generateSigCh(
  double time, const JPetChannel& channel, JPetSigCh::EdgeType edge
) {
  JPetSigCh sigCh;
  sigCh.setTime(1000.*time);
  sigCh.setType(edge);
  sigCh.setChannel(channel);
  sigCh.setRecoFlag(JPetSigCh::Good);
  return sigCh;
}

bool TimeWindowCreator::terminate()
{
  INFO("TimeSlot Creation Ended");
  return true;
}

void TimeWindowCreator::saveSigChs(const vector<JPetSigCh>& sigChVec)
{
  //  for (auto & sigCh : sigChVec) { fOutputEvents->add<JPetSigCh>(sigCh); }
}

void TimeWindowCreator::initialiseHistograms(){



  getStatistics().createHistogram(
    new TH1F("sig_ch_per_time_slot", "Signal Channels Per Time Slot", 50, -0.5, 50.5)
  );
  getStatistics().getHisto1D("sig_ch_per_time_slot")
  ->GetXaxis()->SetTitle("Signal Channels in Time Slot");
  getStatistics().getHisto1D("sig_ch_per_time_slot")
  ->GetYaxis()->SetTitle("Number of Time Slots");

  // Channels
  getStatistics().createHistogram(
    new TH1F("channel_occ", "Channels occupation", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_occ")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_occ")->GetYaxis()->SetTitle("Number of SigCh");

  // 2204 2309
  // high occ: 2198, 2199, 2204, 2205, 2296, 2303, 2304, 2305, 2308, 2309
  getStatistics().createHistogram(
    new TH1F("channel_missing", "Channels missing in configuration", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_missing")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_missing")->GetYaxis()->SetTitle("Counts");

  getStatistics().createHistogram(
    new TH1F("channel_thrnum", "Channels threshold numbers", 4, 0.5, 4.5)
  );
  getStatistics().getHisto1D("channel_thrnum")->GetXaxis()->SetTitle("Channel Threshold Number");
  getStatistics().getHisto1D("channel_thrnum")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_thrnum_leads", "Channels threshold numbers LEADS", 4, 0.5, 4.5)
  );
  getStatistics().getHisto1D("channel_thrnum_leads")->GetXaxis()->SetTitle("Channel Threshold Number");
  getStatistics().getHisto1D("channel_thrnum_leads")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_thrnum_trails", "Channels threshold numbers TRAILS", 4, 0.5, 4.5)
  );
  getStatistics().getHisto1D("channel_thrnum_trails")->GetXaxis()->SetTitle("Channel Threshold Number");
  getStatistics().getHisto1D("channel_thrnum_trails")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_occ_leads", "Channels occupation - Leading channels", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_occ_leads")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_occ_leads")->GetYaxis()->SetTitle("Number of Lading SigCh");

  getStatistics().createHistogram(
    new TH1F("channel_occ_trails", "Channels occupation - Trailing channels", 211, 2099.5, 2310.5)
  );
  getStatistics().getHisto1D("channel_occ_trails")->GetXaxis()->SetTitle("Channel ID");
  getStatistics().getHisto1D("channel_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");

  // SiPMs
  getStatistics().createHistogram(
    new TH1F("pm_occ", "PMs occupation", 111, 399.5, 510.5)
  );
  getStatistics().getHisto1D("pm_occ")->GetXaxis()->SetTitle("PM ID");
  getStatistics().getHisto1D("pm_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("pm_occ_leads", "PMs occupation LEADS", 111, 399.5, 510.5)
  );
  getStatistics().getHisto1D("pm_occ_leads")->GetXaxis()->SetTitle("PM ID");
  getStatistics().getHisto1D("pm_occ_leads")->GetYaxis()->SetTitle("Number of Leading SigCh");

  getStatistics().createHistogram(
    new TH1F("pm_occ_trails", "PMs occupation TRAILS", 111, 399.5, 510.5)
  );
  getStatistics().getHisto1D("pm_occ_trails")->GetXaxis()->SetTitle("PM ID");
  getStatistics().getHisto1D("pm_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");

  getStatistics().createHistogram(
    new TH1F("pm_occ_sides", "PMs occupation of sides", 2, 0.5, 2.5)
  );
  getStatistics().getHisto1D("pm_occ_sides")->GetXaxis()->SetBinLabel(1, "SIDE A");
  getStatistics().getHisto1D("pm_occ_sides")->GetXaxis()->SetBinLabel(2, "SIDE B");
  getStatistics().getHisto1D("pm_occ_sides")->GetYaxis()->SetTitle("Number of SigCh");

  // Matrix position
  getStatistics().createHistogram(
    new TH1F("matrix_occ", "Position in matrix in PMs occupation", 5, -0.5, 4.5)
  );
  getStatistics().getHisto1D("matrix_occ")->GetXaxis()->SetTitle("Matrix position");
  getStatistics().getHisto1D("matrix_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("matrix_occ_leads", "Position in matrix in PMs occupation LEADS", 5, -0.5, 4.5)
  );
  getStatistics().getHisto1D("matrix_occ_leads")->GetXaxis()->SetTitle("Matrix position");
  getStatistics().getHisto1D("matrix_occ_leads")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("matrix_occ_trails", "Position in matrix in PMs occupation TRAILS", 5, -0.5, 4.5)
  );
  getStatistics().getHisto1D("matrix_occ_trails")->GetXaxis()->SetTitle("Matrix position");
  getStatistics().getHisto1D("matrix_occ_trails")->GetYaxis()->SetTitle("Number of SigCh");

  // Scins
  getStatistics().createHistogram(
    new TH1F("scin_occ", "Scins occupation", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("scin_occ")->GetXaxis()->SetTitle("SCIN ID");
  getStatistics().getHisto1D("scin_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("scin_occ_leads", "Scins occupation LEADS",16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("scin_occ_leads")->GetXaxis()->SetTitle("SCIN ID");
  getStatistics().getHisto1D("scin_occ_leads")->GetYaxis()->SetTitle("Number of Leading SigCh");

  getStatistics().createHistogram(
    new TH1F("scin_occ_trails", "Scins occupation TRAILS", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("scin_occ_trails")->GetXaxis()->SetTitle("SCIN ID");
  getStatistics().getHisto1D("scin_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");

  // Slots
  getStatistics().createHistogram(
    new TH1F("slot_occ", "Slots occupation", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("slot_occ")->GetXaxis()->SetTitle("SLOT ID");
  getStatistics().getHisto1D("slot_occ")->GetYaxis()->SetTitle("Number of SigCh");

  getStatistics().createHistogram(
    new TH1F("slot_occ_leads", "Slots occupation LEADS", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("slot_occ_leads")->GetXaxis()->SetTitle("SLOT ID");
  getStatistics().getHisto1D("slot_occ_leads")->GetYaxis()->SetTitle("Number of Leading SigCh");

  getStatistics().createHistogram(
    new TH1F("slot_occ_trails", "Slots occupation TRAILS", 16, 199.5, 215.5)
  );
  getStatistics().getHisto1D("slot_occ_trails")->GetXaxis()->SetTitle("SLOT ID");
  getStatistics().getHisto1D("slot_occ_trails")->GetYaxis()->SetTitle("Number of Trailing SigCh");
}

void TimeWindowCreator::applyTimeCalibration(std::vector<double>& times, int side, int scin, int pm, int thr){

  double offset = fTimeOffsets[{side,scin}][pm-1];

  for(int k=0; k<times.size();++k){
    times[k] += offset;
  }

  std::sort(times.begin(), times.end());
  
}

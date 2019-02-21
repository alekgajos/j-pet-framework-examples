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
 *  @file OPSAnalyzer.cpp
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <JPetWriter/JPetWriter.h>
#include "OPSAnalyzer.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include "JPetOpsEvent.h"
#include <TMath.h>
#include <functional>
using namespace jpet_options_tools;

using namespace std;

OPSAnalyzer::OPSAnalyzer(const char* name): JPetUserTask(name) {}

bool OPSAnalyzer::init()
{

  INFO("Analysis of previously reocnstructed o-Ps->3g decays started.");

  fOutputEvents = new JPetTimeWindow("JPetOpsEvent");

  getStatistics().createHistogram(
				  new TH1F("anh_events_in_tw",
					   "Number of annihilation candidate hits in time window", 10, -0.5, 9.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("dex_events_in_tw",
					   "Number of deexcitation candidate hits time window", 10, -0.5, 9.5)
				  );

  getStatistics().createHistogram(
				  new TH2F("dex_vs_anh_events_in_tw",
					   "Number of deexcitation vs annihilation events"
					   " in a time window; annihilation; deexcitation",
					   10, -0.5, 9.5, 10, -0.5, 9.5)
				  );
  
  // create histograms for annihilation position
  getStatistics().createHistogram(new TH2F("anh_XY",
					   "transverse position of the o-Ps->3g decay point;"
					   "X [cm]; Y [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );

  getStatistics().createHistogram(new TH2F("anh_XY_nocenter",
					   "transverse position of the o-Ps->3g decay point;"
					   "X [cm]; Y [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );
  
  getStatistics().createHistogram(new TH2F("anh_XZ",
					   "position of the o-Ps->3g decay point in XZ;"
					   "Z [cm]; X [cm]",
					   100, -50., 50.,
					   100, -50., 50.
					   )
				  );  
  
  getStatistics().createHistogram(
				  new TH1F("t_dex_anh",
					   "Time between deexcitation and annihilation"
					   ";#Delta t [ns]",
					   20000, -1000., 1000.)
				  );  

  // create histogram for true event angles
  getStatistics().createHistogram(
				  new TH2F("3_hit_angles",
					   "3 Hit true angles difference",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");

  getStatistics().createHistogram(
				  new TH2F("3_hit_angles_2g_band",
					   "3 Hit true angles difference for evts from 2g band in the theta plot",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles_2g_band")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles_2g_band")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");

    getStatistics().createHistogram(
				  new TH2F("3_hit_angles_right_part",
					   "3 Hit angles difference for evts right of the 2g  band in the theta plot",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles_right_part")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles_right_part")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");
  
  return true;
}

bool OPSAnalyzer::exec()
{
  
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n = timeWindow->getNumberOfEvents();

    std::vector<JPetOpsEvent> events = makeOPSEvents(*timeWindow);
    
    /*    
    for(uint i=0;i<n;++i){
      const JPetOpsEvent & event = timeWindow->getEvent<JPetOpsEvent>(i);

      getStatistics().getHisto1D("nhits")->Fill(event.getHits().size());
      auto hits = event.getHits();

      TVector3 anh_point = event.getAnnihilationPoint();
      
      getStatistics().getHisto2D("decay point XY")->Fill(anh_point.X(),
							 anh_point.Y());

      getStatistics().getHisto2D("decay point XZ")->Fill(anh_point.Z(),
							 anh_point.X());



      //
      // calculate true event angles
      //
      
      // remove deexcitation photon candidate
      hits.erase(
		 remove_if(hits.begin(), hits.end(),
			   [](const JPetHit& hit){
			     if( hit.getQualityOfEnergy() > 0.5 ){
			       return true;
			     }
			     return false;
			   }),
		 hits.end());

      assert(hits.size() == 3);

      TVector3 k[3];
      for(int i=0;i<3;++i){
	k[i] = hits.at(i).getPos() - event.getAnnihilationPoint();
      }
      std::vector<double> angles(3, -1.);
      angles[0] = 180.*k[0].Angle(k[1]) / TMath::Pi();
      angles[1] = 180.*k[1].Angle(k[2]) / TMath::Pi();
      angles[2] = 180.*k[0].Angle(k[2]) / TMath::Pi();
      std::sort(angles.begin(), angles.end());
      assert(angles[0] <= angles[1] <= angles[2]);
      
      double sum_2_least = angles[0] + angles[1];
      double diff_2_least = angles[1] - angles[0];
      if(sum_2_least < 180.){
	angles[2] = 360. - angles[2];
      }
      
      getStatistics().getHisto2D("3_hit_angles")->Fill(sum_2_least, diff_2_least);

      // for comparison, calculate sum of detector thetas in the XY plane
      double theta_sum = calcThetaSum(hits);
      if( fabs(theta_sum-180.) < 34. ){
	getStatistics().getHisto2D("3_hit_angles_2g_band")->Fill(sum_2_least, diff_2_least);
      }
      if( theta_sum > 215. ){
	getStatistics().getHisto2D("3_hit_angles_right_part")->Fill(sum_2_least, diff_2_least);
      }
      
    } // end loop over o-Ps events in a time window
*/      
  } else {
    return false;
  }
  return true;
}

double OPSAnalyzer::calcThetaSum(const std::vector<JPetHit>& hits){

  JPetHit firstHit = hits.at(0);
  JPetHit secondHit = hits.at(1);
  JPetHit thirdHit = hits.at(2);

  std::vector<double> angles;
  angles.push_back(firstHit.getBarrelSlot().getTheta());
  angles.push_back(secondHit.getBarrelSlot().getTheta());
  angles.push_back(thirdHit.getBarrelSlot().getTheta());
  std::sort( angles.begin(), angles.begin() +3 );
  float theta_1_2 = angles[1] - angles[0];
  float theta_2_3 = angles[2] - angles[1];
  float theta_3_1 = 360 - theta_1_2 - theta_2_3;
  angles.clear();
  angles.push_back(theta_1_2);
  angles.push_back(theta_2_3);
  angles.push_back(theta_3_1);
  std::sort( angles.begin(), angles.begin() +3 );

  return angles[0]+angles[1];
}

bool OPSAnalyzer::terminate()
{
  INFO("Analysis of o-Ps->3g decays done.");
  return true;
}

std::vector<JPetOpsEvent> OPSAnalyzer::makeOPSEvents(const JPetTimeWindow& time_window){

  vector<JPetOpsEvent> newEventVec;

  int n_events = time_window.getNumberOfEvents();

  int n_prompt = 0;
  int n_anh = 0;

  std::vector<std::reference_wrapper<const JPetHit>> prompt_hits;
  std::vector<std::reference_wrapper<const JPetOpsEvent>> anh_events;
  
  for(int i=0;i<n_events;++i){
    const JPetOpsEvent & event = time_window.getEvent<JPetOpsEvent>(i);
    
    if( event.isTypeOf(JPetEventType::kPrompt)){
      for(auto & hit: event.getHits()){
        if( hit.getQualityOfEnergy() > 1.8 && hit.getQualityOfEnergy() < 2.2){
          prompt_hits.push_back(std::ref(hit));
        }
      }
      n_prompt++;
    }
    if(event.isTypeOf(JPetEventType::k3Gamma)){
      anh_events.push_back(std::ref(event));
      n_anh++;
    }
  }
  getStatistics().getHisto1D("anh_events_in_tw")->Fill(n_anh);
  getStatistics().getHisto1D("dex_events_in_tw")->Fill(n_prompt);
  getStatistics().getHisto2D("dex_vs_anh_events_in_tw")->Fill(n_anh, n_prompt);

  // temporary
  // only consider time windows with exactly 1 prompt and 1 annihilation
  if( n_anh==1 && n_prompt==1 && prompt_hits.size()==1 ){

    TVector3 vertex = anh_events.front().get().getAnnihilationPoint();
    double t_prompt_corr = prompt_hits.front().get().getTime() - 1000.*(prompt_hits.front().get().getPos() - vertex).Mag() / kSpeedOfLight;
    double dt = anh_events.front().get().getAnnihilationTime() - t_prompt_corr;
    getStatistics().getHisto1D("t_dex_anh")->Fill(dt/1000.);
    
    JPetOpsEvent new_event = anh_events.front().get();
    new_event.setHasPrompt(true);
    new_event.setLifeTime(dt);
    newEventVec.push_back(new_event);
  }
  
  return newEventVec;
}
  
  
  /*  
  for(int entry=0; entry<nevents; ++entry){
    const JPetEvent& event = dynamic_cast<const JPetEvent&>(events[entry]);
    
    const auto & hits = event.getHits();
    
    // case 1: dex and annihilation in the same TW
    if( event.isTypeOf(JPetEventType::kPrompt) && 
        event.isTypeOf(JPetEventType::k3Gamma) ){
      
    }
  
    return newEventVec;
  */



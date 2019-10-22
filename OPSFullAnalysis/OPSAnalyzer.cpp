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
#include "OPSAnalysisTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include "JPetOpsEvent.h"
#include <TMath.h>
#include <functional>
using namespace jpet_options_tools;
using namespace ops_analysis_tools;

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

  getStatistics().createHistogram(new TH2F("anh_XY_nocenter_lt_cut",
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

  // technical reconstruction checks
  getStatistics().createHistogram(
				  new TH1F("E_sum", "Sum of 3 photon energies", 100, -5.5, 1094.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("wrong_kinematics",
					   "Did the event have unphysical kinematics?",
					   2, -0.5, 1.5)
				  );

  
  getStatistics().createHistogram(
				  new TH1F("P_sum", "Square of total 3 photon momentum",
					   100, -3., 597.)
				  );


  // 2d plots of energies and angles
  getStatistics().createHistogram(
				  new TH2F("E2_E1",
					   "E2 vs E1;E1 [keV];E2 [keV]",
					   600, -0.5, 599.5,
					   600, -0.5, 599.5)
				  );

  getStatistics().createHistogram(
				  new TH2F("E3_E1",
					   "E3 vs E1;E1 [keV];E3 [keV]",
					   600, -0.5, 599.5,
					   600, -0.5, 599.5)
				  );

  getStatistics().createHistogram(
				  new TH2F("E3_E2",
					   "E3 vs E2;E2 [keV];E3 [keV]",
					   600, -0.5, 599.5,
					   600, -0.5, 599.5)
				  );

  getStatistics().createHistogram(
				  new TH2F("th2_th1",
					   "th2 vs th1;th1 [deg];th2 [deg]",
					   200, -0.5, 189.5,
					   200, -0.5, 189.5)
				  );

  getStatistics().createHistogram(
				  new TH2F("th3_th1",
					   "th3 vs th1;th1 [deg];th3 [deg]",
					   200, -0.5, 189.5,
					   200, -0.5, 189.5)
				  );

  getStatistics().createHistogram(
				  new TH2F("th3_th2",
					   "th3 vs th2;th2 [deg];th3 [deg]",
					   200, -0.5, 189.5,
					   200, -0.5, 189.5)
				  );

  
  
  return true;
}

bool OPSAnalyzer::exec()
{
  
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n = timeWindow->getNumberOfEvents();

    std::vector<JPetOpsEvent> events_tmp = filterAnnihilationCandidates(*timeWindow);
    std::vector<JPetOpsEvent> events = makeOPSEvents(events_tmp);

    for(auto& ev: events){
      fOutputEvents->add<JPetOpsEvent>(ev);        
    }
    
      
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

std::vector<JPetOpsEvent> OPSAnalyzer::makeOPSEvents(const std::vector<JPetOpsEvent>& events){

  vector<JPetOpsEvent> newEventVec;

  int n_events = events.size();

  int n_prompt = 0;
  int n_anh = 0;

  std::vector<std::reference_wrapper<const JPetHit>> prompt_hits;
  std::vector<std::reference_wrapper<const JPetOpsEvent>> anh_events;
  
  for(int i=0;i<n_events;++i){
    const JPetOpsEvent & event = events.at(i);

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

      // fill annihilation points histo
      getStatistics().getHisto2D("anh_XY")->Fill(event.getAnnihilationPoint().Y(),
						 event.getAnnihilationPoint().X());
      if(fabs(event.getAnnihilationPoint().Z()) > 2.0 ){
	getStatistics().getHisto2D("anh_XY_nocenter")->Fill(event.getAnnihilationPoint().Y(),
							    event.getAnnihilationPoint().X());
      }
      
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

    if(fabs(new_event.getAnnihilationPoint().Z()) > 4.0 &&
       fabs(dt) > 20.0 * 1000. // ps
       ){
      getStatistics().getHisto2D("anh_XY_nocenter_lt_cut")->Fill(new_event.getAnnihilationPoint().Y(),
								 new_event.getAnnihilationPoint().X());
    }
  }else if(n_anh==1){ // 3g event did not have a prompt -> no lifetime
    JPetOpsEvent new_event = anh_events.front().get();
    new_event.setHasPrompt(false);
    new_event.setLifeTime(-1.e12);
    newEventVec.push_back(new_event);
  }
  
  return newEventVec;
}


std::vector<JPetOpsEvent> OPSAnalyzer::filterAnnihilationCandidates(const JPetTimeWindow& time_window){
  vector<JPetOpsEvent> newEventVec;

  int n_events = time_window.getNumberOfEvents();

  for(int i=0;i<n_events;++i){
    const JPetOpsEvent & event = time_window.getEvent<JPetOpsEvent>(i);

    if(!event.isTypeOf(JPetEventType::k3Gamma)){ // for pure prompt events,
      // simply rewrite them
      newEventVec.push_back(event);
      continue; 
    }

    Kinematics3g kinematics = calcKinematics(event);
    auto& angles = std::get<0>(kinematics); 
    auto& energies = std::get<1>(kinematics); 
    auto& momenta = std::get<2>(kinematics); 
    
    bool wrong_kinematics = false;
    if( energies[0] <= 0. || energies[1] <= 0. || energies[2] <= 0. ){
      wrong_kinematics = true;
    }

    getStatistics().getHisto1D("wrong_kinematics")->Fill(wrong_kinematics);

    if(wrong_kinematics){ // skip events with unphysical kinematics
      continue;
    }

    // fill angle histograms
    getStatistics().getHisto2D("3_hit_angles")->Fill(angles[0]+angles[1], angles[1]-angles[0]);
    
    // fill energy and momentum histograms
    getStatistics().getHisto1D("E_sum")->Fill(energies[0]+energies[1]+energies[2]);
    getStatistics().getHisto1D("P_sum")->Fill((momenta[0]+momenta[1]+momenta[2]).Mag2());

    getStatistics().getHisto2D("E2_E1")->Fill(energies[0], energies[1]);
    getStatistics().getHisto2D("E3_E2")->Fill(energies[1], energies[2]);
    getStatistics().getHisto2D("E3_E1")->Fill(energies[0], energies[2]);

    // shuffle the angles for the 2D relative plots
    std::random_shuffle(angles.begin(), angles.end());
    
    getStatistics().getHisto2D("th2_th1")->Fill(angles[0], angles[1]);
    getStatistics().getHisto2D("th3_th2")->Fill(angles[1], angles[2]);
    getStatistics().getHisto2D("th3_th1")->Fill(angles[0], angles[2]);
    
    newEventVec.push_back(event);
  }
  
  return newEventVec;
}



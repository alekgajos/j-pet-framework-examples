/**
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
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
 *  @file OPSCandidateFinder.cpp
 */

#include <iostream>
#include "OPSCandidateFinder.h"
#include "OPSAnalysisTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>

using namespace jpet_options_tools;
using namespace ops_analysis_tools;
using namespace std;

OPSCandidateFinder::OPSCandidateFinder(const char* name): JPetUserTask(name) {}

bool OPSCandidateFinder::init()
{
  INFO("Event finding started.");

  fOutputEvents = new JPetTimeWindow("JPetEvent");

  if (isOptionSet(fParams.getOptions(), fEventTimeParamKey)){
    kEventTimeWindow = getOptionAsFloat(fParams.getOptions(), fEventTimeParamKey);
  }else{
    ERROR("Event time window width not provided by the user!");
    return false;
  }

  for(auto key : fTOTcutKeys){
    if (isOptionSet(fParams.getOptions(), key)){
      fTOTcuts.push_back(getOptionAsFloat(fParams.getOptions(), key));
    }else{
      ERROR(Form("TOT cut value (%s) not provided by the user!", key.c_str()));
      return false;
    }
  }
  INFO(Form("Loaded TOT cut values: (%lf, %lf) and (%lf, %lf).",
	    fTOTcuts[0], fTOTcuts[1], fTOTcuts[2], fTOTcuts[3]));  

  if (isOptionSet(fParams.getOptions(), fAngleSumCutKey)){
    fAngleSumCut = getOptionAsFloat(fParams.getOptions(), fAngleSumCutKey);
  }else{
    ERROR("Angles sum cut value not provided by the user!");
    return false;
  }

  /************************************************************************/
  /* Multiplicity of particular hits                                      */
  /************************************************************************/
  getStatistics().createHistogram(
				  new TH1F("hits_per_event", "Number of Hits in event", 20, -0.5, 19.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("anh_hits_per_event",
					   "Number of annihilation candidate hits in Event", 10, -0.5, 9.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("dex_hits_per_event",
					   "Number of deexcitation candidate hits in Event", 10, -0.5, 9.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("other_hits_per_event",
					   "Number of other hits in Event", 10, -0.5, 9.5)
				  );  

  getStatistics().createHistogram(
				  new TH2F("dex_vs_anh_hits_per_event",
					   "Number of deexcitation vs annihilation candidate"
					   " hits in an event; annihilation; deexcitation",
					   10, -0.5, 9.5, 10, -0.5, 9.5)
				  );
  getStatistics().createHistogram(
				  new TH2F("dex_vs_other_hits_per_event",
					   "Number of deexcitation vs other"
					   " hits in an event; other; deexcitation",
					   10, -0.5, 9.5, 10, -0.5, 9.5)
				  );
  getStatistics().createHistogram(
                                  new TH2F("anh_vs_other_hits_per_event",
                                           "Number of annihilation vs other"
                                           " hits in an event; other; annihilation",
					   10, -0.5, 9.5, 10, -0.5, 9.5)
				  );

  /************************************************************************/
  /* Scatter vetoes                                                       */
  /************************************************************************/
  getStatistics().createHistogram(
				  new TH1F("were_same_strip",
					   "were there 2 hits in the same strip?", 2, -0.5, 1.5)
				  );  
  getStatistics().createHistogram(
				  new TH1F("evt_killed_by_2_same_strip",
					   "was the event rejected because of 2 hits in close strips?", 2, -0.5, 1.5)
				  );
  getStatistics().createHistogram(
				  new TH1F("were_close_strips",
					   "were there 2 hits in the same strip?", 2, -0.5, 1.5)
				  );  
  getStatistics().createHistogram(
				  new TH1F("evt_killed_by_2_close_strips",
					   "was the event rejected because of 2 hits in close strips?", 2, -0.5, 1.5)
				  );  

  getStatistics().createHistogram(
				  new TH1F("theta_diffs",
					   "#Delta #theta for hits in event",
					   181, -0.5, 180.5)
				  );

  getStatistics().createHistogram(
				  new TH1F("dvt",
					   "t-d/v;t-d/v [cm]",
					   200, -5., 5.)
				  );

  getStatistics().createHistogram(new TH2F("tot_rel_single_strip",
					   "Relative TOTs of 2 hits found in the same scintillator;"
					   "TOT 1 [ns]; TOT 2 [ns]",
					   1000, 0., 100.,
					   1000, 0., 100.)
				  );
  
  getStatistics().createHistogram(new TH2F("tot_rel_close_strips",
					   "Relative TOTs of 2 hits found in neighbouring strips;"
					   "TOT 1 [ns]; TOT 2 [ns]",
					   1000, 0., 100.,
					   1000, 0., 100.)
				  );
  
  /**************************************************************************/
  /* TOTs                                                                   */
  /**************************************************************************/
  getStatistics().createHistogram(new TH1F("tot4_selected", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));

  getStatistics().createHistogram(new TH1F("tot4_3+hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_2hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_1hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  getStatistics().createHistogram(new TH1F("tot4_any_hits_nocut", "Sum of tot-s at 4 thresholds;TOT [ns]", 1000, 0., 100.));
  
  /************************************************************************/
  /* Relative angles histogram for 3-hit events                           */
  /************************************************************************/
  getStatistics().createHistogram(
				  new TH2F("3_hit_angles",
					   "3 Hit angles difference 1-2, 2-3",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");

  // same as above, but only for passed events
  getStatistics().createHistogram(
				  new TH2F("3_hit_angles_passed",
					   "3 Hit angles difference 1-2, 2-3: passed events",
					   360, -0.5, 359.5,
					   360, -0.5, 359.5)
				  );
  getStatistics().getHisto2D("3_hit_angles_passed")
    ->GetXaxis()->SetTitle("Smallest angle + Second smallest angle [deg]");
  getStatistics().getHisto2D("3_hit_angles_passed")
    ->GetYaxis()->SetTitle("Second smallest angle - Smallest angle [deg]");  


  getStatistics().createHistogram(
                                  new TH1F("refined_events",
                                           "No. refined events in TW",
                                           20, -0.5, 19.5
                                           )
                                  );
  getStatistics().createHistogram(
                                  new TH1F("event_candidates_tw",
                                           "No. event candidates in TW",
                                           20, -0.5, 19.5
                                           )
                                  );
  
  return true;
}

bool OPSCandidateFinder::exec()
{

  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    fTimeWindowHadAnnihilation = false;
    std::vector<JPetEvent> events = findEventCandidates(refineEvents(*timeWindow));
    getStatistics().getHisto1D("event_candidates_tw")->Fill(events.size());
    //    events = vetoScatterings( events );

    // only save anything from a time window
    // if at least one 3g event was identified there
    // otherwise, leave an empty time window
    if(fTimeWindowHadAnnihilation){ 
      saveEvents(events);
    }

  } else {
    return false;
  }
  return true;
}

bool OPSCandidateFinder::terminate()
{
  INFO("Identification of o-Ps->3g candidate events ended.");
  
  return true;
}

std::vector<JPetEvent> OPSCandidateFinder::refineEvents(const JPetTimeWindow& preEvents)
{

  vector<JPetEvent> newEventVec;

  int nevents = preEvents.getNumberOfEvents();
  
  for(int entry=0; entry<nevents; ++entry){
    const JPetEvent& event = dynamic_cast<const JPetEvent&>(preEvents[entry]);
    
    const auto & hits = event.getHits();
    
    int s = 0;
    int nhits = event.getHits().size();
    
    while ( s < nhits ) {

      JPetEvent newEvent;
      newEvent.setEventType(JPetEventType::kUnknown);

      JPetHit startHit = hits[s];

      //if( identifyHitType(startHit, fTOTcuts) != HitCandidateType::None )
      
      newEvent.addHit(startHit);

      int k = 1;
      while ( s + k < nhits ) {
	JPetHit currentHit = hits[s + k];
	if (fabs(currentHit.getTime() - startHit.getTime()) < kEventTimeWindow) {
	  newEvent.addHit(currentHit);
	  k++;
	} else {
	  break;
	}
      }
      s += k;
      
      newEventVec.push_back(newEvent);
    }
  }

  getStatistics().getHisto1D("refined_events")->Fill(newEventVec.size());
  return newEventVec;
}

std::vector<JPetEvent> OPSCandidateFinder::findEventCandidates(const std::vector<JPetEvent> events){

  std::vector<JPetEvent> newEventVec;

  for (const auto & event : events) {

    int n_anh_hits = 0;
    int n_prompt_hits = 0;
    int n_other_hits = 0;

    getStatistics().getHisto1D("hits_per_event")->Fill(event.getHits().size());
    
    // count hits of particular types
    for(const auto& hit : event.getHits()){
      HitCandidateType type = identifyHitType(hit, fTOTcuts);
      if(type==HitCandidateType::Annihilation){
        n_anh_hits++;
      }else if(type==HitCandidateType::Prompt){
        n_prompt_hits++;
      }else{
        n_other_hits++;
      } 
    }

    // fill control histograms
    getStatistics().getHisto1D("dex_hits_per_event")->Fill(n_prompt_hits);
    getStatistics().getHisto1D("anh_hits_per_event")->Fill(n_anh_hits);
    getStatistics().getHisto1D("other_hits_per_event")->Fill(n_other_hits);

    getStatistics().getHisto2D("anh_vs_other_hits_per_event")->Fill(n_other_hits, n_anh_hits);
    getStatistics().getHisto2D("dex_vs_other_hits_per_event")->Fill(n_other_hits, n_prompt_hits);
    getStatistics().getHisto2D("dex_vs_anh_hits_per_event")->Fill(n_anh_hits, n_prompt_hits);

    fillTOThistos(event);
    
    JPetEvent new_event = event;
    
    // set event type flags
    if( n_anh_hits >= 3 ){
      new_event.setEventType(JPetEventType::k3Gamma);
      fTimeWindowHadAnnihilation = true; // used to decide whether to storea anything from the time window or not
    }
    
    if( n_prompt_hits >= 1 ){
      new_event.addEventType(JPetEventType::kPrompt);
    }
    
    // filter events
    if( new_event.isTypeOf(JPetEventType::k3Gamma) || 
        new_event.isTypeOf(JPetEventType::kPrompt)){

      // temporary
      // remove hits other than annihilation and prompt
      std::vector<JPetHit> new_hits;
      for(const auto& hit : event.getHits()){
        HitCandidateType type = identifyHitType(hit, fTOTcuts);
        if(type==HitCandidateType::Annihilation || type==HitCandidateType::Prompt){
          JPetHit new_hit(hit);
          double qual = -1.0;
          if(type==HitCandidateType::Annihilation) qual = 1.0;
          else if(type==HitCandidateType::Prompt) qual = 2.0;
          new_hit.setQualityOfEnergy(qual);
          new_hits.push_back(new_hit);
          getStatistics().getHisto1D("tot4_selected")->Fill(hit.getEnergy());
        }
      }
      new_event.setHits(new_hits);

      newEventVec.push_back(new_event);
    }

  }

  return newEventVec;
}

void OPSCandidateFinder::saveEvents(const std::vector<JPetEvent>& events)
{
  for (const auto & event : events) {

    // study angles for 3-hit events
    //    if( analyseThreeHitEvent(event) ){
    fOutputEvents->add<JPetEvent>(event);
    //      }
  }
  
}


bool OPSCandidateFinder::analyseThreeHitEvent(const JPetEvent& event){

  vector<JPetEvent> eventVector2;
  vector<JPetHit> hits = event.getHits();

  if(hits.size() > 3){
    
    hits.erase(
	       remove_if(hits.begin(), hits.end(),
			 [](const JPetHit& hit){
			   if( hit.getQualityOfEnergy() > 0.5 ){
			     return true;
			   }
			   return false;
			 }),
	       hits.end());
    
  } 
  
  JPetHit firstHit = event.getHits().at(0);
  JPetHit secondHit = event.getHits().at(1);
  JPetHit thirdHit = event.getHits().at(2);

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
  getStatistics().getHisto2D("3_hit_angles")->Fill(angles[0]+angles[1], angles[1] - angles[0]);

  // only take o-Ps -> 3g decay candidates
  if(angles[0]+angles[1] > fAngleSumCut){
    getStatistics().getHisto2D("3_hit_angles_passed")->Fill(angles[0]+angles[1], angles[1] - angles[0]);
    return true;
  }

  return false;
}


std::vector<JPetEvent> OPSCandidateFinder::vetoScatterings(const std::vector<JPetEvent> events){

  std::vector<JPetEvent> returnEvents;
  
  for(auto & event : events){

    if(!event.isTypeOf(JPetEventType::k3Gamma)){
      //      continue;
    }
    
    auto hits = event.getHits();
    std::set<std::vector<JPetHit>::iterator> to_erase;
    if(hits.size()<3)continue;
    double min_dvt = -10000.;

    /**********************************************************************/
    /* Rejecting hits from the same strip                                 */
    /**********************************************************************/
    for(auto i=hits.begin(); i!=hits.end();i++){
      for(auto j=i+1; j!=hits.end();j++){
    
	auto & hit1 = *i;
	auto & hit2 = *j;

	// reject hits in the same scintillators
	if( hit1.getBarrelSlot() == hit2.getBarrelSlot() ){

	  // study the relative TOT of such hits on the same module
	  getStatistics().getHisto2D("tot_rel_single_strip")->Fill(hit1.getEnergy(), hit2.getEnergy());
          
	  to_erase.insert(i);
	  to_erase.insert(j);
	  break;
	}
      }
    }

    getStatistics().getHisto1D("were_same_strip")->Fill(!to_erase.empty());
    
    // actually remove the vetoed hits
    for(const auto & it: to_erase){
      hits.erase(it);
    }
    if(hits.size() < 3){
      getStatistics().getHisto1D("evt_killed_by_2_same_strip")->Fill(1);
      continue;
    }
    getStatistics().getHisto1D("evt_killed_by_2_same_strip")->Fill(0);
    
    
    /*********************************************************************/
    /* Rejecting annihilation candidate hits with close angles           */
    /*********************************************************************/
    to_erase.clear();
    for(auto i=hits.begin(); i!=hits.end();i++){
      for(auto j=i+1; j!=hits.end();j++){
    
	auto & hit1 = *i;
	auto & hit2 = *j;
        
        double d_theta = fabs(hit1.getBarrelSlot().getTheta() - hit2.getBarrelSlot().getTheta());
        if( d_theta > 180. ){
          d_theta = 360.0 - d_theta;
        }
        getStatistics().getHisto1D("theta_diffs")->Fill(d_theta);

        if( d_theta < 8.0 ){

          // study the relative TOT of such close hits
	  getStatistics().getHisto2D("tot_rel_close_strips")->Fill(hit1.getEnergy(), hit2.getEnergy());

	  to_erase.insert(i);
	  to_erase.insert(j);
	  break;          
        }
      }
    }
    getStatistics().getHisto1D("were_close_strips")->Fill(!to_erase.empty());

    // actually remove the vetoed hits
    for(const auto & it: to_erase){
      hits.erase(it);
    }
    if(hits.size() < 3){
      getStatistics().getHisto1D("evt_killed_by_2_close_strips")->Fill(1);
      continue;
    }
    getStatistics().getHisto1D("evt_killed_by_2_close_strips")->Fill(0);
    
    /****************************************************************/
    /* Identifying scatterings based on the d - vt                  */
    /****************************************************************/
    for(auto i=hits.begin(); i!=hits.end();i++){
      for(auto j=i+1; j!=hits.end();j++){

	auto & hit1 = *i;
	auto & hit2 = *j;
        
        double d = (hit1.getPos() - hit2.getPos()).Mag();
        double dt = fabs(hit1.getTime() - hit2.getTime()) / 1000.;
        double dvt = dt - d/kSpeedOfLight;

        if( fabs(dvt) < fabs(min_dvt) ){
          min_dvt = dvt;
        }
	  	 
      }      
    } 

    getStatistics().getHisto1D("dvt")->Fill(min_dvt);	    

    if( min_dvt > -1.8 ){
      continue;
    }
    
    JPetEvent new_event = event;
    new_event.setHits(hits);
    returnEvents.push_back(new_event);
      
  } // end loop over events
  
  return returnEvents;
}


void OPSCandidateFinder::fillTOThistos(const JPetEvent & event) {

    // Filling of histograms
    for(auto & hit: event.getHits()){
      double tot = hit.getEnergy();
      
      // fill histos for all events
      getStatistics().getHisto1D("tot4_any_hits_nocut")->Fill(tot);

      // fill histos only 1-hit events
      if(event.getHits().size()==1){
	getStatistics().getHisto1D("tot4_1hits_nocut")->Fill(tot);
      }

      // fill histos only for 2-hit events
      if(event.getHits().size()==2){
	getStatistics().getHisto1D("tot4_2hits_nocut")->Fill(tot);
      }

      // fill histos for 3+ hits events
      if(event.getHits().size()>=3){
	getStatistics().getHisto1D("tot4_3+hits_nocut")->Fill(tot);
      }
    }
}

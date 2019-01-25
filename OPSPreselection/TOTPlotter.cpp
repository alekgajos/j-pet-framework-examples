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
 *  @file TOTPlotter.cpp
 */

#include <iostream>
#include <JPetWriter/JPetWriter.h>
#include "TOTPlotter.h"
#include <JPetOptionsTools/JPetOptionsTools.h>

using namespace jpet_options_tools;

using namespace std;

TOTPlotter::TOTPlotter(const char* name): JPetUserTask(name) {}

JPetHit TOTPlotter::calculateTOT(const JPetHit& hit){

  double tot = 0.;

  for(auto& phys_sig: {hit.getSignalA(), hit.getSignalB()}){

    const JPetRawSignal& raw_sig = phys_sig.getRecoSignal().getRawSignal();
    for(auto& p: raw_sig.getTimesVsThresholdValue(JPetSigCh::Leading)){
      // fill multiplicity of particular threshold values at leading edge
      getStatistics().getHisto1D("thr_mult")->Fill(p.second.first);
    }

    // calculate single-signal TOT
    auto leading_points = raw_sig.getTimesVsThresholdNumber(JPetSigCh::Leading);
    auto trailing_points = raw_sig.getTimesVsThresholdNumber(JPetSigCh::Trailing);
      
    for(int i=1;i<5;i++){
      auto lead_search = leading_points.find(i);
      auto trail_search = trailing_points.find(i);
      if (lead_search != leading_points.end()
          && trail_search != trailing_points.end()){
        tot += (trail_search->second - lead_search->second) / 1000.; // in ns
      }
    }
    
  }
  
  JPetHit new_hit = hit;
  new_hit.setEnergy(tot);
  return new_hit;
}

bool TOTPlotter::init()
{

  INFO("TOT plotting started");

  fOutputEvents = new JPetTimeWindow("JPetHit");

  for(int i=1; i <= getParamBank().getBarrelSlotsSize(); ++i){

    getStatistics().createHistogram(
				    new TH1F(Form("tot_strip_good_%d", getParamBank().getBarrelSlot(i).getID()),
                                             "TOT for good hits; TOT [ns]; counts", 1000., 0., 100.)
				    );
    getStatistics().createHistogram(
				    new TH1F(Form("tot_strip_corrupted_%d", getParamBank().getBarrelSlot(i).getID()),
                                             "TOT for corrupted hits; TOT [ns]; counts", 1000., 0., 100.)
				    );

  }

  getStatistics().createHistogram(new TH1F("thr_mult", "threshold value multiplicity in raw signals",
                                           100, 0., 500.)
                                  );
  
  return true;
}

bool TOTPlotter::exec()
{
  auto timeWindow = dynamic_cast<const JPetTimeWindow*>(fEvent);

  uint n = timeWindow->getNumberOfEvents();

  for(uint i=0;i<n;++i){

    const JPetHit & hit =  dynamic_cast<const JPetHit&>(timeWindow->operator[](i));
    JPetHit new_hit = calculateTOT(hit);
    
    double tot = new_hit.getEnergy();

    if(hit.getRecoFlag() == JPetHit::Good){
      getStatistics().getHisto1D(Form("tot_strip_good_%d", new_hit.getBarrelSlot().getID()))->Fill(tot);
    }else{
      getStatistics().getHisto1D(Form("tot_strip_corrupted_%d", new_hit.getBarrelSlot().getID()))->Fill(tot);
    }

    fOutputEvents->add<JPetHit>(new_hit);
  }

  return true;
}

bool TOTPlotter::terminate()
{
  INFO("Done plotting TOT-s");
  return true;
}


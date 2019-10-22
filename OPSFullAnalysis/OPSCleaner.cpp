/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
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
 *  @file OPSCleaner.cpp
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <functional>
#include "OPSCleaner.h"
#include "OPSAnalysisTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>

using namespace jpet_options_tools;
using namespace ops_analysis_tools;

using namespace std;


OPSCleaner::OPSCleaner(const char* name): JPetUserTask(name) {}

void OPSCleaner::bookHisto(TH1* h){

  std::string base_name = std::string(h->GetName());
  for(int step=0;step<kNsteps;++step){
    std::string name = base_name + std::string(Form("_%d", step));
    getStatistics().createHistogram(h->Clone(name.c_str()));
  }
}

TH1F* OPSCleaner::getHisto1D(std::string name, int step){
  
  std::string full_name = name + std::string(Form("_%d", step));
  return getStatistics().getHisto1D(full_name.c_str());
}

TH2F* OPSCleaner::getHisto2D(std::string name, int step){
  
  std::string full_name = name + std::string(Form("_%d", step));
  return getStatistics().getHisto2D(full_name.c_str());
}

bool OPSCleaner::init()
{
  INFO("Cleaning of o-Ps->3g candidate sample started.");

  fOutputEvents = new JPetTimeWindow("JPetOpsEvent");

  if (isOptionSet(fParams.getOptions(), fAngleSumCutKey)){
    fAngleSumCut = getOptionAsFloat(fParams.getOptions(), fAngleSumCutKey);
  }else{
    ERROR("Angles sum cut value not provided by the user!");
    return false;
  }

  // initialize counters
  for(int i=0;i<10;++i){
    fEventCouters[i] = 0;
  }

  /************************************************************************/
  /* New histograms, to be filled after every step of the cleaning        */
  /************************************************************************/
  bookHisto(new TH1F("min_angle_2d",
                     "Minimal #theta angle difference between hit scintillators (XY);"
                     "#theta_{MIN} [deg]",
                     180, -0.5, 179.5)
            );

  bookHisto(new TH1F("min_angle_3d",
                     "Minimal angle between photons 3D;"
                     "#theta^{3D}_{MIN} [deg]",
                     180, -0.5, 179.5)
            );

  bookHisto(
            new TH2F("E2_E1",
                     "E2 vs E1;E1 [keV];E2 [keV]",
                     600, -0.5, 599.5,
                     600, -0.5, 599.5)
            );
  
  bookHisto(
            new TH2F("E3_E1",
                     "E3 vs E1;E1 [keV];E3 [keV]",
                     600, -0.5, 599.5,
                     600, -0.5, 599.5)
            );
  
  bookHisto(
            new TH2F("E3_E2",
                     "E3 vs E2;E2 [keV];E3 [keV]",
                     600, -0.5, 599.5,
                     600, -0.5, 599.5)
            );
  
  bookHisto(
            new TH2F("th2_th1",
                     "th2 vs th1;th1 [deg];th2 [deg]",
                     200, -0.5, 189.5,
                     200, -0.5, 189.5)
            );

  bookHisto(
            new TH2F("th3_th1",
                     "th3 vs th1;th1 [deg];th3 [deg]",
                     200, -0.5, 189.5,
                     200, -0.5, 189.5)
            );
  
  bookHisto(
            new TH2F("th3_th2",
                     "th3 vs th2;th2 [deg];th3 [deg]",
                     200, -0.5, 189.5,
                     200, -0.5, 189.5)
            );
  
  bookHisto(
            new TH2F("theta_angles",
                     "Theta angle differences between hit scintillators;"
                     "Smallest angle + Second smallest angle [deg];"
                     "Second smallest angle - Smallest angle [deg]",
                     360, -0.5, 359.5,
                     360, -0.5, 359.5)
            );

  bookHisto(
            new TH2F("3_hit_angles",
                     "3 Hit 3D angles difference;"
                     "Smallest angle + Second smallest angle [deg];"
                     "Second smallest angle - Smallest angle [deg]",
                     360, -0.5, 359.5,
                     360, -0.5, 359.5)
            );

  bookHisto(
            new TH2F("anh_XY",
                     "transverse position of the o-Ps->3g decay point;"
                     "Y [cm]; X [cm]",
                     100, -50.5, 49.5,
                     100, -50.5, 49.5
                     )
            );
  
  bookHisto(
            new TH2F("anh_XZ",
                     "position of the o-Ps->3g decay point in XZ;"
                     "Z [cm]; X [cm]",
                     100, -50., 50.,
                     100, -50., 50.
                     )
            );  

  bookHisto(new TH1F("anh_R",
                     "Transverse radius of the o-Ps->3g decay point; R [cm]",
                     100, 0., 50.
                     )
            );

  bookHisto(new TH1F("anh_R_jacobian",
                     "Transverse radius of the o-Ps->3g decay point - Jacobian weighted; R [cm]",
                     100, 0., 50.
                     )
            );
  
  bookHisto(
            new TH1F("lifetime",
                     "Time between deexcitation and annihilation;"
                     "#Delta t [ns]",
                     2000, -260., 260.)
            );  

  bookHisto(
            new TH1F("Sk1",
                     "|S|#bullet |k_{1}|;|S|#bullet |k_{1}|",
                     1000, -1.1, 1.1)
            );

  bookHisto(
            new TH1F("Sk1xk2",
                     "|S|#bullet |k_{1}#times k_{2}|;|S|#bullet|k_{1}#times k_{2}|",
                     1000, -1.1, 1.1)
            );

  bookHisto(
            new TH1F("Sk1.Sk1xk2",
                     "(|S|#bullet |k_{1}|)(|S|#bullet |k_{1}#times k_{2}|);"
                     "(|S|#bullet |k_{1}|)(|S|#bullet |k_{1}#times k_{2}|)",
                     1000, -1.1, 1.1)
            );  

  bookHisto(
            new TH1F("lors_d_min",
                     "Min dist of 2#gamma vtx to 3#gamma vtx; d_{MIN}^{LOR} [cm]",
                     400, 0., 200.)
            );  

  bookHisto(
            new TH1F("lors_d_max",
                     "Max dist of 2#gamma vtx to 3#gamma vtx; d_{MAX}^{LOR} [cm]",
                     400, 0., 200.)
            );  

  bookHisto(
            new TH1F("lors_d_total",
                     "Total dist of 2#gamma vertices to 3#gamma vtx; d_{TOTAL}^{LOR} [cm]",
                     400, 0., 200.)
            );  

  bookHisto(
            new TH2F("lors_d_max_vs_min",
                     "Max vs min dist of 2#gamma vertices to 3#gamma vtx;"
                     "d_{MIN}^{LOR} [cm];"
                     "d_{MAX}^{LOR} [cm]",
                     400, 0., 200.,
                     400, 0., 200.)
            );

  bookHisto(
            new TH2F("lors_d_total_vs_min",
                     "Total vs min dist of 2#gamma vertices to 3#gamma vtx;"
                     "d_{MIN}^{LOR} [cm];"
                     "d_{TOTAL}^{LOR} [cm]",
                     400, 0., 200.,
                     400, 0., 200.)
            );

  bookHisto(
            new TH2F("lors_d_total_vs_max",
                     "Total vs max dist of 2#gamma vertices to 3#gamma vtx;"
                     "d_{MIN}^{LOR} [cm];"
                     "d_{TOTAL}^{LOR} [cm]",
                     400, 0., 200.,
                     400, 0., 200.)
            );

  bookHisto(
            new TH2F("dvts",
                     "d-ct;|d-ct|_{1} [cm]; |d-ct|_{2} [cm]",
                     400, 0., 100.,
                     400, 0., 100.)
            );  

return true;
}

bool OPSCleaner::exec()
{

  if (auto time_window = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {

    uint n_events = time_window->getNumberOfEvents();
    for(int i=0;i<n_events;++i){

      const JPetOpsEvent & event = time_window->getEvent<JPetOpsEvent>(i);

      EventInfo evt_info;
      evt_info.event_odd = (i%2==1);
      evt_info.kinematics = calcKinematics(event);
      evt_info.angles_xy = calcSmallestAnglesXY(event);
      auto& angles = std::get<0>(evt_info.kinematics); 
      evt_info.angles_3d = std::make_pair(angles[0], angles[1]);
      // shuffle the angles for the 2D relative plots
      std::random_shuffle(angles.begin(), angles.end());

      auto& energies = std::get<1>(evt_info.kinematics); 
      
      // calculate the operators
      evt_info.operators = calcOperators(event, evt_info.kinematics);

      // Calculate the LOR distances
      evt_info.distances = calcLORdistances(event);

      // Perform scatter tests
      evt_info.scatter_tests = calcScatterTests(event);
      
      /********************************************************************/
      /* Study of particular cuts starts here                             */
      /********************************************************************/
      fillHistos(event, evt_info, 0);

      // cut on minimal 3D angle between photons' momenta
      if( evt_info.angles_3d.first > 35.0 ){
        fillHistos(event, evt_info, 1);
      }

      if( evt_info.scatter_tests.at(0) > 15.0 ){
        fillHistos(event, evt_info, 2);

        if( evt_info.angles_3d.second + evt_info.angles_3d.first > 190.0){
          fillHistos(event, evt_info, 3);
        }

        if( energies[0] > 20.0 &&
            energies[1] > 20.0 &&
            energies[2] > 20.0
            ){
          fillHistos(event, evt_info, 4);
        }

        if( event.getLifeTime()/1000. > 20.0 && event.getLifeTime()/1000. < 150. ){
          fillHistos(event, evt_info, 5);
        }
        
      }

      if( energies[0] > 20.0 &&
          energies[1] > 20.0 &&
          energies[2] > 20.0
          ){
        fillHistos(event, evt_info, 6);
      }

      if( evt_info.scatter_tests.at(0) > 15.0 &&
          evt_info.angles_3d.second + evt_info.angles_3d.first > 190.0 &&
          energies[0] > 20.0 &&
          energies[1] > 20.0 &&
          energies[2] > 20.0 &&
          event.getLifeTime()/1000. > 20.0 &&
          event.getLifeTime()/1000. < 150.
          ){

        fillHistos(event, evt_info, 7);
      }

    } // end loop over events

  } else {
    return false;
  }
  return true;
}

bool OPSCleaner::terminate()
{
  INFO("Cleaning of o-Ps->3g candidate sample finished.");
  
  return true;
}

void OPSCleaner::fillHistos(const JPetOpsEvent& event, EventInfo evt_info, int step){

  auto& angles = std::get<0>(evt_info.kinematics); 
  auto& energies = std::get<1>(evt_info.kinematics); 

  // fill angle histograms
  getHisto1D("min_angle_2d", step)->Fill(evt_info.angles_xy.first);
  getHisto1D("min_angle_3d", step)->Fill(evt_info.angles_3d.first);
  
  getHisto2D("3_hit_angles", step)->Fill(evt_info.angles_3d.second + evt_info.angles_3d.first,
                                          evt_info.angles_3d.second - evt_info.angles_3d.first);
  getHisto2D("theta_angles", step)->Fill(evt_info.angles_xy.second + evt_info.angles_xy.first,
                                          evt_info.angles_xy.second - evt_info.angles_xy.first);

  if(evt_info.event_odd){
    getHisto2D("th2_th1", step)->Fill(angles[0], angles[1]);
    getHisto2D("th3_th2", step)->Fill(angles[1], angles[2]);
    getHisto2D("th3_th1", step)->Fill(angles[0], angles[2]);
  }else{
    getHisto2D("th2_th1", step)->Fill(angles[1], angles[0]);
    getHisto2D("th3_th2", step)->Fill(angles[2], angles[1]);
    getHisto2D("th3_th1", step)->Fill(angles[2], angles[0]);
  }
    
  // fill energy histograms
  if(evt_info.event_odd){
    getHisto2D("E2_E1", step)->Fill(energies[0], energies[1]);
    getHisto2D("E3_E2", step)->Fill(energies[1], energies[2]);
    getHisto2D("E3_E1", step)->Fill(energies[0], energies[2]);
  }else{
    getHisto2D("E2_E1", step)->Fill(energies[1], energies[0]);
    getHisto2D("E3_E2", step)->Fill(energies[2], energies[1]);
    getHisto2D("E3_E1", step)->Fill(energies[2], energies[0]);
  }
  
  // lifetime
  if(event.hasPrompt()){
    getHisto1D("lifetime", step)->Fill(event.getLifeTime() / 1000.);
  }

  // annihilation point location
  getHisto2D("anh_XZ", step)->Fill(event.getAnnihilationPoint().Z(), event.getAnnihilationPoint().X());

  if( fabs(event.getAnnihilationPoint().Z()) > 4.0 ){
    getHisto2D("anh_XY", step)->Fill(event.getAnnihilationPoint().Y(), event.getAnnihilationPoint().X());

    double r = event.getAnnihilationPoint().Perp();
    getHisto1D("anh_R", step)->Fill(r);
    getHisto1D("anh_R_jacobian", step)->Fill(r, 1./r);
  }
    
  // distances
  getHisto1D("lors_d_min", step)->Fill(evt_info.distances.at(0));
  getHisto1D("lors_d_max", step)->Fill(evt_info.distances.at(2));
  getHisto1D("lors_d_total", step)->Fill(evt_info.distances.at(3));
  getHisto2D("lors_d_max_vs_min", step)->Fill(evt_info.distances.at(0),
                                              evt_info.distances.at(2));
  getHisto2D("lors_d_total_vs_min", step)->Fill(evt_info.distances.at(0),
                                                evt_info.distances.at(3));
  getHisto2D("lors_d_total_vs_max", step)->Fill(evt_info.distances.at(2),
                                                evt_info.distances.at(3));
  
  // operators
  getHisto1D("Sk1", step)->Fill(evt_info.operators.at(0));
  getHisto1D("Sk1xk2", step)->Fill(evt_info.operators.at(1));
  getHisto1D("Sk1.Sk1xk2", step)->Fill(evt_info.operators.at(2));

  // scatter tests
  getHisto2D("dvts", step)->Fill(evt_info.scatter_tests.at(0), evt_info.scatter_tests.at(1));
  
}





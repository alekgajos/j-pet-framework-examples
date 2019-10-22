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
 *  @file OPSAnalysisTools.cpp
 */

#include "OPSAnalysisTools.h"
#include "JPetLoggerInclude.h"
#include "../LargeBarrelAnalysis/EventCategorizerTools.h"
#include <TVector3.h>
#include <TMath.h>
#include <cmath>
#include <cassert>
#include <deque>

namespace ops_analysis_tools{

  namespace{

    /**
     * @brief Calculates the relative angle between i-th and j-th photon in an event
     *  (numbering starting from 1) in radians.
     */
    double relativeAngle(std::vector<TVector3> momenta, int i, int j){
      return  momenta[i-1].Angle(momenta[j-1]);
    }
  }

  const double kSpeedOfLight = 29.9792458;
  
  /**
     @brief Identify if the hit could correspond to annihilation photon,
     prompt photon or none of them.

     Optionally, this function also discriminates corrupted hits 
     by identifying then as None
     
     @param tot_cuts vector of 4 real numbers describing TOT cut boundaries:
     tot_cuts[0] - lower annihilation photon TOT cut
     tot_cuts[1] - upper annihilation photon TOT cut
     tot_cuts[2] - lower prompt photon TOT cut
     tot_cuts[3] - upper prompt photon TOT cut
  */
  HitCandidateType identifyHitType(const JPetHit& hit, std::vector<double>& tot_cuts){
    double tot = hit.getEnergy();

    // if hit was corrupted, do not identify its type
    if( hit.getRecoFlag() == JPetHit::Corrupted ){
      return HitCandidateType::None;
    }

    // check annihilation TOT cuts
    if( tot > tot_cuts[0] && tot < tot_cuts[1] ){
      return HitCandidateType::Annihilation;
    }
    if( tot > tot_cuts[2] && tot < tot_cuts[3] ){
      return HitCandidateType::Prompt;
    }
    return HitCandidateType::None;
  }

  /**
   * @brief Calculates and returns (smallest angle, second smallest angle) in the transverse plane. Uses theta angle positions of scintillator strips w.r.t. detector center.
   */
  SmallestAngles calcSmallestAnglesXY(const JPetOpsEvent& event){

    std::deque<float> theta_angles;
    
    for(auto& hit: event.getHits()){
      if( hit.getQualityOfEnergy() > 0.7 && hit.getQualityOfEnergy() < 1.3 ){ // annihilation hit

	double theta = hit.getBarrelSlot().getTheta();
	// already collect the angles in sorted order
	if( theta_angles.empty() ){
	  theta_angles.push_front(theta);
	}else if( theta <= theta_angles.front() ){
	  theta_angles.push_front(theta);	    
	}else if( theta >= theta_angles.back()){
	  theta_angles.push_back(theta);	    
	}else{
	  theta_angles.insert(theta_angles.begin()+1, theta);
	}
	
      }
    }

    if( theta_angles.size() != 3 ){
      WARNING(Form("The number of annihlation candidate"
		   " hits is %lu instead of expected 3.", theta_angles.size() ));
      return std::make_pair(-1000., -1000.);
    }
    
    std::deque<float> theta_diffs;
    for(int i=0;i<3;++i){
      for(int j=i+1;j<3;++j){
	double d_theta = theta_angles[j] - theta_angles[i];
	if(i==0 && j==2){
	  d_theta = 360.0 - d_theta;
	}
	if( theta_diffs.empty() ){
	  theta_diffs.push_front(d_theta);
	}else if( d_theta <= theta_diffs.front() ){
	  theta_diffs.push_front(d_theta);	    
	}else if( d_theta >= theta_diffs.back()){
	  theta_diffs.push_back(d_theta);	    
	}else{
	  theta_diffs.insert(theta_diffs.begin()+1, d_theta);
	}
      }
    }

    
    
    return std::make_pair(theta_diffs[0], theta_diffs[1]);
  }

  /**
   * @brief Calculates photons' momentum vectors, energies and angles between the momenta in 3D
   * 
   * returns:
   *  * a vector of 3 sorted relative angles in ascending order
   *  * a vector of 3 photon energies - unsorted
   *  * a vector of photons' momenta in sorted in order of ascending magnitude
   *
   * Note that therefore the returned vectors of angles, energies 
   * and momenta will be ordered differently!
   *
   * Do not rely on the momenta vectors if any of the calculated energies is negative!
   *
   */
  Kinematics3g calcKinematics(const JPetOpsEvent& event){

    std::vector<TVector3> momenta;
    momenta.reserve(3);
    std::vector<double> angles(3, -1000);
    std::vector<double> energies(3, -1000);
    
    for(auto& hit: event.getHits()){
      if( hit.getQualityOfEnergy() > 0.7 && hit.getQualityOfEnergy() < 1.3 ){ // annihilation hit
	momenta.push_back(hit.getPos()-event.getAnnihilationPoint());
      }
    }    

    if( momenta.size() != 3 ){
      WARNING(Form("The number of annihlation candidate"
		   " hits is %lu instead of expected 3.", momenta.size() ));
      return std::make_tuple(angles, energies, momenta);
    }

    // calculate the angles (in radians)
    double th12 = relativeAngle(momenta, 1, 2);
    double th23 = relativeAngle(momenta, 2, 3);
    double th31 = relativeAngle(momenta, 3, 1);

    // calculate energies of the photons
    double m_e = 510.998928; // electron mass in keV
    energies[0] = -2.*m_e*(-cos(th31)+cos(th12)*cos(th23))/((-1+cos(th12))*(1.+cos(th12)-cos(th23)-cos(th31)));
    energies[1] = -2.*m_e*(cos(th12)*cos(th31)-cos(th23))/((-1+cos(th12))*(1.+cos(th12)-cos(th23)-cos(th31)));
    energies[2] = 2.*m_e*(1+cos(th12))/(1.+cos(th12)-cos(th23)-cos(th31));
    
    // scale the momentum vectors according to the energy
    for(int i=0;i<3;++i){
      momenta[i].SetMag(energies[i]);
    }

    // sort the momentum vectors in the order of ascending energy
    std::sort(momenta.begin(), momenta.end(),
    	      [](const TVector3& a, const TVector3& b){
    		return a.Mag() < b.Mag();
    	      }
    	      );

    // sort the angles in ascending order and transform to degrees
    angles[0] = th12*TMath::RadToDeg();
    angles[1] = th23*TMath::RadToDeg();
    angles[2] = th31*TMath::RadToDeg();
    std::sort(angles.begin(), angles.end());
    
    // correct the maximum angle (this correction has no effect on energies' calculation)
    double epsilon = 1.e-4;
    if(fabs(angles[2] - (angles[0]+angles[1])) < epsilon){
      angles[2] = 360. - (angles[0] + angles[1]);
    }
    
    return std::make_tuple(angles, energies, momenta);
  }


  /**
   * @brief Calculates values of the symmetry-odd operators
   *
   * Takes Kinematics3g as a parameter which must be obtained in advance 
   * for the event with the `calcKinematics` function.
   *
   * Contents of the returned vector are:
   * [0] - S*k1
   * [1] - S*(k1xk2)
   * [2] - (S*k1)*(S*(k1xk2))
   *
   */ 
  Operators calcOperators(const JPetOpsEvent& event, Kinematics3g kinematics){
    Operators operators;
    
    TVector3 S = event.getAnnihilationPoint().Unit();
    
    auto& momenta = std::get<2>(kinematics); 
    TVector3& k1 = momenta.at(2);
    TVector3& k2 = momenta.at(1);
    
    // S*k1
    operators.push_back(
                        S.Dot(k1.Unit())
                        );
  
    // S*(k1xk2)
    operators.push_back(S.Dot( k1.Cross(k2).Unit() ) );

    // (S*k1)*(S*(k1xk2))
    operators.push_back( operators.at(0) * operators.at(1) );
    
    return operators;
  }


  /**
   * @brief Calculates distances of hypothetical 2g annihilation points
   * from the 3g annihilation point for all 3 LOR-s conceiveable in a 3-hit event.
   *
   * Contents of the returned vector are:
   * [0] - smallest distance of hypothetical 2g vertex from o-Ps->3g vertex
   * [1] - second smallest distance of hypothetical 2g vertex from o-Ps->3g vertex
   * [2] - largest distance of hypothetical 2g vertex from o-Ps->3g vertex
   * [3] - sum of the above distances
   *
   */ 
  Distances calcLORdistances(const JPetOpsEvent& event){

    Distances distances;

    // loop over all pairs among three hits
    double d_sum = 0.;
    for(int i=0;i<3;++i){
      for(int j=i+1;j<3;++j){
        TVector3 vtx_2g = EventCategorizerTools::calculateAnnihilationPoint(event.getHits().at(i),
                                                                            event.getHits().at(j));        
        double d = (vtx_2g - event.getAnnihilationPoint()).Mag();
        distances.push_back(d);
        d_sum += d;
      }
    }

    // sort the distances
    std::sort(distances.begin(), distances.end());

    // include the sum for convenience
    distances.push_back(d_sum);
       
    return distances;
  }

  /**
   * @brief Calculates discrepancies between inter-hit distance and hyopothetical TOF
   * for all hit pairs in a 3-hit event.
   *
   * Returns a vector of 3 values of |d-c*dt| [cm] in ascending order.
   *
   */ 
  ScatterTests calcScatterTests(const JPetOpsEvent& event){
    ScatterTests dvts;

    for(int i=0;i<3;++i){
      for(int j=i+1;j<3;++j){
        
        double d = (event.getHits().at(i).getPos() - event.getHits().at(j).getPos()).Mag();
        double dt = fabs(event.getHits().at(i).getTime() - event.getHits().at(j).getTime()) / 1000.;
        double dvt = fabs(d - dt*kSpeedOfLight);
        dvts.push_back(dvt);
      }
    }
    std::sort(dvts.begin(), dvts.end());

    return dvts;
  }
  
};

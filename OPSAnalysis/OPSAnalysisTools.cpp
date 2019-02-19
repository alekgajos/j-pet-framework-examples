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

namespace ops_analysis_tools{

  /**
     @brief Identify if the hit could correspond to annihilation photon,
     prompt photon or none of them.
     
     @param tot_cuts vector of 4 real numbers describing TOT cut boundaries:
     tot_cuts[0] - lower annihilation photon TOT cut
     tot_cuts[1] - upper annihilation photon TOT cut
     tot_cuts[2] - lower prompt photon TOT cut
     tot_cuts[3] - upper prompt photon TOT cut
  */
  HitCandidateType identifyHitType(const JPetHit& hit, std::vector<double>& tot_cuts){
    double tot = hit.getEnergy();
    // check annihilation TOT cuts
    if( tot > tot_cuts[0] && tot < tot_cuts[1] ){
      return HitCandidateType::Annihilation;
    }
    if( tot > tot_cuts[2] && tot < tot_cuts[3] ){
      return HitCandidateType::Prompt;
    }
    return HitCandidateType::None;
  }
  
};

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
 *  @file OPSAnalysisTools.h
 */

#ifndef OPSANALYSISTOOLS_H
#define OPSANALYSISTOOLS_H

#include <JPetHit/JPetHit.h>
#include "JPetOpsEvent.h"
#include <vector>
#include <tuple>
#include <utility>
#include <TVector3.h>

/**
   @brief Tools for identification and analysis of o-Ps->3gamma events
 */
namespace ops_analysis_tools{

  using Kinematics3g = std::tuple<std::vector<double>, std::vector<double>, std::vector<TVector3>>;
  using SmallestAngles = std::pair<double, double>;
  using Operators = std::vector<double>;
  using Distances = std::vector<double>;
  using ScatterTests = std::vector<double>;
  
  enum HitCandidateType{
    None,
    Annihilation,
    Prompt,
  };

  HitCandidateType identifyHitType(const JPetHit& hit, std::vector<double>& tot_cuts);

  SmallestAngles calcSmallestAnglesXY(const JPetOpsEvent& event);
  //  std::pair<double, double> calcSmallestAngles(const JPetOpsEvent& event);
  Kinematics3g calcKinematics(const JPetOpsEvent& event);
  Operators calcOperators(const JPetOpsEvent& event, Kinematics3g kinematics);
  Distances calcLORdistances(const JPetOpsEvent& event);
  ScatterTests calcScatterTests(const JPetOpsEvent& event);
};

#endif /* OPSANALYSISTOOLS_H */

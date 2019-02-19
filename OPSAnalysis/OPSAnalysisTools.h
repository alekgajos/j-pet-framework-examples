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
#include <vector>

/**
   @brief Tools for identification and analysis of o-Ps->3gamma events
 */
namespace ops_analysis_tools{

  enum HitCandidateType{
    None,
    Annihilation,
    Prompt,
  };

  HitCandidateType identifyHitType(const JPetHit& hit, std::vector<double>& tot_cuts);
    
};

#endif /* OPSANALYSISTOOLS_H */

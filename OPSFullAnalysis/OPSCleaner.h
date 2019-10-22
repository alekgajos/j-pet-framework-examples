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
 *  @file OPSCleaner.h
 */

#ifndef OPSCLEANER_H 
#define OPSCLEANER_H 

#include <vector>
#include <map>
#include <string>
#include <JPetUserTask/JPetUserTask.h>
#include <JPetHit/JPetHit.h>
#include <JPetEvent/JPetEvent.h>
#include "JPetOpsEvent.h"
#include "OPSAnalysisTools.h"

#ifdef __CINT__
#	define override
#endif

struct EventInfo{
  bool event_odd;
  ops_analysis_tools::SmallestAngles angles_xy;
  ops_analysis_tools::SmallestAngles angles_3d;
  ops_analysis_tools::Kinematics3g kinematics;
  ops_analysis_tools::Operators operators;
  ops_analysis_tools::Distances distances;
  ops_analysis_tools::ScatterTests scatter_tests;
};

class OPSCleaner : public JPetUserTask{
public:
  OPSCleaner(const char * name);
  virtual ~OPSCleaner(){}
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;
protected:
  void bookHisto(TH1* h);
  TH1F* getHisto1D(std::string name, int step);
  TH2F* getHisto2D(std::string name, int step);
  void fillHistos(const JPetOpsEvent& event, EventInfo evt_info, int step);
  
  const std::string fAngleSumCutKey = "OPSCleaner_angles_sum_cut_float";
  float fAngleSumCut;
  
  const double kSpeedOfLight = 29.9792458; // cm  / ns
  const int kNsteps = 10;
  
  int fEventCouters[10];

};
#endif /*  !OPSCLEANER_H */











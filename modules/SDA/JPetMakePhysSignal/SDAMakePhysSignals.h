/**
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file SDAMakePhysSignals.h
 *  @brief Dummy producer of JPetPhysSignal objects
 *  Reads a TTree of RecoSignals and transforms them into JPetPhysSignal objects
 *  PhysSignals have charge value of their RecoSignals set as fPhe
 */

#ifndef _JPETANALYSISMODULE_SDAMAKEPHYSSIGNALS_H_
#define _JPETANALYSISMODULE_SDAMAKEPHYSSIGNALS_H_

#include "JPetUserTask/JPetUserTask.h"
#include <TCanvas.h>

class SDAMakePhysSignals : public JPetUserTask {
public:
  SDAMakePhysSignals(const char *name);
  virtual ~SDAMakePhysSignals();
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;

private:
};

#endif

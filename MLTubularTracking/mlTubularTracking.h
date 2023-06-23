// Copyright (c) Fraunhofer MEVIS, Germany. All rights reserved.
// -----------------------------------------------------------------------
// 
// Copyright (c) 2001-2022, Fraunhofer MEVIS, Bremen, Germany
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Fraunhofer MEVIS nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY FRAUNHOFER MEVIS ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL FRAUNHOFER MEVIS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//----------------------------------------------------------------------------------

#ifndef MLTUBULARTRACKING_H
#define MLTUBULARTRACKING_H

#include "MLTubularTrackingSystem.h"
#include "VesselTemplateTracker.h"
#include "VesselTemplate.h"
#include "VesselGraph.h"

#include <mlGeneralGaussFunctions.h>

// Newmat includes
#include <ThirdPartyWarningsDisable.h>
#if defined(MLAB_CMAKE_BUILDSYSTEM)
#include "newmat/newmatap.h"
#else
#include "newmatap.h"
#endif
#include <ThirdPartyWarningsRestore.h>

#include <mlModuleIncludes.h>
#include <mlXMarkerList.h>
#include <mlVirtualVolume.h>
#include <mlBitImage.h>
#include <mlGraph.h>

ML_START_NAMESPACE

//! Tracks tubular structures
class MLTUBULARTRACKING_EXPORT TubularTracking : public Module {
public:

  TubularTracking (void);
  virtual ~TubularTracking (void) override;

  //! Handle field changes of the field \c field.
  virtual void handleNotification (Field *field) override;

  //! The mask at input 1 is optional. If there is no mask, this function will
  //! attach a dummy input to input 1.
  virtual INPUT_HANDLE handleInput(int inIndex, INPUT_STATE /*state*/) const override;

private:

  //! Implements interface for the runtime type system of the ML.
  ML_MODULE_CLASS_HEADER(TubularTracking)

  //! Output base field of type MLVesselGraph
  BaseField* _outputGraphFld;
  Graph   _outputGraph;

  //! Minimum vessel radius.
  DoubleField *_minRadiusFld;

  //! Maximum vessel radius.
  DoubleField *_maxRadiusFld;

  //! Initial vessel radius.
  DoubleField *_initRadiusFld;
    
  //! Use initial vessel radius.
  BoolField *_useInitRadiusFld;
    
  //! Maximum angle relative to the previous tracking step (in degrees).
  DoubleField *_maxAngleFld;

  //! Maximum number of search angles.
  IntField *_nbrSearchAnglesFld;

  //! Score threshold that determines when tracking should be terminated.
  DoubleField *_pruningThresholdFld;

  //! Score threshold that determines when tracking should be terminated.
  DoubleField *_terminationThresholdFld;

  //! Maximum number of tracking steps.
  IntField *_maxNbrStepsFld;

  //! Maximum tracking length in millimeters.
  DoubleField *_maxLengthFld;

  //! How deep should the search tree be?.
  IntField *_searchDepthFld;

  //! Minimum distance between two branchings in millimeters.
  DoubleField *_minBranchingDistanceFld;

  //! Tracking step length. The step length is calculated
  //! by multiplying this factor with the current vessel radius.
  DoubleField *_stepLengthFactorFld;

  //! Size factor for the (Gaussian) weight window used for fitting
  //! the vessel model.
  DoubleField *_windowSizeFactorFld;

  //! Base field containing an XMarkerList with two initialization.
  //! points or a point and a vector.
  BaseField *_inputInitPointsFld;

  //! Single hypothesis or multiple hypothesis tracking.
  BoolField*  _useMultipleHypothesesTrackingFld;

  //! Toggle field for (dis-)allowing branchings in the tracking process.
  BoolField*  _allowBranchingFld;

  //! Toggle field for growing bidirectional instead of in just the direction 
  //! given by the user.
  BoolField *_growBidirectionalFld;

  //! Toggle field for enabling/disabling tracking termination after a certain number of steps.
  BoolField* _toggleMaxStepsFld;

  //! Toggle field for enabling/disabling tracking termination after a certain length.
  BoolField* _toggleMaxLengthFld;

  //! Base field containing an XMarkerList with tracked points.
  BaseField *_outputTrackedPointsFld;

  //! List containing initial points.
  XMarkerList _initPointsXMarkerList;

  //! List containing the tracked points.
  XMarkerList _trackedPointsXMarkerList;

  //! Button for starting tracking.
  NotifyField*  _updateButtonFld;

  //! Button for clearing previous tracking result.
  NotifyField*  _clearButtonFld;

  //! Virtual volume for input data image.
  TVirtualVolume<float> *_virtInputVolume;
 
  //! Virtual volume for input mask image.
  TVirtualVolume<MLuint8> *_virtMaskVolume;

  //! Segmented binary image.
  BitImage _binaryImg; 

  //! Variable that indicates if output data has changed. Used to suppress unnecessary touch().
  bool _outputChanged;

  bool _getInputXmarkers();

  //! Get block of image data
  void _getImageData(TSubImage<float> &imData, TSubImage<MLuint8> &maskData, const MedicalImageProperties &imgProps, double worldCenterPoint[3], int blockSize,  DIM trackDim) const;

  //! A vector of pairs containing an XMarker and a boost graph with a vessel
  //! tree tracked starting from the XMarker position.
  std::vector<std::pair<XMarker, boostGraph> > _trackedBoostGraphs;
  
  //! ImageVector of pointers to the algorithm parameter Fields. This
  //! vector is used to check whether a change to *any* of the parameters
  //! influencing the algorithm has been made. In such case is the stored
  //! tracking result in _trackedBoostGraphs invalid..
  std::vector<Field*> _parameterFlds;

  Vector3 estimateOrientation(const Vector3& pos, const MedicalImageProperties& imgProps, TVirtualVolume<float>& virtDataVol, 
                              double radius, DIM trackDim);
  
};

ML_END_NAMESPACE

#endif // __mlTubularTracking_H

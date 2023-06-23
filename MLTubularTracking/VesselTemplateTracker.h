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
#ifndef __VesselTemplateTracker_H
#define __VesselTemplateTracker_H

#include "VesselTracker.h"
#include "ImPatch.h"
#include "VesselTemplate.h"


class VesselTemplateTracker : public VesselTracker
{
public:
  // When the data pointer in maskData is NULL, no mask is supplied.
  VesselTemplateTracker(const ML_NAMESPACE::TSubImage<float>& imData, const ML_NAMESPACE::TSubImage<MLuint8>& maskData,
                        const ML_NAMESPACE::MedicalImageProperties& imProps, DIM TRACKDIM);

  virtual ~VesselTemplateTracker() override {};

  //! Mark voxels that have been visited, i.e, voxelize a vessel model.
  virtual void setVisitedVoxels(const VesselData& vd, ML_NAMESPACE::BitImage& binaryImg) override;

  //! Sets the width factor of the Gaussian weight window. sigma = vesselRadius*windowSizeFactor.
  void setWindowSizeFactor(double windowSizeFactor) {_windowSizeFactor = windowSizeFactor;}

private:

  //! Levenberg-Marquardt optimization of the parameters.
  virtual void _fitModel(VesselData& vd, const VesselData& vd_prev) override;

  //! Calculate score of a vessel template, i.e., how well it fits the image
  virtual void _calcScore(VesselData& vd) override;

  //! Determines a suitable template size in voxels. The half size is returned for
  //! each dimension (x,y,z) is returned. That is, if n is the half side, 2*n+1 is
  //! the full size.
  virtual void _getPatchSize(int voxelHalfSize[3], double radius) override;

  //! Linear least squares fit for contrast and mean parameters.
  void _lsFit(VesselTemplate& vt, const ImPatch& T, const ImPatch& I, const ImPatch& W);

  void _calcScore(VesselTemplate& vt, const ImPatch& T, const ImPatch& I, const ImPatch& W);

  //! A Gaussian weightfunction.
  void _weightFunction(ImPatch& W, const VesselData& vd);

  double _getWeightThreshold();

  //! Returns a suitable size of the window function.
  double _getWindowSigma(double radius);

  //! Get an image patch adapted to the vessel data in VesselData.
  void _getImPatch(ImPatch& I, const VesselData& vd);



private:

    VesselTemplate _vt;

    //! Width factor of the Gaussian weight window. sigma = vesselRadius*windowSizeFactor.
    double  _windowSizeFactor;

    // ImPatches for template, weight and image data
    ImPatch _T;
    ImPatch _W;
    ImPatch _I;
    ImPatch _trial_T;

    // Vector holding residuals, i.e., remaining image noise
    // when the fitted template has been substracted from the
    // image patch.
    std::vector<double> _wresiduals;
};

#endif

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
#ifndef __VesselTemplate_H
#define __VesselTemplate_H

// Get rid of the min() / max () macro definitions in Visual Studio
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#include "VesselData.h"
#include "VesselProfile.h"
#include "ImPatch.h"
#include <vector>
#include <assert.h>

#define MAXSIZE 100

class VesselTemplate : public VesselData, public VesselProfile
{
public:
  //! Standard constructor.
  VesselTemplate(void);

  //! Standard destructor.
  virtual ~VesselTemplate() {};

  //! Constructor for converting a VesselData object into a
  //! VesselTemplate object.
  VesselTemplate(const VesselData& vd);

  //! Get template value for the input coordinate (x,y,z).
  //! The weight W acts as a mask here, the template is calculated for all non-zero W.
  void getValue(ImPatch& T,   const ImPatch& W);
  void getValue2D(ImPatch& T, const ImPatch& W);
  void getValue3D(ImPatch& T, const ImPatch& W);

  //! Get template value and derivatives 2D.
  void  getDerivatives(std::vector<double>& dT_dr,
                       std::vector<double>& dT_du,
                       std::vector<double>& dT_dtheta,
                       const std::vector<double>& T,
                       double xstart, double dx, int xsize,
                       double ystart, double dy, int ysize);

  //! Get template value and derivatives 3D.
  void  getDerivatives(std::vector<double>& dT_dr,
    std::vector<double>& dT_du1,
    std::vector<double>& dT_du2,
    std::vector<double>& dT_dtheta,
    std::vector<double>& dT_dphi,
    const std::vector<double>& T,
    double xstart, double dx, int xsize,
    double ystart, double dy, int ysize,
    double zstart, double dz, int zsize);

  //! Set contrast. Use protected member in VesselData class for this.
  void setContrast(double contrast) {_p1 = contrast;}
  //! Get contrast.
  double getContrast(void) const {return _p1;}

  //! Set mean. Use protected member in VesselData class for this.
  void setMean(double mean) {_p2 = mean;}
  //! Get mean.
  double getMean(void) const {return _p2;}

};
#endif

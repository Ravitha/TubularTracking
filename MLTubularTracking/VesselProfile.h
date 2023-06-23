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
#ifndef __VesselProfile_H
#define __VesselProfile_H

#include <math.h>
#include "VesselData.h"

// Forward declaration of VesselData class.
class VesselData;

#if 0
// Gaussian profile
class VesselProfile {
public:

  // Constructor. Default is a radius=1 and beta=8
  VesselProfile(void) : _1_div_r2(1.0), _2_log2_div_r3(2*log((double)2.0)), log2_div_r2(log((double)2.0)), _lowLimit(0.0), _highLimit(10.0) {}


  inline double getProfile(double d2) const
  {
    return pow((double)2.0,-d2*_1_div_r2);
  }


  inline double getProfile_dr(double d2, double p) const  {
    return _2_log2_div_r3*d2*p;
  }


  inline double getProfile_dx2(double /*d2*/, double p) const
  {
    return -log2_div_r2*p;
  }


  // For d2<lowLimit the vessel profile can be approximated with 1
  // For d2>highLimit the vessel profile can be approximated with 0
  void  getProfileLimits(double& lowLimit, double& highLimit) const {lowLimit=_lowLimit; highLimit=_highLimit;}


  void updateProfileParameters(const VesselData* const vd)
  {
    const double r = vd->getRadius();

    _1_div_r2 = 1.0/(r*r);
    _2_log2_div_r3 = 2.0*log((double)2.0)/(r*r*r);
    log2_div_r2 = log((double)2.0)/(r*r);


    // Inverse formula:
    // d2 = -r^2/log(2)*log(y)
    const double yLow = 0.03f;
    const double yHigh = 1-yLow;
    _highLimit = -r*r/log((double)2.0)*log(yLow);
    _lowLimit = -r*r/log((double)2.0)*log(yHigh);

  }

private:

  //! Internal variables that are used repeatedly.
  double _1_div_r2;
  double _2_log2_div_r3;
  double log2_div_r2;

  //! For d2<_lowLimit the profile value >0.99 and can be approximated with 1.
  double _lowLimit;

  //! For d2>_highLimit the profile value <0.01 and can be approximated with 0.
  double _highLimit;
};


#endif




class VesselProfile {
public:

  // Constructor. Default is a radius=1 and beta=8
  VesselProfile(void) : _rpowbeta(1), _beta_div_rpowbeta_div2(4), _beta_div_radius(8), _lowLimit(0), _highLimit(10) {}


  inline double getProfile(double d2) const
  {
    double dpowbeta_2 = d2*d2;
    dpowbeta_2 *= dpowbeta_2;                 // For beta = 8
    return _rpowbeta/(dpowbeta_2 + _rpowbeta);
  }


  inline double getProfile_dr(double x2, double p) const  {
    return _beta_div_radius*p*(1-p);
  }


  inline double getProfile_dx2(double x2, double p) const
  {
    return -p*p*_beta_div_rpowbeta_div2*x2*x2*x2;
  }


  // For d2<lowLimit the vessel profile can be approximated with 1
  // For d2>highLimit the vessel profile can be approximated with 0
  void  getProfileLimits(double& lowLimit, double& highLimit) const {lowLimit=_lowLimit; highLimit=_highLimit;}


  void updateProfileParameters(const VesselData* const vd)
  {
    const double radius = vd->getRadius();
    const double beta = 8.0;
    _rpowbeta = pow(radius,beta);
    _beta_div_rpowbeta_div2 = beta/_rpowbeta/2;
    _beta_div_radius = beta/radius;

    // For which d2 is r^beta/(d2^(beta/2) + r^beta) equal y
    // Inverse formula:
    // d2 = r^2 * ((1-y)/y)^(2/beta)
    const double yLow = 0.03f;
    const double yHigh = 1-yLow;
    _highLimit = radius*radius * pow((1-yLow )/yLow ,2.0/beta);
    _lowLimit  = radius*radius * pow((1-yHigh)/yHigh,2.0/beta);

  }

private:

  //! Internal variables that are used repeatedly.
  double _rpowbeta;
  double _beta_div_rpowbeta_div2;
  double _beta_div_radius;

  //! For d2<_lowLimit the profile value >0.99 and can be approximated with 1.
  double _lowLimit;

  //! For d2>_highLimit the profile value <0.01 and can be approximated with 0.
  double _highLimit;
};


#endif


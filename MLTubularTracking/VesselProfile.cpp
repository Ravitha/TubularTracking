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
#include "VesselProfile.h"

//! Calculate profile value according to the formula
//! p(d2) = r^beta/(d2^(beta/2) + r^beta).
double VesselProfile::getProfile(double d2) const
{
    double dpowbeta_2 = d2*d2;
    dpowbeta_2 *= dpowbeta_2;                 // For beta = 8
    return _rpowbeta/(dpowbeta_2 + _rpowbeta);
}


//! Derivative of profile with respect to radius parameter.
//! Input is squared distance x2 and profile value p(x2).
double VesselProfile::getProfile_dr(double x2, double p) const
{
  return _beta_div_radius*p*(1-p);
}

//! Derivative of profile with respect to radius parameter.
//! Input is squared distance x2 and profile value p(x2).
double VesselProfile::getProfile_dx2(double x2, double p) const
{
  return -p*p*_beta_div_rpowbeta_div2*x2*x2*x2;
}


void VesselProfile::updateProfileParameters(const VesselData* const vd) 
{
  const double radius = vd->getRadius();
  const int beta = 8;
  _rpowbeta = pow(radius,beta);
  _beta_div_rpowbeta_div2 = beta/_rpowbeta/2;
  _beta_div_radius = beta/radius;

  // For which d2 is r^beta/(d2^(beta/2) + r^beta) equal y
  // Inverse formula:
  // d2 = r^2 * ((1-y)/y)^(2/beta)
  const double yLow = 0.03f;
  const double yHigh = 1-yLow;
  _highLimit = radius*radius * pow((1-yLow)/yLow,2.0/(float)beta);
  _lowLimit = radius*radius * pow((1-yHigh)/yHigh,2.0/(float)beta);
 
}



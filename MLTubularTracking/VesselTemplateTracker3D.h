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
#ifndef __VesselTemplateTracker3D_H
#define __VesselTemplateTracker3D_H

#include "VesselTracker.h"
#include "VesselTemplate.h"

template<class TM>    // TODO: Another solution should be possible
class VesselTemplateTracker3D : public VesselTracker
{
public:
  VesselTemplateTracker3D(const ML_NAMESPACE::TSubImage<double>& imData, 
    const ML_NAMESPACE::MedicalImageProperties& imProps, DIM TRACKDIM) : VesselTracker(imData, imProps, TRACKDIM) {};

private:

  //! Levenberg-Marquardt optimization of the parameters.
  virtual void _fitModel(VesselData& vd, const VesselData& vd_prev) const;

  //! Calculate score of a vessel template, i.e., how well it fits the image
  virtual void _calcScore(VesselData& vd) const;

  //! Linear least squares fit for contrast and mean parameters.
  void _lsFit(TM& vt, const std::valarray<double>& T, const std::valarray<double>& I, const std::valarray<double>& W2) const;

  void _calcScore(TM& vt, const ImPatch& I, const std::valarray<double>& vesselTemplateValues) const;

  //! A Gaussian weightfunction.
  void _weightFunction(std::valarray<double>& W, const VesselData& vd, const ImPatch& I, bool square) const;

  //! Returns a suitable size of the window function.
  double _getWindowFWHM(double radius) const;

  //! Determines a suitable template size in voxels. The half size is returned for
  //! each dimension (x,y,z) is returned. That is, if n is the half side, 2*n+1 is
  //! the full size.
  void _getTemplateSize(int voxelHalfSize[3], double radius) const;

  void _getImPatch(ImPatch& I, const VesselData& vt) const;

};


// --------------------------------------------------------------------------------------------------
//! Levenberg-Marquardt method for solving the non-linear least squares problem
//! min || W(x)[I(x) - k*T(x;r,v,x0) - m] ||^2
//! where I(x) is an image patch, T(x;r,v,x0) is the tubular template and W(x) is a
//! weight function or a window centered (approximately) above x0.
//!
//! The algorithm iterates the following two steps:
//! 1) Levenberg-Marquardt method for finding the non-linear parameters r (radius)
//!    v (direction) and x0 (center point). The linear parametes k and m are kept fixed.
//! 2) Ordinary least squares for finding the linear parameters k (contrast) and 
//!    m (mean voxel intensity). The non-linear parameters r,v and x0 are kept fixed.
// --------------------------------------------------------------------------------------------------
template<class TM> void VesselTemplateTracker3D<TM>::_fitModel(VesselData& vd, const VesselData& vd_prev) const
{ 
  // Create a VesselTemplate object out of the VesselData object
  TM vt(vd);

   // Get position of the previous template
  double x0_prev,y0_prev,z0_prev;
  vd_prev.getCenterPoint(x0_prev,y0_prev,z0_prev);
  
  // Get initial position of the current template
  double x0,y0,z0;
  vt.getCenterPoint(x0,y0,z0);

  // Calculate step length between the templates. This step length will be preserved.
  double stepLength = sqrt((x0_prev-x0)*(x0_prev-x0) + (y0_prev-y0)*(y0_prev-y0) + (z0_prev-z0)*(z0_prev-z0));

  // Get image data
  ImPatch I; 
  _getImPatch(I, vt);
  const unsigned int nbrOfVoxels = I.voxelValues.size();

  // Calculate sum of voxel values and the weight function W2 containing squared weights. 
  // These quantities areused repeatedly below. The image is weighted with W already here.
  //unsigned int i = 0;
  //double sW2 = W2.sum();
  //double sW2I = (W2*I.voxelValues).sum();

  // Get template
  std::valarray<double> vesselTemplateValues(nbrOfVoxels);       
  if (_TRACKDIM == _2D) {
    vt.getValue(vesselTemplateValues, I.x, I.y);
  }
  else if (_TRACKDIM == _3D) {
    vt.getValue(vesselTemplateValues, I.x, I.y, I.z);
  }

  // Outer optimization loop that iterates
  // 1) Levenberg-Marquardt method for finding the non-linear parameters r (radius)
  //    v (direction) and x0 (center point). The linear parametes k and m are kept fixed.
  // 2) Ordinary least squares for finding the linear parameters k (contrast) and 
  //    m (mean voxel intensity). The non-linear parameters r,v and x0 are kept fixed.
  // 
  // We do a maximum of 10 iterations and as long as we find Levenberg-Marquardt
  // steps that improve the fit.
  int outerIts = 0;
  bool foundLMStep = true;
  while( (foundLMStep) && (outerIts < 10))    // TODO: Introduce convergence criterion
  {
    // Increment iteration counter
    ++outerIts;

    // Get appropriate weight function.
    // TODO: Try moving the weight calculation outside the loop.
    std::valarray<double> W(nbrOfVoxels);
    _weightFunction(W, vd, I, false);            

    vt.DEBUG_W = W;
    vt.DEBUG_I = I.voxelValues;

    // Linear fit to find contrast and mean
    _lsFit(vt,vesselTemplateValues,I.voxelValues,W);

    // Get template parameters
    double contrast = vt.getContrast();
    double mean = vt.getMean();
    double radius = vt.getRadius();
    vt.getCenterPoint(x0,y0,z0);
    double u1x, u1y, u1z;
    double u2x, u2y, u2z;
    vt.getOrthoDirection(u1x,u1y,u1z,u2x,u2y,u2z);

    // Calculate residuals and loss function.
    std::valarray<double> wResiduals = W*(contrast*vesselTemplateValues + mean - I.voxelValues);     
    double lossFunction = (wResiduals*wResiduals).sum();

    // Levenberg-Marquardt iterations to find optimal radius (r), direction (v), and center point (x0).
    // Build matrices J'*W2*J and J'*W2*r where J is the Jacobian matrix, i.e, a matrix with the
    // partial derivatives, and r is the residuals (i.e, residuals). In this case
    // J = [dT/dr, dT/du1, dT/du2, dT/dtheta, dT/dphi] stretched to vectors.
    // u1 and u2 are two orthogonal direction to the template direction v. These are used to determine
    // the center point x0 (it makes no sense to translate the center point along the direction v).
    // Theta and phi are the angles of v in a spherical coordinate system.

    const unsigned int NBRFITPARAMS = (_TRACKDIM == _3D)? 5 : 3;
    SymmetricMatrix JW2J(NBRFITPARAMS); JW2J = 0.0;
    ColumnVector JW2r(NBRFITPARAMS); JW2r = 0.0;
    
    if(_TRACKDIM == _2D) {
      std::valarray<double> wj1,wj2,wj3;
      vt.getDerivatives(wj1, wj2, wj3, I.x, I.y);
      wj1 *= W; wj2 *= W; wj3 *= W;
      JW2J(1,1) = (wj1*wj1).sum(); JW2J(1,2) = (wj1*wj2).sum(); JW2J(1,3) = (wj1*wj3).sum();
      JW2J(2,2) = (wj2*wj2).sum(); JW2J(2,3) = (wj2*wj3).sum();
      JW2J(3,3) = (wj3*wj3).sum();
      JW2r(1) = (wj1*wResiduals).sum(); JW2r(2) = (wj2*wResiduals).sum(); JW2r(3) = (wj3*wResiduals).sum();   
    }
    else if (_TRACKDIM == _3D) {
      std::valarray<double> wj1,wj2,wj3,wj4,wj5;
      vt.getDerivatives(wj1, wj2, wj3, wj4, wj5, I.x, I.y, I.z);
      wj1 *= W; wj2 *= W; wj3 *= W; wj4 *= W; wj5 *= W;
      JW2J(1,1) = (wj1*wj1).sum(); JW2J(1,2) = (wj1*wj2).sum(); JW2J(1,3) = (wj1*wj3).sum(); JW2J(1,4) = (wj1*wj4).sum(); JW2J(1,5) = (wj1*wj5).sum();
      JW2J(2,2) = (wj2*wj2).sum(); JW2J(2,3) = (wj2*wj3).sum(); JW2J(2,4) = (wj2*wj4).sum(); JW2J(2,5) = (wj2*wj5).sum();
      JW2J(3,3) = (wj3*wj3).sum(); JW2J(3,4) = (wj3*wj4).sum(); JW2J(3,5) = (wj3*wj5).sum(); 
      JW2J(4,4) = (wj4*wj4).sum(); JW2J(4,5) = (wj4*wj5).sum();
      JW2J(5,5) = (wj5*wj5).sum();
      JW2r(1) = (wj1*wResiduals).sum(); JW2r(2) = (wj2*wResiduals).sum(); JW2r(3) = (wj3*wResiduals).sum(); 
      JW2r(4) = (wj4*wResiduals).sum(); JW2r(5) = (wj5*wResiduals).sum();
    }
    JW2J *= contrast*contrast;
    JW2r *= contrast;  

    // Find maximum of the diagonal of JW2J. This is needed as we below (in the Levenberg-Marquardt algorithm)
    // will add a small number to the diagonal (as a regularization), and we need to find out how small "small" is.
    // Specifically a first guess of "small" will be 1e-3*max(diag(J'*W2*J))
    double maxVal = -1;
    unsigned int i = 0;
    for(i=1; i<=NBRFITPARAMS; ++i) {
      if(JW2J(i,i) > maxVal) {maxVal = JW2J(i,i);}
    }
    // Check if maxVal is less or equal to zero. This means problems
    // as the diagonal of JW2J should contain positive elements.
    if(maxVal<=0) {break;}

    // Value to be added to the diagonal
    double mu = 1e-3*maxVal;               
    // Factor to modify mu with in case Levenberg-Marquardt is rejected/accepted
    double ny = 2.0;

    // Perform Levenberg-Marquardt iterations
    foundLMStep = false;
    int lmIts = 0;
    while ((!foundLMStep) && (lmIts < 10)) {

      // Increment iteration counter
      ++lmIts;

      // Solve the equation system (JW2J + mu)*hlm = JW2I
      SymmetricMatrix A = JW2J;
      ColumnVector hlm;
      for(i=1; i<=NBRFITPARAMS; ++i) {A(i,i) += mu;}
      try {
         hlm = A.i() * (-JW2r);
      }
      catch(BaseException) {
        ML_PRINT_WARNING("VesselTracker::_detectBranching", ML_PROGRAMMING_ERROR,  BaseException::what());
        break;
      }

      // Calculate norm of hlm. If the norm is too small we don't bother to go on
      double lmNorm = 0;
      for(i=1; i<=NBRFITPARAMS; ++i) {lmNorm += hlm(i)*hlm(i);}
      if(sqrt(lmNorm)<1e-3) {break;}
      
      // TODO: If norm is too big, don't bother to test

      // Create a new test vessel template
      TM t_vt = vt;

      // Update radius of the test template. 
      // The radius is the first parameter in both 2D and 3D.
      double t_radius = radius + hlm(1);
      if(t_radius < getMinRadius()) { t_radius = getMinRadius();}
      if(t_radius > getMaxRadius()) { t_radius = getMaxRadius();}
      t_vt.setRadius(t_radius);

      // Update center point by moving in the directions u1 and u2 that are orthogonal to the current
      // template direction v. There is a not-so-elegant step length solution here: for the user-defined
      // start template the step length is 0 because no previous template exists. For this special case
      // there is no restriction on the template movement.
      if(_TRACKDIM == _2D) {
        const double sx = x0 + hlm(2)*u1x - x0_prev;
        const double sy = y0 + hlm(2)*u1y - y0_prev;
        const double norm = sqrt(sx*sx + sy*sy);
        const double t_x = (stepLength < 1e-8)? x0 + hlm(2)*u1x : x0_prev + sx/norm*stepLength;
        const double t_y = (stepLength < 1e-8)? y0 + hlm(2)*u1y : y0_prev + sy/norm*stepLength;
        t_vt.setCenterPoint(t_x,t_y);

        // Rotatate direction according to the angles dtheta and dphi.
        t_vt.rotate(hlm(3));
      }
      else if(_TRACKDIM == _3D) {
        const double sx = x0 + hlm(2)*u1x + hlm(3)*u2x - x0_prev;
        const double sy = y0 + hlm(2)*u1y + hlm(3)*u2y - y0_prev;
        const double sz = z0 + hlm(2)*u1z + hlm(3)*u2z - z0_prev;
        const double norm = sqrt(sx*sx + sy+sy + sz+sz);
        
        const double t_x = (stepLength < 1e-8)? x0 + hlm(2)*u1x + hlm(3)*u2x : x0_prev + sx/norm*stepLength;
        const double t_y = (stepLength < 1e-8)? y0 + hlm(2)*u1y + hlm(3)*u2y : y0_prev + sy/norm*stepLength;
        const double t_z = (stepLength < 1e-8)? z0 + hlm(2)*u1z + hlm(3)*u2z : z0_prev + sz/norm*stepLength;
        t_vt.setCenterPoint(t_x,t_y,t_z);

        // Rotatate direction according to the angles dtheta and dphi.
        t_vt.rotate(hlm(4),hlm(5));
      }

      // Get trial template
      std::valarray<double> t_vesselTemplateValues(nbrOfVoxels);       
      if (_TRACKDIM == _2D) {
        t_vt.getValue(t_vesselTemplateValues, I.x, I.y);
      }
      else if (_TRACKDIM == _3D) {
        t_vt.getValue(t_vesselTemplateValues, I.x, I.y, I.z);
      }

      // Calculate the loss function for the trial vessel template
      double t_lossFunction = 0;    
      for(i=0; i<nbrOfVoxels; ++i) {
        double wr = W[i]*(contrast*t_vesselTemplateValues[i] + mean - I.voxelValues[i]);  
        t_lossFunction += wr*wr;
      }

      // Compare the loss function with the one after the linear fit to see if
      // the adjustment of the non-linear parameters have brought a reduction
      // of the loss function. If so, we accept the adjustments.
      if(t_lossFunction < lossFunction)
      {
        // Set current template to the test template
        vt = t_vt;
        // Copy the template values so that we don't have to calculate them again below
        vesselTemplateValues = t_vesselTemplateValues;
        // Yes, we found a step.
        foundLMStep = true;
        // Let's go on to the linear step.
        break;   
      }
      // In case no reduction in the loss function we increse the regularisation
      // factor that is added to the diagonal an make a new try.
      else
      {
        mu = mu*ny;
        ny = 2*ny;
      }
    } // End while Levenberg-Marquardt iterations
  } // End while out loop over linear fit and non-linear Levenberg-Marquardt iterations

  // Finally, calculate the score for the tubular template
  _calcScore(vt,I,vesselTemplateValues);

  // Set fitted parameters.
  vd = vt;
}


// --------------------------------------------------------------------------------------------------
//! Weighted linear least squares to update k and m parameters.
// --------------------------------------------------------------------------------------------------
template<class TM> void VesselTemplateTracker3D<TM>::_lsFit(TM& vt, const std::valarray<double>& T, const std::valarray<double>& I, const std::valarray<double>& W) const
{
  // Get number of voxels. It is assumed that the vectors containing the template values T and the weights in W2
  // are of equal length, so this is not checked.
  const unsigned int nbrOfVoxels = I.size();  

  // The normal equation for solving this problem involves
  // only 2x2 matrices, therefore, we do this "by hand". The normal equations are
  // [m_opt, k_opt] = inv(X'*W2*X)*X'*W2*I, where X = [1 T] stretched to vectors and the image data I also
  // strectched to a vector. Below, the components in inv(X'*W2*X) and X'*W2*I are calculated separately
  // and then multiplied together. 
  double sW2=0, sW2T = 0,sW2T2 = 0, sW2TI = 0, sW2I = 0;
  for(unsigned int i=0; i<nbrOfVoxels; ++i) {
    const double wi = W[i]*I[i];
    const double wt = W[i]*T[i];
    sW2   += W[i]*W[i];
    sW2I  += W[i]*wi;   
    sW2T  += W[i]*wt;
    sW2T2 += wt*wt;   
    sW2TI += wi*wt;             
  }
  const double factor =  1/(sW2*sW2T2-sW2T*sW2T);

  // Note: We now have inv(X'*W2*X) = 1/factor [sW2T2, -sW2T; -sW2T, sW2]
  // and X'*W2*I =[sW2I; sW2TI]

  // Set fitted parameters m and k
  vt.setMean(factor*(sW2T2*sW2I - sW2T*sW2TI));
  vt.setContrast(factor*(sW2*sW2TI - sW2T*sW2I));
}

// --------------------------------------------------------------------------------------------------
//! Convenience _calcScore function that calls the overloaded _calcScore(...) above.
//! This function simply gets the image data and the template values and 
//! then calls the overloaded function.
// --------------------------------------------------------------------------------------------------
template<class TM> void VesselTemplateTracker3D<TM>::_calcScore(VesselData& vd) const
{
  // Create a VesselTemplate object out of the VesselData object
  TM vt(vd);

  // Get image data
  ImPatch I; 
  _getImPatch(I, vt);

  // Get number of voxels in patch
  unsigned int nbrOfVoxels = I.voxelValues.size();

  // Get template
  std::valarray<double> T(nbrOfVoxels);       
  if (_TRACKDIM == _2D) {
    vt.getValue(T, I.x, I.y);
  }
  else if (_TRACKDIM == _3D) {
    vt.getValue(T, I.x, I.y, I.z);
  }

  // Call overloaded _calcScore() function.
  _calcScore(vt,I,T);

  // Set score
  vd.setScore(vt.getScore());
}



// --------------------------------------------------------------------------------------------------
//! Calculate template score. This is based on the template contrast, i.e., how different from
//! zero the template contrast parameter is. The score equals the t-statistic
//!  t = contrast/std(contrast) or more generally c'beta/sqrt(var(c'beta))
//! where beta=[mean,contrast]' are the image model parameters 
//! Image = contrast*Template + mean + noise
//! and c is a so-called contrast vector, in our case c = [0,1]'.
//! The larger the score, the more different from zero the contrast is.
//! The variance in the denominator is calculated using this very nice formula:
//! var(c'beta) = sigma2*c'*inv(X'*W^2*X)*X'*W^4*X*inv(X'*W^2*X)*c.
//! X = [1 T] is the design matrix with a column of ones and the tubular template stretched to a column.
//! sigma2 is the noise variance which is estimated from the residuals r = I - X*beta. 
//! W is a diagonal matrix with a Gaussian weight window. This window focusses the score measurement to
//! the center of the tubular template and weights down areas far away from the center.
//! The estimate of sigma2 is weighted so that
//! sigma2 = sum(wi*ri^2)/sum(wi)
//! Note that the weighted estimate of the variance is ad-hoc, but it is a correct estimator.
// --------------------------------------------------------------------------------------------------
template<class TM> void VesselTemplateTracker3D<TM>::_calcScore(TM& vt, const ImPatch& I, const std::valarray<double>& vesselTemplateValues) const
{
  // Get number of voxels in patch
  const unsigned int nbrOfVoxels = I.voxelValues.size();

  // Get template parameters
  const double contrast = vt.getContrast();
  const double mean = vt.getMean();
  //double x0,y0,z0;
  //vt.getCenterPoint(x0,y0,z0);

  // The expression c'*inv(X'*W^2*X)*X'*W^4*X*inv(X'*W^2*X)*c was entered in
  // the Matlab symbolic toolbox, which gave the formula shown further down.
  // The components in this formula are calculated below.
  // sigma2 estimator
  double sigma2 = 0;

  // Get appropriate weight function
  std::valarray<double> W(nbrOfVoxels);
  _weightFunction(W, vt, I, false);

  // Sums resulting from the matrix multiplications
  double sW=0, sW2=0, sW2T=0, sW2T2=0, sW4=0, sW4T=0, sW4T2=0;
  for(unsigned int i=0; i<nbrOfVoxels; ++i) {

    // Pre-calculate powers of the weight
    double W2 = W[i]*W[i];
    double W4 = W2*W2;
    double W2T = W2*vesselTemplateValues[i];

    // Increment sigma2
    // TODO: Should be weighted with W2 here?
    double r = I.voxelValues[i] - contrast*vesselTemplateValues[i] - mean;
    sigma2 += W[i]*r*r;

    // Increment sums
    sW    += W[i];
    sW2   += W2;
    sW2T  += W2T;
    sW2T2 += W2T*vesselTemplateValues[i];           
    sW4   += W4;
    sW4T  += W4*vesselTemplateValues[i];
    sW4T2 += W2T*W2T;
  }
  sigma2 /= sW;


  // TODO: Should perhaps only calculate contrast/sigma as score!! 
  // The following expression for c'*inv(X'*W^2*X)*X'*W^4*X*inv(X'*W^2*X)*c was found using
  // the symbolic toolbox in Matlab.
  double tmp = sW2T2*sW2-sW2T*sW2T;
  double score = contrast/sqrt(sigma2*(sW2*sW2*sW4T2-2*sW2*sW2T*sW4T+sW2T*sW2T*sW4)/(tmp*tmp));

  // Finally, set the score
  vt.setScore(score);
}


// --------------------------------------------------------------------------------------------------
//! Calculates an anisotropic Gaussian weight function.  
//! TODO: For an isotropic kernel an outer product may be faster.
// --------------------------------------------------------------------------------------------------
template<class TM> void VesselTemplateTracker3D<TM>::_weightFunction(std::valarray<double>& W, const VesselData& vd, const ImPatch& I, bool square) const
{

  // Get suitable Full-Width at Half-Maximum of weight window.
  const double worldFWHM = _getWindowFWHM(vd.getRadius());
  
  // Converions factor 2.35 between sigma in a Gaussian weight function
  // and FWHM.
  const double worldSigma2 = (worldFWHM*worldFWHM)/(2.35*2.35);

  // The weight function is anisotropic and broader orthogonal to the
  // template direction (i.e. in the plane spanned by u1 and u2) 
  // than in the template direction given by v. Moreover, if the weight 
  // function is not going to be squared, we multiply sigma2 with 2 here 
  // to save multiplications when calculating the Gaussian exp(-x^2/(2*sigma2)).
  const double worldSigma2_u = (square)? worldSigma2   : 2*worldSigma2;
  const double worldSigma2_v = (square)? worldSigma2/4 : worldSigma2;

  // Get number of voxels and prepare valarray holding the weights.
  const unsigned int nbrOfVoxels = I.x.size();
  W.resize(nbrOfVoxels);
  
  // Get template center point.
  double x0,y0,z0;
  vd.getCenterPoint(x0,y0,z0);

  // Get template direction.
  double vx,vy,vz;
  vd.getDirection(vx, vy, vz);

  // Get orthogonal directions.
  double u1x,u1y,u1z,u2x,u2y,u2z;
  vd.getOrthoDirection(u1x,u1y,u1z,u2x,u2y,u2z);

  // Transform voxel coordinates into the coordinate system spanned
  // by the vessel template direction and the orthogonal direction.
  // Then calculate an anisotropic Gaussian function in this coordinate system.
  if(_TRACKDIM == _2D) {
    for(unsigned int i=0; i<nbrOfVoxels; ++i) {
      const double cx = I.x[i]-x0;            // Centered x
      const double cy = I.y[i]-y0;            // Centered y
      const double tv = cx*vx + cy*vy;        // v - coordinate
      const double tu = cx*u1x + cy*u1y;      // u - coordinate
      W[i] = exp(-tv*tv/worldSigma2_v - tu*tu/worldSigma2_u);
    }
  }
  else if (_TRACKDIM == _3D) {
    for(unsigned int i=0; i<nbrOfVoxels; ++i) {
      const double cx = I.x[i]-x0;            // Centered x
      const double cy = I.y[i]-y0;            // Centered y
      const double cz = I.z[i]-z0;            // Centered y
      const double tv = cx*vx + cy*vy + cz*vz;          // v - coordinate
      const double tu1 = cx*u1x + cy*u1y + cz*u1z;      // u1 - coordinate
      const double tu2 = cx*u2x + cy*u2y + cz*u2z;      // u1 - coordinate
      W[i] = exp(-tv*tv/worldSigma2_v - (tu1*tu1+tu2*tu2)/worldSigma2_u); 
    }
  }
}


// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
template<class TM> inline  double VesselTemplateTracker3D<TM>::_getWindowFWHM(double worldRadius) const
{
  return 5*worldRadius;   // One diameter on each side: d+d+d = 6*r (Subject to tuning)
}


// --------------------------------------------------------------------------------------------------
//! Sets the size of the vessel template as a function of the vessel radius
// --------------------------------------------------------------------------------------------------
template<class TM> inline void VesselTemplateTracker3D<TM>::_getTemplateSize(int voxelHalfSize[3], double worldRadius) const
{
  
  // Get Full-Width at Half-Maximum of weight window in millimeters.
  const double worldFWHM = _getWindowFWHM(worldRadius);
  
  // The template size must be big enough to encompass the weight window.
  // Weights above 0.2 must be included.
  const double weightThreshold = 0.2;
  // The following formula calculates the x coordinate for which
  // a Gaussian function attains the value y:
  // x = sigma*sqrt(log(1/y^2))
  // Note: there are two solutions, but only the positive one is considered.
  const double sigma = worldFWHM/2.35;
  const double xTh = sigma*sqrt(log(1/(weightThreshold*weightThreshold)));
  
  // Deduce template size in voxels
  voxelHalfSize[0] = (int)ceil(xTh / _worldVoxelSize[0]);         
  voxelHalfSize[1] = (int)ceil(xTh / _worldVoxelSize[1]);
  voxelHalfSize[2] = (_TRACKDIM==_3D)? (int)ceil(xTh/ _worldVoxelSize[2]) : 0;
}



// --------------------------------------------------------------------------------------------------
//! Get a patch of image data for fitting the vessel template vt.
// --------------------------------------------------------------------------------------------------
template<class TM> void VesselTemplateTracker3D<TM>::_getImPatch(ImPatch& I, const VesselData& vd) const
{
  // Get template parameters
  const double worldRadius = vd.getRadius();
  double x0,y0,z0;
  vd.getCenterPoint(x0,y0,z0);

  // Find the voxel coordinates for the current center point
  double voxelCenter[3] = {x0, y0, z0};
  _imProps.transformToVoxelCoord(&voxelCenter[0], &voxelCenter[1], &voxelCenter[2]);

  // Find the closest grid point (with a voxel value) in voxel coordinates
  // MeVisLab has its voxel origin at the upper left corner of the upper left voxel. Hence,
  // we need to apply a floor() operation to find integer indexes into the array holding the image.
  double voxelClosestVoxel[3] = {floor(voxelCenter[0]), floor(voxelCenter[1]), floor(voxelCenter[2])};   

  // Find the closest grid point (with a voxel value) in world coordinates
  // We have to add 0.5 here as MeVisLab has the center of the voxels at 0.5, 1.5, 2.5 etc.
  double worldClosestVoxel[3] = {voxelClosestVoxel[0]+0.5, voxelClosestVoxel[1]+0.5, voxelClosestVoxel[2]+0.5};
  _imProps.transformToWorldCoord(&worldClosestVoxel[0], &worldClosestVoxel[1], &worldClosestVoxel[2]);

  // Get a suitable patch size.
  int voxelHalfSize[3];
  _getTemplateSize(voxelHalfSize, worldRadius);
  int voxelFullSize[3] = {2*voxelHalfSize[0] + 1, 2*voxelHalfSize[1] + 1, 2*voxelHalfSize[2] + 1};
  const unsigned int nbrOfVoxels = voxelFullSize[0]*voxelFullSize[1]*voxelFullSize[2];

  // Resize the variables in I to have the correct size
  I.voxelValues.resize(nbrOfVoxels);
  I.x.resize(nbrOfVoxels);
  I.y.resize(nbrOfVoxels);
  I.z.resize(nbrOfVoxels);

  // Fill image patch I with voxel values and world coordinates.
  int counter = 0;
  double worldCoordz = worldClosestVoxel[2] - _worldVoxelSize[2]*voxelHalfSize[2];
  for(int z=-voxelHalfSize[2]; z<=voxelHalfSize[2]; ++z) {
    double worldCoordy = worldClosestVoxel[1] - _worldVoxelSize[1]*voxelHalfSize[1];
    for(int y=-voxelHalfSize[1]; y<=voxelHalfSize[1]; ++y) {
      double worldCoordx = worldClosestVoxel[0] - _worldVoxelSize[0]*voxelHalfSize[0];
      for(int x=-voxelHalfSize[0]; x<=voxelHalfSize[0]; ++x) {
        // TODO: use strides here!
        I.voxelValues[counter] = _imData.getImgVal(voxelClosestVoxel[0] + x, voxelClosestVoxel[1] + y, voxelClosestVoxel[2] + z);  
        I.x[counter] = worldCoordx;
        I.y[counter] = worldCoordy;
        I.z[counter] = worldCoordz;   
        ++counter;
        worldCoordx += _worldVoxelSize[0];
      }
      worldCoordy += _worldVoxelSize[1];
    }
    worldCoordz += _worldVoxelSize[2];
  }
}


#endif
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

#include "mlTubularTracking.h"

ML_START_NAMESPACE

//! Implements code for the runtime type system of the ML
ML_MODULE_CLASS_SOURCE(TubularTracking, Module);


TubularTracking::TubularTracking() : Module(2, 0),
_virtInputVolume(NULL),
_virtMaskVolume(NULL),
_outputChanged(false),
_trackedBoostGraphs(0) {
  handleNotificationOff();

  // Add base output for vessel graph
  _outputGraphFld = addBase("outputGraph");
  _outputGraphFld->addAllowedType<Graph>();
  _outputGraphFld->setBaseValue(&_outputGraph);

  // Add minimum radius field (millimeters)
  _minRadiusFld = addDouble("minRadius", 1);
  //_minRadiusFld->setDoubleMaxValue();
  _parameterFlds.push_back(_minRadiusFld);

  // Add maximum radius field (millimeters)
  _maxRadiusFld = addDouble("maxRadius", 3);
  //_maxRadiusFld->setDoubleMaxValue(50);
  _parameterFlds.push_back(_maxRadiusFld);

  // Add initial radius field (millimeters)
  _initRadiusFld = addDouble("initRadius", 2);
  _parameterFlds.push_back(_initRadiusFld);

  // Add flag field whether to use initial radius
  _useInitRadiusFld = addBool("useInitRadius", false);
  _parameterFlds.push_back(_useInitRadiusFld);
    
  // Add number of search angles field
  _nbrSearchAnglesFld = addInt("nbrOfSearchAngles", 3);
  _parameterFlds.push_back(_nbrSearchAnglesFld);

  //! Maximum angle relative to the previous tracking step (in degrees)
  _maxAngleFld = addDouble("maxAngle", 60);
  _parameterFlds.push_back(_maxAngleFld);

  // Threshold for a single tracking step.
  _pruningThresholdFld = addDouble("pruningThreshold", 1);
  _parameterFlds.push_back(_pruningThresholdFld);

  // Threshold for the average score of a tracking path of
  // length _searchDepth. Should be higher than _scoreThreshold.
  _terminationThresholdFld = addDouble("terminationThreshold", 3);
  _parameterFlds.push_back(_terminationThresholdFld);

  // Add maximum number of tracking steps field
  _maxNbrStepsFld = addInt("maxNbrSteps", 100);
  _parameterFlds.push_back(_maxNbrStepsFld);

  // Maximum length in millimeter to track.
  _maxLengthFld = addDouble("maxLength", 100);
  _parameterFlds.push_back(_maxLengthFld);

  // Set search depth in tracking tree
  _searchDepthFld = addInt("searchDepth", 3);
  _parameterFlds.push_back(_searchDepthFld);

  // Minimum distance between two branchings in millimeters.
  _minBranchingDistanceFld = addDouble("minBranchingDistance", 5);
  _parameterFlds.push_back(_minBranchingDistanceFld);

  // Tracking step length factor
  _stepLengthFactorFld = addDouble("stepLengthFactor", 1);
  _parameterFlds.push_back(_stepLengthFactorFld);

  // Size of the Gaussian window used for fitting
  _windowSizeFactorFld = addDouble("windowSizeFactor", 5);
  _parameterFlds.push_back(_windowSizeFactorFld);


  // Add toggle field that enables/disables vessel branching.
  _allowBranchingFld = addBool("allowBranching", false);
  _parameterFlds.push_back(_allowBranchingFld);

  // Add toggle field for growing bidirectionally.
  _growBidirectionalFld = addBool("growBidirectional", false);
  _parameterFlds.push_back(_growBidirectionalFld);

  //! Single hypothesis or multiple hypothesis tracking.
  _useMultipleHypothesesTrackingFld = addBool("useMultipleHypotheses", false);
  _parameterFlds.push_back(_useMultipleHypothesesTrackingFld);

  // Add toggle field that enables/disables termination after a certain number of tracking steps.
  _toggleMaxStepsFld = addBool("toggleMaxSteps", false);
  _parameterFlds.push_back(_toggleMaxStepsFld);

  // Add toggle field that enables/disables termination after a certain tracking length.
  _toggleMaxLengthFld = addBool("toggleMaxLength", false);
  _parameterFlds.push_back(_toggleMaxLengthFld);

  // Add base input for initialization points
  _inputInitPointsFld = addBase("inputInitPoints", NULL);
  _inputInitPointsFld->addAllowedType<XMarkerList>();

  // Add base input for tracked points
  _outputTrackedPointsFld = addBase("outputTrackedPoints");
  _outputTrackedPointsFld->addAllowedType<XMarkerList>();
  _outputTrackedPointsFld->setBaseValue(&_trackedPointsXMarkerList);

  // Add button for starting/updating calculations.
  _updateButtonFld = addNotify("update");

  //! Add for clearing previous tracking result.
  _clearButtonFld = addNotify("clear");

  //! Clear binary image.
  _binaryImg.clear(false);

  handleNotificationOn();
}

TubularTracking::~TubularTracking() {
  // Delete allocated images and virtual volumes
  if (_virtInputVolume) {
    delete _virtInputVolume;
    _virtInputVolume = NULL;
  }
  if (_virtMaskVolume) {
    delete _virtMaskVolume;
    _virtMaskVolume = NULL;
  }
  _trackedPointsXMarkerList.clearList();
  _outputGraph.clearGraph();
  _trackedBoostGraphs.clear();
}

//----------------------------------------------------------------------------------
// Small helper functions to clamp field values. 
//----------------------------------------------------------------------------------
inline void TT_DOUBLE_MIN_CLAMP(DoubleField* df, double minVal)
{
  if (df && (df->getDoubleValue() < minVal)){
    df->setDoubleValue(minVal);
  }
}
inline void TT_DOUBLE_CLAMP(DoubleField* df, double minVal, double maxVal)
{
  TT_DOUBLE_MIN_CLAMP(df, minVal);
  if (df && (df->getDoubleValue() > maxVal)){
    df->setDoubleValue(maxVal);
  }
}
inline void TT_MLINT_CLAMP(IntField* intf, MLint minVal, MLint maxVal)
{
  if (intf){
    if (intf->getIntValue() < minVal){
      intf->setIntValue(minVal);
    }
    if (intf->getIntValue() > maxVal){
      intf->setIntValue(maxVal);
    }
  }
}

void TubularTracking::handleNotification (Field *field) {
  // Clamp some fields against limits or intervals.
  TT_DOUBLE_MIN_CLAMP(_minRadiusFld, 0);
  TT_DOUBLE_MIN_CLAMP(_maxRadiusFld, 0);
  TT_MLINT_CLAMP(_nbrSearchAnglesFld, 0, 20);
  TT_DOUBLE_CLAMP(_maxAngleFld, 0, 85);
  TT_DOUBLE_CLAMP(_pruningThresholdFld, -1000, 1000);
  TT_DOUBLE_CLAMP(_terminationThresholdFld, -1000, 1000);
  TT_MLINT_CLAMP(_maxNbrStepsFld, 0, 100000);
  TT_DOUBLE_CLAMP(_maxLengthFld, 0, 1000000);
  TT_MLINT_CLAMP(_searchDepthFld, 1, 10);
  TT_DOUBLE_CLAMP(_minBranchingDistanceFld, 0, 1000);
  TT_DOUBLE_CLAMP(_stepLengthFactorFld, 0.01, 5);
  TT_DOUBLE_MIN_CLAMP(_windowSizeFactorFld, 3);
  
  // Check if one of the algorithm parameters have been changed or
  // if a new image is attached. The stored tracking result is then
  // invalid and we need to retrack everything on the next update.
  // First, check if a parameter has been changeed.
  bool parameterChanged = false;
  std::vector<Field*>::iterator it = _parameterFlds.begin();
  for(; it !=_parameterFlds.end(); ++it) {
    if(*it == field) {
      parameterChanged = true;
      break;
    }
  }

  // Has input images changed?
  bool inputChanged = (field == getInputImageField(0)) || (field == getInputImageField(1));

  // If the clear button has been presses, just clear outputs but keep
  // the already tracked vessels in _trackedBoostGraphs.
  if( (field ==_clearButtonFld) || parameterChanged || inputChanged) {

    // Clear output data structures. These will be filled the next time
    // the "Update" button of this module is pressed.
    _trackedPointsXMarkerList.clearList();
    _outputGraph.clearGraph();
    _binaryImg.reset(); 

    // If an input parameter of image has changed, the internal
    // tracking graph is also invalid and should be cleard. The
    // alternative is that the "Clear" button has been pressed,
    // then the internal structure is still valid.
    if(parameterChanged || inputChanged) {
      _trackedBoostGraphs.clear();
    }

    // If output has changed we send notifications, otherwise not.
    if(_outputChanged) {
      _outputTrackedPointsFld->touch(); // Output XMarker list with tracked points
      _outputGraphFld->touch();         // Ouput vessel graph
      _outputChanged = false;
    }
  }

  // If a new data image is attached, create new virtual volume
  if( field == getInputImageField(0) ) {

    // Delete previously allocated virtual volume
    if (_virtInputVolume)  { delete _virtInputVolume; _virtInputVolume = NULL; }

    // Check that valid input image is available
    PagedImage *inImg0 = getUpdatedInputImage(0);
    if(inImg0) {
      // Create typed virtual volume for input to access voxels
      _virtInputVolume = new TVirtualVolume<float>(inImg0);   
    }
    // If there is no image available, we return without processing.
    else { return; }
  }


  // If a new mask image is attached, create new virtual volume
  if(field == getInputImageField(1)) {

    // Delete previously allocated virtual volume
    if (_virtMaskVolume)   { delete _virtMaskVolume;  _virtMaskVolume = NULL;  }

    // Check if valid mask image is available. The mask image is optional,
    // we can go on without it.
    PagedImage *inImg1 = getUpdatedInputImage(1);
    if(inImg1) {
      _virtMaskVolume  = new TVirtualVolume<MLuint8>(inImg1);
    }
  }

  // Track vessels
  if(field == _updateButtonFld) {

    // Do we have input data? Otherwise no tracking.
    PagedImage *inImg = getUpdatedInputImage(0);
    if(!inImg) {return;}

    // If _virtMaskVolume is not NULL, there is a real mask at input 1.
    // Check that this mask has same extent as image at input 0.
    if(_virtMaskVolume) {
      PagedImage *maskImg = getUpdatedInputImage(1);
      if(maskImg->getImageExtent() != inImg->getImageExtent()) {
        ML_PRINT_WARNING("TubularTracking", ML_BAD_DIMENSION, "Data and mask images must have the same size.");
        return;
      }
    }

    // Do we have 2D or 3D input data
    DIM inputDimension = (inImg->getImageExtent().z > 1)? _3D : _2D;

    // Resize the binary image containing a mask of the segmented vessels.
    _binaryImg.setExtent(inImg->getImageExtent());
    _binaryImg.clear(false);

    // Clear previous entries in the internal XMarkerLists 
    _trackedPointsXMarkerList.clearList();
    _initPointsXMarkerList.clearList(); 

    // Clear output ML graph
    _outputGraph.clearGraph();

    // Get intial tracking points as an XMarkerList
    if(!_getInputXmarkers()) {return;}


    // Check if we have already tracked from any of the input XMarkers, in such
    // case we don't have to track from this point again, the result is already
    // stored in the _trackedBoostGraphs vector. To do this, remove all entries 
    // in the current _trackedBoostGraphs vector that do not match  markers in 
    // the input XMarker list.
    std::vector<std::pair<XMarker, boostGraph> >::iterator graphIterator = _trackedBoostGraphs.begin();
    while(graphIterator!=_trackedBoostGraphs.end())
    {
      bool foundMatch = false;
      XMarkerList::iterator markerIterator = _initPointsXMarkerList.begin();
      while(markerIterator != _initPointsXMarkerList.end())
      {
        // Is the input XMarker equal to the marker in the _trackedBoostGraphs?
        // Then erase the XMarker from the input list, the result is already there.
        if((*graphIterator).first == *markerIterator) {
          foundMatch = true;
          markerIterator = _initPointsXMarkerList.erase(markerIterator);
        }
        else {
          ++markerIterator;
        }
      }

      // If we went through the list of input XMarkers and didn't find a match
      // for the current entry in the _trackedBoostGraphs, the user is not interested
      // in the stored result any longer. Remove it!
      if(!foundMatch) {
        graphIterator = _trackedBoostGraphs.erase(graphIterator);   
      }
      else {
        ++graphIterator;
      }
    }
    // At this point, the _initPointsXMarkerList only contains Xmarkers from where we 
    // previously haven't tracked and the _trackedBoostGraphs has been cleaned from
    // entries of which the user is no longer interested

    // Clear segmentation binary image 
    _binaryImg.clear(0);   

    // Create a typed subimages containing blocks of indata and maskdata
    TSubImage<float> imData;
    TSubImage<MLuint8> maskData;    // When the data pointer in maskData is NULL, no mask is supplied.

    // Get world coordinates system and other properties
    MedicalImageProperties imgProps;
    imgProps.setImageProperties(*inImg);

    // Fill binary image with already segmented data in _trackedBoostGraphs
    for(unsigned int k=0; k<_trackedBoostGraphs.size(); ++k)
    {
      // Get a more convenient name for the current boost graph
      const boostGraph& currentGraph = _trackedBoostGraphs[k].second;

      // We need a dummy instantiation of the VesselTemplateTracker because
      // the function we need is a virtual method in this class.
      // Not a beautiful solution, but works.
      // We do this already here because the binary image showing where the
      // tracking has been is required when tracking new vessels so as to
      // avoid segmenting the same areas again.
      VesselTemplateTracker dummyTracker(imData, maskData, imgProps, inputDimension);

      // Each tubular model in the vertices in currentGraph should be voxelized.
      // The typedefs etc below are very messy, blame the boost library for this!
      // All we want to do is to loop over all the vertices in the graph.
      typedef boost::property_map<boostGraph, boost::vertex_index_t>::type indexMap;
      indexMap index = boost::get(boost::vertex_index, currentGraph);
      typedef boost::graph_traits<boostGraph>::vertex_iterator vertex_iter;
      std::pair<vertex_iter, vertex_iter> vp;
      for (vp = boost::vertices(currentGraph); vp.first != vp.second; ++vp.first) {
        unsigned int vertexIdx = index[*vp.first];
        dummyTracker.setVisitedVoxels(currentGraph[vertexIdx],_binaryImg);
      }
    }

    double initRadius = _useInitRadiusFld->getBoolValue()
      ? _initRadiusFld->getDoubleValue()
      : 0.5*(_maxRadiusFld->getDoubleValue() + _minRadiusFld->getDoubleValue());
      
    // Loop over all input xmarkers. Each xmarker point results in a root in the ML graph.
    const MLssize_t numXMarkers = mlrange_cast<MLssize_t>(_initPointsXMarkerList.size());
    for(MLssize_t markerNbr=0;  markerNbr<numXMarkers; ++markerNbr) 
    {

      // Get XMarker from input list
      const XMarker* const currentMarker = static_cast<XMarker*>(_initPointsXMarkerList.getItemAt(markerNbr));

      // If no direction vector is given in the XMarker, we estimate the direction from
      // the image data using the function estimateOrientation()
      Vector3 vesselDirection;
      bool userVectorAvailable = true;
      if((currentMarker->vec).norm2() > 1e-5) {
        vesselDirection = currentMarker->vec;
      }
      else { 
        // Estimate image structure orientation
        vesselDirection = estimateOrientation((currentMarker->pos).getVec3(), imgProps, *_virtInputVolume, 
          initRadius, inputDimension);

        // No user direction available
        userVectorAvailable = false;
      }

      // Normalize vessel direction
      vesselDirection.normalize();

      // Do we need to grow bidirectionally?
      const bool growBidirectional = !userVectorAvailable || _growBidirectionalFld->isOn();
        
      // Create an initial VesselData structure from the first input XMarker
      VesselData firstModel;
      firstModel.setCenterPoint(currentMarker->x(),currentMarker->y(), currentMarker->z());
      firstModel.setDirection(vesselDirection[0], vesselDirection[1], vesselDirection[2]); 
      firstModel.setRadius(initRadius);

      // Declare a boost graph containing the tubular models tracked from this start point
      boostGraph treeGraph;

      // Decalre a vector containing the root models, could be 1 or 2 depending on whether
      // we are tracking uni-directional or bi-directional
      std::vector<Vertex> rootVertices;

      // Add the start-model to the boost tree graph
      Vertex firstVertex = boost::add_vertex(treeGraph);
      rootVertices.push_back(firstVertex);
      treeGraph[firstVertex] = firstModel;

      // If no user vector for the tracking direction was added, we track in two directions. We also
      // do this if the user has indicated the Bidirectional tracking option.
      if (growBidirectional) {

        // Create a vessel model segment that is the same as the initial model but with reversed direction.
        VesselData secondModel;
        secondModel.setCenterPoint(currentMarker->x(),currentMarker->y(), currentMarker->z());
        secondModel.setDirection(-vesselDirection[0], -vesselDirection[1], -vesselDirection[2]); 
        secondModel.setRadius(initRadius);

        // Add a new vertex to the boost graph containing the secondModel.
        // Connect this new vertex with the first vertex containing the initial model. The second vertex
        // will be seen as a non-terminated leaf that should be tracked from below.
        Vertex secondVertex = boost::add_vertex(treeGraph);
        treeGraph[secondVertex] = secondModel;
        //rootVertices.push_back(secondVertex);                   // TODO: Remove this line and the for loop below is unnecessary
        boost::add_edge(firstVertex,secondVertex,treeGraph);
      }

      // Set the maximum number of steps and the maximum length to track in millimeters.
      // If the user has defined values we use them, otherwise we set large values.
      const std::pair<int,double> maxTrackingLength(_toggleMaxStepsFld->isOn()  ? _maxNbrStepsFld->getIntValue() : 30000,
        _toggleMaxLengthFld->isOn() ? _maxLengthFld->getDoubleValue() : 1000000.0);

      // Max image block size we are allowed to use.
      const int maxBlockSize = 200;

      // Loop over 1 or 2 root models depending on if we are tracking bi-directionally or not.
      for(unsigned int rootVertexNbr = 0; rootVertexNbr < rootVertices.size(); ++rootVertexNbr) {   // TODO: This loop can be removed!

        // Image block size in voxels. Image blocks with side of this size will be supplied
        // to the tracking algorithm.
        int blockSize = 100;

        // Number of tracked steps and tracked path length in millimeters
        std::pair<int, double> trackingLength(0,0);

        // Initial vertex is set to the root vertex
        const Vertex rootVertex = rootVertices[rootVertexNbr];
        Vertex initVertex = rootVertex;

        bool trackingDone = false;
        while(!trackingDone)
        {
          // Fill image area with data from input image
          double x,y,z;
          treeGraph[initVertex].getCenterPoint(x, y, z);
          double currentPoint[3] = {x,y,z};
          _getImageData(imData, maskData, imgProps, currentPoint, blockSize,inputDimension);

          // Set up tracker
          VesselTemplateTracker tracker(imData, maskData, imgProps, inputDimension);
          tracker.setMinRadius(_minRadiusFld->getDoubleValue());
          tracker.setMaxRadius(_maxRadiusFld->getDoubleValue());
          tracker.setMaxSearchAngle(_maxAngleFld->getDoubleValue()/360*2*3.14159265358);  
          tracker.setTerminationThreshold(_terminationThresholdFld->getDoubleValue());
          tracker.setMinBranchingDistance(_minBranchingDistanceFld->getDoubleValue()); 
          tracker.setStepLength(_stepLengthFactorFld->getDoubleValue());
          tracker.setWindowSizeFactor(_windowSizeFactorFld->getDoubleValue());
          // Tracker setup for multiple hypothesis tracking
          if(_useMultipleHypothesesTrackingFld->isOn()) {
            tracker.setSearchDepth(_searchDepthFld->getIntValue());
            tracker.setNbrSearchAngles(_nbrSearchAnglesFld->getIntValue());
            tracker.setPruningThreshold(_pruningThresholdFld->getDoubleValue());
            if(_allowBranchingFld->isOn()) {tracker.allowBranching();} else {tracker.disallowBranching();}
          }
          // Tracker setup for single hypothesis tracking
          else {
            tracker.setSearchDepth(1);              // Search depth 1 for single hypothesis
            tracker.setNbrSearchAngles(0);          // 0 search angles => 2*0+1 = 1 prediction along the tangential line
            tracker.disallowBranching();            // Branching cannot be detected in single hypothesis mode
            tracker.setPruningThreshold(-100000.0); // Set pruning threshold low so that the prediction is always fitted
          }


          // Track
          handleNotificationOff();
          const int prevNbrSteps = trackingLength.first;
          unsigned int retCode = tracker.track(treeGraph, initVertex, trackingLength, maxTrackingLength, _binaryImg);
          const int nbrStepsTaken = trackingLength.first - prevNbrSteps;
          handleNotificationOn();

          //std::cout << "Terminatation reason: " << tracker.getTerminationMessage(retCode) << std::endl;

          // Have taken the maximum number of steps or tracked the maximum length, aborting
          if(  (MAX_NBR_STEPS_TAKEN == retCode) || (MAX_LENGTH_REACHED == retCode)) {
            break;
          }

          // The image data block supplied to the tracking algorithm is not big enough to carry out the tracking,
          // try increasing the block size. If not allowed, abort.        
          if(INSUFFICIENT_IMAGE_DATA == retCode) {
            blockSize += 50;
            if(blockSize > maxBlockSize) {
              ML_PRINT_WARNING("TubularTracking", ML_BAD_PARAMETER, "Block size for tracking too big.");
              break;
            }
          }

          // If tracking is terminating because we reach the border in just a few steps (<10)
          // we should increase the image block size.
          if( (TOO_CLOSE_TO_BORDER == retCode) && (nbrStepsTaken < 10) && (blockSize + 50 <= maxBlockSize)) {
            blockSize += 50;
          }

          // Search the resulting _treeGraph for leafs that are no terminators, that is, branches
          // where the tracking was terminated because we hit the allocated image border. We need
          // to extract these leafs retrack from these points with an updated image block. 
          // The ActiveLeafRecorder class is also used within the VesselTracker class, and this code
          // was copied from there. The VertexScorePairs is a vertex index and score pair defined
          // in VesselGraph.h. Here we do not need the score or the leafDepth, just the vertex indexes.
          VertexScorePairs leafs;
          int leafDepth = 0;
          boost::breadth_first_search(treeGraph, rootVertex, boost::visitor(ActiveLeafRecorder(&leafs, &leafDepth)));
          if(leafs.size()>0) {
            initVertex = leafs[0].first;
          }
          else {
            trackingDone = true;
          }
        } // end-while !trackingDone from the current input XMarker
      } // end-for over roots (for bi-directional tracking)

      // Add XMarker and tracked boostGraph to the _trackedBoostGraphs vector
      _trackedBoostGraphs.push_back(std::make_pair(*currentMarker,treeGraph));

    }  // end-for over input XMarkers


    // ==================
    // DONE WITH TRACKING
    // ==================

    // Convert all boost graphs in _trackedBoostGraphs to output XMarkers and a Vessel Graph
    for(unsigned int k=0; k<_trackedBoostGraphs.size(); ++k)
    {

      // Get a more convenient name for the current boost graph
      const boostGraph& currentGraph = _trackedBoostGraphs[k].second;

      // Convert tracking graph, which is a boost graph, to an ML Graph.
      boost::depth_first_search(currentGraph, boost::visitor(boostToML(_outputGraph,currentGraph)));
      // Clean graph. Joins edges: x-------x-------x   to  x---------------x
      joinEdgesAtRoots(_outputGraph);

      // Copy points from boost graph to XMarkerList
      typedef boost::property_map<boostGraph, boost::vertex_index_t>::type indexMap;
      indexMap index = boost::get(boost::vertex_index, currentGraph);

      typedef boost::graph_traits<boostGraph>::vertex_iterator vertex_iter;
      std::pair<vertex_iter, vertex_iter> vp;
      for (vp = boost::vertices(currentGraph); vp.first != vp.second; ++vp.first) {
        size_t vertexIdx = index[*vp.first];
        double px,py,pz,vx,vy,vz;
        currentGraph[vertexIdx].getCenterPoint(px,py,pz);
        currentGraph[vertexIdx].getDirection(vx, vy, vz);
        const double radius = currentGraph[vertexIdx].getRadius();
        vx *= radius; vy *= radius; vz *= radius; // Scale direction with the radius
        _trackedPointsXMarkerList.appendItem(XMarker(Vector6(px,py,pz,0,0,0),Vector3(vx,vy,vz)));
        //std::cout<< " Radius: " << currentGraph[vertexIdx].getRadius() << ". Score: " << currentGraph[vertexIdx].getScore() << std::endl;
      }
    }

    // Free allocated image data
    ML_TRY {
      imData.free();
      maskData.free();
    }
    ML_CATCH_RETHROW;

    // Notify attached fields that new output is available
    _outputTrackedPointsFld->touch();
    _outputGraphFld->touch();

    // Set flag that output has changed.
    _outputChanged = true;

  }
}


bool TubularTracking::_getInputXmarkers()
{
  // Check if anything at all is attached to the base input
  if(!_inputInitPointsFld->isValidValue()) {
    return false;
  }

  // Check that input at base field is an XMarkerList
  if (!ML_BASE_IS_A(_inputInitPointsFld->getBaseValue(), XMarkerList))
  {     
    //ML_PRINT_WARNING("TubularTracking::_getInputXmarkers()", ML_BAD_PARAMETER, "XMarkerList expected at base input.");
    return false;
  }

  // Get pointer to the input XMarkerList
  XMarkerList* inputXMarkerList = static_cast<XMarkerList*>(_inputInitPointsFld->getBaseValue());

  if(inputXMarkerList->getSize() == 0) {
    //ML_PRINT_WARNING("TubularTracking::_getInputXmarkers()",  ML_BAD_PARAMETER, "Input XMarkerList with tracking intialization points is empty.");
    return false;
  }

  // Check that the XMarkerList contains correct information and copy the init
  // points to an internal XMarkerList
  const MLssize_t numXMarkers = mlrange_cast<MLssize_t>(inputXMarkerList->getSize());
  for (MLssize_t i=0; i < numXMarkers; ++i) {
    XMarker *currentMarker = static_cast<XMarker*>(inputXMarkerList->getItemAt(i)); 
    _initPointsXMarkerList.appendItem(*currentMarker);
  }
  return true;
}


Vector3 TubularTracking::estimateOrientation(const Vector3& pos, const MedicalImageProperties& imgProps, TVirtualVolume<float>& virtDataVol, 
                                             double radius, DIM trackDim)
{

  // Calculate Hessian for 3 different scales [hme: where? I only see one scale...]

  // Get voxel size
  const Vector3 voxelSize = imgProps.getVoxelSize();

  // Get closest voxel position
  double closestVoxel_x = pos[0];
  double closestVoxel_y = pos[1];
  double closestVoxel_z = pos[2];
  const Vector3 vc = imgProps.mapWorldToVoxel(Vector3(closestVoxel_x, closestVoxel_y, closestVoxel_z));
  closestVoxel_x = vc[0];
  closestVoxel_y = vc[1]; 
  closestVoxel_z = vc[2];
  
  closestVoxel_x = floor(closestVoxel_x);
  closestVoxel_y = floor(closestVoxel_y);
  closestVoxel_z = floor(closestVoxel_z);

  // Find size of image so that we don't read outside the border
  ImageVector imExt = virtDataVol.getVirtualVolume().getExtent();

  // For each scale
  std::vector<double> kernelx;
  std::vector<double> kernely;
  std::vector<double> kernelz;

  // Get sigma in voxels
  const double sigma_x = radius/voxelSize[0]; 
  const double sigma_y = radius/voxelSize[1]; 
  const double sigma_z = (trackDim == _2D)? 0.0 : radius/voxelSize[2]; 

  // Array for keeping 2D (indexes 0,1,2) or 3D (all indexes) Hessian matrix.
  double hessian[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  if(trackDim == _2D) {
    // Calculate Hessian components
    for(unsigned int hessianComponent=0; hessianComponent<3; ++hessianComponent) {

      get1DGaussKernel(sigma_x, kernelx, hessian2DOrder[hessianComponent][0]);
      get1DGaussKernel(sigma_y, kernely, hessian2DOrder[hessianComponent][1]);
      const int halfSize_x = ((int)kernelx.size() - 1)/2;
      const int halfSize_y = ((int)kernely.size() - 1)/2;

      for(int y=-halfSize_y; y<=halfSize_y; ++y) {
        int x = -halfSize_x;
        virtDataVol.setCursorPosition(ImageVector(closestVoxel_x + x, closestVoxel_y + y, closestVoxel_z, 0, 0, 0));
        for( ; x<=halfSize_x; ++x) {
          // Check first that position is within image, otherwise we have a crash in .getCursorVoxel()
          const float imVal = ((closestVoxel_x + x >= 0) && (closestVoxel_x + x < imExt.x) &&
                               (closestVoxel_y + y >= 0) && (closestVoxel_y + y < imExt.y)) ?
                              virtDataVol.getCursorValue() : 0.0; 
          hessian[hessianComponent] += kernelx[halfSize_x + x]*kernely[halfSize_y + y]*imVal;
          virtDataVol.moveCursorX();
        } // end-for over x
      } // end-for over y 
    } // end-for over 3 Hessian components 

    

    // Eigen-decomposition of Hessian matrix, want eigenvector associated with
    // the eigenvalue with the smallest magnitude. This vector points out the
    // direction in which the image varies the least. On a tube-shaped object,
    // this should be the tube direction.
    SymmetricMatrix H(2);
    H(1,1) = hessian[0]; H(1,2) = hessian[1]; 
    H(2,2) = hessian[2];
    DiagonalMatrix D(2);
    Matrix E(2,2);
    eigenvalues(H,D,E);

    // The eigenvalues are sorted so that D(1)<D(2). Find the eigenvalue with
    // the smallest magnitude and return the associated eigenvector.
    if( std::abs(D(1)) < std::abs(D(2))) {
      return Vector3(E(1,1), E(2,1), 0.0);
    }
    else { // Magnitude of D(2) is smaller than magnitude of D(1), return eigenvector 2.
      return Vector3(E(1,2), E(2,2), 0.0);
    }

    // mlDebug("Hessian matrix: " << hessian[0] << " " << hessian[1] << " " << hessian[2]);
    // mlDebug("l1: " << D(1,1) << ". l2: " << D(2,2));
    // mlDebug("e1: ("<<E(1,1)<<", " << E(2,1) << ")");
    // mlDebug("e2: ("<<E(1,2)<<", " << E(2,2) << ")");



  } //end-if trackDim == _2D
  else if (trackDim == _3D) {

    // Calculate Hessian components
    for(unsigned int hessianComponent=0; hessianComponent<6; ++hessianComponent) {

      get1DGaussKernel(sigma_x, kernelx, hessian3DOrder[hessianComponent][0]);
      get1DGaussKernel(sigma_y, kernely, hessian3DOrder[hessianComponent][1]);
      get1DGaussKernel(sigma_z, kernelz, hessian3DOrder[hessianComponent][2]);
      const int halfSize_x = ((int)kernelx.size() - 1)/2;
      const int halfSize_y = ((int)kernely.size() - 1)/2;
      const int halfSize_z = ((int)kernelz.size() - 1)/2;

      for(int z=-halfSize_z; z<=halfSize_z; ++z) {
        for(int y=-halfSize_y; y<=halfSize_y; ++y) {
          int x = -halfSize_x;
          virtDataVol.setCursorPosition(ImageVector(closestVoxel_x + x, closestVoxel_y + y, closestVoxel_z + z, 0, 0, 0));
          const double yzKernelValue = kernely[y + halfSize_y]*kernelz[z + halfSize_z];
          for( ; x<=halfSize_x; ++x) {
            // Check first that position is within image, otherwise we have a crash in .getCursorVoxel()
            const float imVal = ((closestVoxel_x + x >= 0) && (closestVoxel_x + x < imExt.x) &&
                                 (closestVoxel_y + y >= 0) && (closestVoxel_y + y < imExt.y) &&
                                 (closestVoxel_z + z >= 0) && (closestVoxel_z + z < imExt.z)) ? 
                                  virtDataVol.getCursorValue() : 0.0; 
            hessian[hessianComponent] += kernelx[x + halfSize_x]*yzKernelValue*imVal;
            virtDataVol.moveCursorX();
          } // end-for over x
        } // end-for over y 
      } // end-for over z
    } // end-for over 6 Hessian components 

    // Eigen-decomposition of Hessian matrix, want eigenvector associated with
    // the eigenvalue with the smallest magnitude. This vector points out the
    // direction in which the image varies the least. On a tube-shaped object,
    // this should be the tube direction.
    SymmetricMatrix H(3);
    H(1,1) = hessian[0]; H(1,2) = hessian[1]; H(1,3) = hessian[2];
    H(2,2) = hessian[3]; H(2,3) = hessian[4];
    H(3,3) = hessian[5];
    DiagonalMatrix D(3);
    Matrix E(3,3);
    eigenvalues(H,D,E);

    // The eigenvalues are sorted so that D(1)<D(2)<D(3). Find the eigenvalue with
    // the smallest magnitude and return the associated eigenvector.
    if( (std::abs(D(1)) < std::abs(D(2))) && (std::abs(D(1)) < std::abs(D(3))) ) {
      return Vector3(E(1,1), E(2,1), E(3,1));
    }
    else if( (std::abs(D(2)) < std::abs(D(1))) && (std::abs(D(2)) < std::abs(D(3))) ) {
      return Vector3(E(1,2), E(2,2), E(3,2));
    }
    else {
      return Vector3(E(1,3), E(2,3), E(3,3));
    }
  } //end-if trackDim == _3D

  ML_PRINT_WARNING("TubularTracking::estimateOrientation", 
    ML_PROGRAMMING_ERROR,
    "Should have returned from function before reaching this point - let the module author know!");

  return Vector3(0);
}


void TubularTracking::_getImageData(TSubImage<float> &imData, TSubImage<MLuint8> &maskData, const MedicalImageProperties &imgProps, double worldCenterPoint[3], int blockSize, DIM trackDim) const
{
  // Free previously allocated image data
  ML_TRY{
    imData.free();
    maskData.free();
  }
  ML_CATCH_RETHROW;

  // Convert the centerPoint to the voxel coordinate system
  double voxelCenterPoint[3] = {worldCenterPoint[0], worldCenterPoint[1], worldCenterPoint[2]};
  const Vector3 vc = imgProps.mapWorldToVoxel(Vector3(voxelCenterPoint[0], voxelCenterPoint[1], voxelCenterPoint[2]));
  voxelCenterPoint[0] = vc[0];
  voxelCenterPoint[1] = vc[1]; 
  voxelCenterPoint[2] = vc[2];

  // Round up and down respectively to get voxel index
  int voxelBox_v1[3];
  int voxelBox_v2[3];
  voxelBox_v1[0] = (int)floor(voxelCenterPoint[0] - blockSize/2);
  voxelBox_v1[1] = (int)floor(voxelCenterPoint[1] - blockSize/2);
  voxelBox_v1[2] = (trackDim==_2D) ? (int)voxelCenterPoint[2] : (int)floor(voxelCenterPoint[2] - blockSize/2);
  voxelBox_v2[0] = (int)floor(voxelCenterPoint[0] + blockSize/2 + 1);
  voxelBox_v2[1] = (int)floor(voxelCenterPoint[1] + blockSize/2 + 1);
  voxelBox_v2[2] = (trackDim==_2D) ? (int)voxelCenterPoint[2] : (int)floor(voxelCenterPoint[2] + blockSize/2 + 1);

  // Get suitable box for tracking
  // MeVisLab has its voxel origin at the upper left corner of the upper left voxel. Hence,
  // we need to apply a floor() operation to find integer indexes into the array holding the image.
  ImageVector v1(voxelBox_v1[0], voxelBox_v1[1], voxelBox_v1[2], 0, 0, 0);
  ImageVector v2(voxelBox_v2[0], voxelBox_v2[1], voxelBox_v2[2], 0, 0, 0);
  SubImageBox imBox(v1,v2);

  // Copy image data
  imData.setBox(imBox);             
  ML_TRY{
    imData.allocate(ML_FATAL_MEMORY_ERROR);
    imData.fill(0.0f);                     // Fill with zeros (This is necessary!!)
    _virtInputVolume->copyToSubImage(imData);     // Copy input image data to subimage
  }
  ML_CATCH_RETHROW;

  // Copy mask data if a mask image is attached.
  if(_virtMaskVolume) {
    maskData.setBox(imBox);             
    ML_TRY{
      maskData.allocate(ML_FATAL_MEMORY_ERROR);
      maskData.fill(0);                     // Fill with zeros (This is necessary!!)
      _virtMaskVolume->copyToSubImage(maskData);     // Copy input image data to subimage
    }
    ML_CATCH_RETHROW;
  }
  // If no mask is supplied by the user at input 1, we set the data pointer in maskData to NULL.
  else{
    maskData.setData(NULL);
  }
}

Module::INPUT_HANDLE TubularTracking::handleInput(int inIndex, Module::INPUT_STATE /*state*/) const {
  // If input 1 is not connected, connect to dummy input.
  if(inIndex == 1) {  
    return Module::ALLOW_INVALID_INPUT;
  }
  else {
    return Module::INVALIDATE;
  }
}


ML_END_NAMESPACE

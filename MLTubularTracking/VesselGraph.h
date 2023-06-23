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
#ifndef __VesselTrackingGraph_H
#define __VesselTrackingGraph_H

#include "VesselData.h"

// MeVisLab includes
#include "mlGraph.h"
#include "mlVesselNode.h"
#include "mlVesselEdge.h"
#include "mlSkeleton.h"


// Includes from boost library
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>


// Create a graph type where each vertex contains a templateData and each
// edge a double (for example the length). // TODO: How to define empty data structure on edges // FIXME "double"?
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VesselData, int> boostGraph;
typedef boost::graph_traits<boostGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<boostGraph>::edge_descriptor Edge;
typedef boost::graph_traits<boostGraph>::out_edge_iterator itOutEdge;
typedef std::vector<std::pair<Vertex,double> >  VertexScorePairs;

// A boost graph visitor that simply gathers all vertexes in a list.
// Is typically used to find all vertexes in a subgraph.
class GetVertexes : public boost::default_bfs_visitor
{
public:
  GetVertexes(std::vector<Vertex>* vertexes) : _vertexes(vertexes) {_vertexes->clear(); }

  void discover_vertex(Vertex v, const boostGraph& /*g*/)
  {
    _vertexes->push_back(v);
  }

private:
  std::vector<Vertex>* _vertexes;
};


// Declaration
void copyBoostGraph(boostGraph& targetGraph, const Vertex targetRoot, const boostGraph& sourceGraph, const Vertex sourceRoot);

// The following Matlab code is used to generate the
// "hot" colorscale:
// >> h = hot(50);h=h(11:end-10,:);imagesc([1:30]),colormap(h)
const int NBROFCOLORS = 30;
const double hot[NBROFCOLORS][3] = {{0.6111, 0, 0}, {0.6667, 0, 0},
{0.7222,         0,         0},
{0.7778,         0,         0},
{0.8333,         0,         0},
{0.8889,         0,         0},
{0.9444,         0,         0},
{1.0000,         0,         0},
{1.0000,    0.0556,         0},
{1.0000,    0.1111,         0},
{1.0000,    0.1667,         0},
{1.0000,    0.2222,         0},
{1.0000,    0.2778,         0},
{1.0000,    0.3333,         0},
{1.0000,    0.3889,         0},
{1.0000,    0.4444,         0},
{1.0000,    0.5000,         0},
{1.0000,    0.5556,         0},
{1.0000,    0.6111,         0},
{1.0000,    0.6667,         0},
{1.0000,    0.7222,         0},
{1.0000,    0.7778,         0},
{1.0000,    0.8333,         0},
{1.0000,    0.8889,         0},
{1.0000,    0.9444,         0},
{1.0000,    1.0000,         0},
{1.0000,    1.0000,    0.0714},
{1.0000,    1.0000,    0.1429},
{1.0000,    1.0000,    0.2143},
{1.0000,    1.0000,    0.2857}};

// VesselEdge is an ML Graph edge
// VesselNode is an ML Graph junction vertex
// Vertex is a boost graph vertex
// Edge is a boost graph edge
// NOTE: It is important that a depth-first search is done.
class boostToML : public boost::default_dfs_visitor
{
public:
  boostToML(ML_NAMESPACE::Graph& mlGraph, const boostGraph& g)
  : _scoreAccessor(mlGraph.getSkeletonPropertyManager(), "Score", -1.0),
    _nbrNodes(0),
    _mlGraph(mlGraph),
    _currentEdge(NULL)
  {
    // Get max and min scores of the tubular segments in order to color the vessel graph.
    typedef boost::property_map<boostGraph, boost::vertex_index_t>::type indexMap;
    indexMap index = boost::get(boost::vertex_index, g);
    typedef boost::graph_traits<boostGraph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;

    // Loop over the vertexes in the boost graph and extract min and max scores
    _minScore = 1000000;
    _maxScore = -1000000;
    vp = boost::vertices(g);
    _initScore = g[index[*vp.first]].getScore();
    for (; vp.first != vp.second; ++vp.first)
    {
      unsigned int vertexIdx = index[*vp.first];
      double score = g[vertexIdx].getScore();
      if(score > _maxScore) { _maxScore=score;}
      if(score < _minScore) { _minScore=score;}
    }
  }

    // copy constructor needed because _scoreAccessor is not copyable
  boostToML(const boostToML &other)
  : boost::default_dfs_visitor(other),
    _minScore(other._minScore),
    _maxScore(other._maxScore),
    _initScore(other._initScore),
    _scoreAccessor(other._mlGraph.getSkeletonPropertyManager(), "Score", -1.0),
    _nbrNodes(other._nbrNodes),
    _mlGraph(other._mlGraph),
    _currentEdge(other._currentEdge),
    _edgeQueue(other._edgeQueue)
  {
  }

  void discover_vertex(Vertex v, const boostGraph& g)
  {
    // Get coordinates from boost graph
    double x,y,z;
    g[v].getCenterPoint(x,y,z);

    // Get vessel radius from boost graph.
    double radius = g[v].getRadius();

    // Create a skeleton to be added to the ML Graph edge
    ML_NAMESPACE::Vector3 vec3Point(x,y,z);

    // Set color of vessel.
    double score = g[v].getScore();
    int colorIndex = (int)floor((score - _minScore)/(_initScore) * (double)NBROFCOLORS);
    if(colorIndex<0) { colorIndex=0;}
    if(colorIndex>=NBROFCOLORS) {colorIndex=NBROFCOLORS-1;}

    // Add a skeleton to the current edge
    if(_currentEdge != NULL) {
      ML_NAMESPACE::Skeleton* skeleton = _currentEdge->createSkeleton(vec3Point, radius, radius);
      skeleton->setRGBA(hot[colorIndex][0], hot[colorIndex][1], hot[colorIndex][2], 1);
      skeleton->setLabel(1);
      _scoreAccessor.set(*skeleton, score);
    }


    // If vertex is the root node, a junction or a terminator vertex (i.e, a graph leaf), a new node
    // must be created in the ML VesselGraph.
    const int nbrOutEdges = (int)boost::out_degree(v, g);
    if( (_nbrNodes == 0) || (nbrOutEdges != 1) )
    {
      // Add node to ML Graph
      ML_NAMESPACE::VesselNode *newNode = _mlGraph.createNode(vec3Point);

      // Close the current edge. For the boost root vertex there is no
      // current edge and the pointer _currentEdge should be NULL
      if(_currentEdge != NULL) {
        _currentEdge->setSucc(newNode);
      }
      // For the boost vertex we should instead make a ML Graph root.
      else {
        _mlGraph.insertRoot(newNode);
      }

      // Create as many ML Graph edges as there are output edges from the boost graph root
      for(int i=0; i<nbrOutEdges; ++i) {
        ML_NAMESPACE::VesselEdge *e = _mlGraph.createEdge(newNode, newNode);
        if(_currentEdge != NULL) {
          ML_NAMESPACE::Skeleton* skeleton = _currentEdge->createSkeleton(vec3Point, radius, radius);
          skeleton->setRGBA(hot[colorIndex][0], hot[colorIndex][1], hot[colorIndex][2], 1);
          skeleton->setLabel(1);
          _scoreAccessor.set(*skeleton, score);
        }
        _edgeQueue.push_back(e);
      }

      // Open up new ML Graph edge
      if(!_edgeQueue.empty())
      {
        _currentEdge = _edgeQueue.back();
        _edgeQueue.pop_back();
      }
      else
      {
        _currentEdge = NULL;
      }

      // Increment nodes
      ++_nbrNodes;
    }
  }

private:
  double _minScore, _maxScore, _initScore;
  ml::DefaultPropertyAccessor<double> _scoreAccessor;
  int _nbrNodes;
  ML_NAMESPACE::Graph& _mlGraph;
  ML_NAMESPACE::VesselEdge *_currentEdge;
  std::deque<ML_NAMESPACE::VesselEdge*> _edgeQueue;
};




//! The FindVertexPaths class is used to find the paths from a root vertex to all
//! other vertexes connected to the root, i.e., the path of vertexes one has to pass
//! from the root to the vertexes. The paths for the vertexes are stored in a C++ map<>
//! structure, where the key is the leaf vertex and the associated value is a vector<>
//! of vertexes (that is, the path). See the graph boost documentation for more
//! information on breadth-first-search (bfs) visitors.
class FindVertexPaths : public boost::default_bfs_visitor
{
public:

  //! Constructor.
  FindVertexPaths(std::map<Vertex, std::vector<Vertex> >* vertexPaths, const Vertex& root) :
  _vertexPaths(vertexPaths) // Set private member variable pointer to the input map<> structure.
  {
    //! Add the path to the root vertex to _vertexPaths
    std::vector<Vertex> rootVertex(1);
    rootVertex[0] = root;
    _vertexPaths->insert(std::make_pair(root,rootVertex));
  }

  //! When a new edge is found in the graph, update the path
  //! between source vertex and target vertex.
  void tree_edge(Edge e, const boostGraph& g) {
    // Get path to source vertex.
    std::vector<Vertex> path = (*_vertexPaths)[boost::source(e,g)];
    //! Add target vertex to the path.
    path.push_back(boost::target(e,g));
    //! Add new entry in the map<> structure where the target
    //! vertex is associated with the updated path.
    _vertexPaths->insert(std::make_pair(boost::target(e,g),path));
  }

private:
  std::map<Vertex, std::vector<Vertex> >* _vertexPaths;
};



//! This class searches the track graph and find the leafs
//! (i.e., endpoints). In addition, the score along the paths
//! to the leafs are calculated.
class ActiveLeafRecorder : public boost::default_bfs_visitor
{
public:
  ActiveLeafRecorder(std::vector<std::pair<Vertex,double> >* leafs, int* depth) :
      _leafs(leafs), _depth(depth), _depthSet(false) {}

  void discover_vertex(Vertex v, const boostGraph& g) {

    // Are the _vertexSumScore and _vertexDepth empty? In
    // such case is v the start vertex and it should be added
    // to these maps.
    if(_vertexSumScore.size() == 0) {
      _vertexSumScore.insert(std::make_pair(v,g[v].getScore()));
      _vertexDepth.insert(std::make_pair(v,0));
    }

    // Get the sumScore of this vertex
    double sumScore = _vertexSumScore[v];

    // If the current vertex has no output edges and the terminator
    // flag is not marked, it is an active leaf that should be recorded.
    if( ((int)boost::out_degree(v, g) == 0) && (!g[v].isTerminator())) {

      // If this is the first leaf we encounter, _depthSet equals false,
      // and we record the depth of the leaf. If we have encountered
      // leafs previously, _depthSet is true, and we check if the previous
      // depths are the same as for the depth of the current leaf. This
      // can be used to check that the graph is correctly built.
      if(!_depthSet) {
        *_depth = _vertexDepth[v];
        _depthSet = true;
      }
      else if (*_depth != _vertexDepth[v]) {
        *_depth = -1;                         // Warning: Unequal depths, _depth is set to -1
      }

      // Add leaf to leaf-list´.
      _leafs->push_back(std::make_pair(v,sumScore));
    }

    // Otherwise, the vertex is an intermediate vertex and we then propagate
    // the sumScore of this vertex to the connected vertexes
    else {

      // Get iterator over the out edges of the current vertex
      std::pair<itOutEdge, itOutEdge> itOutEdgePair = boost::out_edges(v,g);

      // Loop over the out edges
      for( ; itOutEdgePair.first != itOutEdgePair.second; ++itOutEdgePair.first) {

        // Get depth of this vertex
        unsigned int depth = _vertexDepth[v];

        // Get the adjacent vertex connected via this edge
        Vertex connectedVertex = boost::target(*itOutEdgePair.first,g);

        // Get score attached to this vertex
        double score = g[connectedVertex].getScore();

        // Set the sumScore of the adjacent vertex to the sum of the sumScore
        // of the current vertex plus the score attaced to the edge.
        _vertexSumScore.insert(std::make_pair(connectedVertex, sumScore + score));

        // Set the depth of the adjacent vertex to the current depth + 1
        _vertexDepth.insert(std::make_pair(connectedVertex,depth+1));
      }
    }
  }

private:
  VertexScorePairs* _leafs;
  int* _depth;
  bool _depthSet;
  std::map<Vertex, double> _vertexSumScore;
  std::map<Vertex, int> _vertexDepth;
};


ML_START_NAMESPACE
void joinEdgesAtRoots(Graph& mlGraph);
ML_END_NAMESPACE

#endif

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
#include "VesselGraph.h"
#include <map>

void copyBoostGraph(boostGraph& targetGraph, const Vertex targetInsertVertex, const boostGraph& sourceGraph, const Vertex sourceRoot)
{

  // Check if the target graph is empty. In such case the targetInsertVertex will not be used
  bool targetIsEmpty = (boost::num_vertices(targetGraph) == 0)? true : false;

  // Find all source vertexes
  std::vector<Vertex> sourceVertexes(0);
  boost::breadth_first_search(sourceGraph, sourceRoot, boost::visitor(GetVertexes(&sourceVertexes)));

  // Create vertexes in target graph
  unsigned int nbrOfVertexes = sourceVertexes.size();
  std::map<Vertex,Vertex> vertexMapping;
  for(unsigned int i=0; i<nbrOfVertexes; ++i) {
    // Add new vertex to target graph
    Vertex newVertex = boost::add_vertex(targetGraph);

    // Copy data from source vertex to target vertex
    targetGraph[newVertex] = sourceGraph[sourceVertexes[i]];

    // Associate target and source vertexes
    vertexMapping[sourceVertexes[i]] = newVertex;
  }

  // The newly inserted vertexes are not yet connected to the other vertexes in the target graph.
  // If the target graph was not empty from the beginning, add this connection
  if(!targetIsEmpty) {
    boost::add_edge(targetInsertVertex,vertexMapping[sourceRoot], targetGraph);
  }

  // Add edges between the newly created vertexes in the target graph
  for(unsigned int i=0; i<nbrOfVertexes; ++i) {

      // Get iterator over the out edges of the current vertex
      std::pair<itOutEdge, itOutEdge> itOutEdgePair = boost::out_edges(sourceVertexes[i],sourceGraph);

      // Loop over the out edges
      for( ; itOutEdgePair.first != itOutEdgePair.second; ++itOutEdgePair.first)
      {

        // Get prev and succ vertexes in the source graph
        Vertex prevVertex = boost::source(*itOutEdgePair.first,sourceGraph);
        Vertex succVertex = boost::target(*itOutEdgePair.first,sourceGraph);

        // Add a corresponding edge in the target graph
        boost::add_edge(vertexMapping[prevVertex], vertexMapping[succVertex], targetGraph);
      }
  }
}


ML_START_NAMESPACE
void joinEdgesAtRoots(Graph& mlGraph)
{

  // Iterator over roots
  const size_t nbrOfRoots = mlGraph.numRoots();

  // ImageVector of graph nodes that are going to be promoted to root nodes
  std::vector<long> newRootsId;

  // ImageVector if graph nodes that are going to be deleted
  std::vector<long> deleteEdgesId;

  // Loop over roots. Don't delete any root nodes
  // within this loop, we may mess up the root indexes
  // and the loop then. Instead, the unnecessary roots
  // are deleted in a separate loop below.
  for(size_t j=0; j<nbrOfRoots; ++j) {

    // Get root as a node
    VesselNode *currentNode = mlGraph.getRoot(j);

    // Only remove if two edges are going out
    if(currentNode->getEdgeNum() != 2) {
      continue;
    }

    /*
    // Get edges leading out from root and the connected nodes
    //             root                         root
    //              *                            *
    //             / \
    //     edge 1 /   \ edge 2      to
    //           *     *                      *-->--*
    //        node 1  node 2             node 1     node 2
     */
    VesselEdge* edge1 = static_cast<VesselEdge*>(currentNode->getDepEdge(0));
    VesselEdge* edge2 = static_cast<VesselEdge*>(currentNode->getDepEdge(1));
    VesselNode* node1 = edge1->succNode();
    VesselNode* node2 = edge2->succNode();

    // Get skeletons
    std::vector<const Skeleton*> skeletons1;
    std::vector<const Skeleton*> skeletons2;
    edge1->skeletonsDescending(skeletons1);      // Note ascending order, from node 1 to root
    edge2->skeletonsAscending(skeletons2);     // Note descending order, from root to node 2

    // Create a new edge and add the combined skeletons
    VesselEdge* combinedEdge = mlGraph.createEdge(node1, node2);
    for(size_t i = 0; i <skeletons1.size(); ++i) {
      combinedEdge->createSkeleton(skeletons1[i]);
    }
    for(size_t i = 0; i <skeletons2.size(); ++i) {
      combinedEdge->createSkeleton(skeletons2[i]);
    }

    // Note ID number of node1, it will be promoted to a root node
    newRootsId.push_back(node1->getId());

    // Note ID number of edge1 and edge2, they will be deleted later.
    deleteEdgesId.push_back(edge1->getId());
    deleteEdgesId.push_back(edge2->getId());
  }

 // Add new root nodes
 size_t i = 0;
 for(i=0; i<newRootsId.size(); ++i) {
    mlGraph.insertRoot(mlGraph.getIdNode(newRootsId[i]));
  }

  // Delete edges to obsolete root nodes. The nodes will be deleted
  // automatically when all edges from the have been deleted.
  for(i=0; i<deleteEdgesId.size(); ++i) {
     mlGraph.removeEdge(mlGraph.getIdEdge(deleteEdgesId[i]));
  }
}
ML_END_NAMESPACE

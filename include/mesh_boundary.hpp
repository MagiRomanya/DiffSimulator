#ifndef MESH_BOUNDARY_H_
#define MESH_BOUNDARY_H_

#include <vector>
#include "edge.hpp"
#include "linear_algebra.hpp"

void mesh_boundary(const std::vector<Scalar>& vertices,
                   const std::vector<unsigned int>& indices,
                   std::vector<Edge> &internalEdges,
                   std::vector<Edge> &externalEdges);

#endif // MESH_BOUNDARY_H_

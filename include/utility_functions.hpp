#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <raylib.h>
#include <vector>

#include "linear_algebra.hpp"

Camera3D create_camera(unsigned int FPS=60);

void translate_vertices(std::vector<Scalar>& vert, const Vec3& translation);

void rotate_vertices_arround_axis(std::vector<Scalar>& vert, const Vec3& rotation);

Mesh GenMeshPlaneNoGPU(float width, float length, int resX, int resZ);

#endif // UTILITY_FUNCTIONS_H_

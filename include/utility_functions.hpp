#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <iostream>
#include <raylib.h>
#include <vector>
#include <fstream>

#include "linear_algebra.hpp"
#include "physics_state.hpp"

Camera3D create_camera(unsigned int FPS=60);

void translate_vertices(std::vector<Scalar>& vert, const Vec3& translation);

void rotate_vertices_arround_axis(std::vector<Scalar>& vert, const Vec3& rotation);

Mesh GenMeshPlaneNoGPU(float width, float length, int resX, int resZ);

inline void save_current_state(const PhysicsState& state, unsigned int frame_index) {
  std::ofstream myfile;
  myfile.open("example.csv");
  if (myfile.is_open()) {
    // Fill the before states with no information
    for (unsigned int i=0; i < frame_index; i++) {
      for (unsigned int j=0; j < state.q.size()*2-1; j++) {
        myfile << 0.0 << ",";
      }
      myfile << 0.0 << "\n";
    }
    // Write the state
    for (unsigned int i = 0; i < state.q.size(); i++) {
      myfile << state.q[i] << ",";
    }
    for (unsigned int i = 0; i < state.q_dot.size(); i++) {
      myfile << state.q_dot[i];
      if (i < state.q_dot.size() - 1) {
        myfile << ",";
      }
    }
    myfile << std::endl;
    myfile.close();
    std::cout << "file_saved_successfully" << std::endl;
  }
  else {
    std::cerr << "Unable to open save state file!" << std::endl;
  }
}

#endif // UTILITY_FUNCTIONS_H_

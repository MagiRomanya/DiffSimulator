#ifndef SIMULABLE_GENERATORS_H_
#define SIMULABLE_GENERATORS_H_

#include <vector>

#include "linear_algebra.hpp"
#include "simulable.hpp"
#include "simulation.hpp"


/*
 * Construct a simulable from a mesh with mass-spring interactions
 * The vertices positions must be in world coordinates
 */
Simulable generate_mass_spring(Simulation* simulation,
                          const std::vector<Scalar>& vertices,
                          const std::vector<unsigned int>& indices,
                          Scalar node_mass,
                          Scalar stiffness,
                          Scalar bend_stiffness);

Simulable generate_mass_spring_pair(Simulation* simulation, Scalar node_mass, Scalar stiffness_value);

#endif // SIMULABLE_GENERATORS_H_

#ifndef SIMULABLE_GENERATORS_H_
#define SIMULABLE_GENERATORS_H_

#include <vector>

#include "linear_algebra.hpp"
#include "simulable.hpp"
#include "simulation.hpp"


/*
 * Construct a simulable from a mesh with mass-spring interactions
 * The vertices positions must be in world coordinates.
 *
 * The differentiable parameters are the stiffness, node stiffness and one garbage parameter.
 */
Simulable generate_mass_spring(Simulation* simulation,
                          const std::vector<Scalar>& vertices,
                          const std::vector<unsigned int>& indices,
                          Scalar node_mass,
                          Scalar stiffness,
                          Scalar bend_stiffness);

/*
 * Equally generates a simulable from a mesh using mass-springs.
 * The differentibale parameters are a stiffness value for each stiffness spring and a bend value for each bend spring.
 */
Simulable generate_mass_spring(Simulation* simulation,
                          const std::vector<Scalar>& vertices,
                          const std::vector<unsigned int>& indices,
                          Scalar node_mass,
                          const std::vector<Scalar>& stiffness_values,
                          const std::vector<Scalar>& bend_stiffness_values);

Simulable generate_mass_spring_pair(Simulation* simulation, Scalar node_mass, Scalar stiffness_value);


#endif // SIMULABLE_GENERATORS_H_

#ifndef SIMULATION_PARAMETERS_H_
#define SIMULATION_PARAMETERS_H_

#include <memory>
#include <vector>

#include "linear_algebra.hpp"

struct SimulationParameters {
    /* Values here do not change throughout the simulation */
    Scalar TimeStep = 0.1;

    // Initial conditions
    Vector q0;
    Vector q_dot0;

    // Initial state jacobian
    std::vector<Triplet> dq0_dp_triplets;
    std::vector<Triplet> dq_dot0_dp_triplets;

    // Frozen i.e. not dynamic dofs
    std::vector<unsigned int> frozen_dof;

    std::vector<Triplet> mass;

    // Differentiable parameters
    Vector p;
};

struct Parameter {
    Scalar value;
    unsigned int index;
};

#endif // SIMULATION_PARAMETERS_H_

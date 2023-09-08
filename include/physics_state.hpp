#ifndef PHYSICS_STATE_H_
#define PHYSICS_STATE_H_

#include <vector>

#include "linear_algebra.hpp"
#include "simulation_parameters.hpp"

struct PhysicsState {
    /* Values here change each time-step */
    Vector q;                   // Generalized positions
    Vector q_dot;               // Generalized velocities
    Scalar time;
};

struct EnergyDerivatives {
    EnergyDerivatives(unsigned int nDoF, unsigned int nParameters) {
        energy = 0;
        force.setZero(nDoF);
        df_dx_triplets.clear();
        df_dv_triplets.clear();
        df_dp.setZero(nDoF, nParameters);
    }
    /* Variables here can be computed from the physics state and the interactions. */
    Scalar energy;

    Vector force;

    // Force jacobians
    std::vector<Triplet> df_dx_triplets;
    std::vector<Triplet> df_dv_triplets;

    Mat df_dp;
};

/*
 * Gets the position vector in the physics state starting at index.
 * The index argument is relative to the state
*/
inline Vec3 get_particle_position(const PhysicsState* state, unsigned int index) {
    return Vec3(state->q[index], state->q[index+1], state->q[index+2]);
}

inline Vec3 get_particle_velocity(const PhysicsState* state, unsigned int index) {
    return Vec3(state->q_dot[index], state->q_dot[index+1], state->q_dot[index+2]);
}

/* Integration function
 * Calculates the next step in the simulation given the precomputed energy derivatives.
 */
void integrate(const SimulationParameters& sim_param, PhysicsState* state, const EnergyDerivatives& f);

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, SparseMatrix* mat);

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vector* vec);

#endif // PHYSICS_STATE_H_

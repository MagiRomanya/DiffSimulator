#include "linear_algebra.hpp"
#include "simulable_generators.hpp"
#include "simulation_parameters.hpp"


Simulable generate_free_particle(Simulation* simulation, Scalar mass, const Vec3& x0, const Vec3& v0) {
    SimulationParameters* sim_parameters = &simulation->simulation_parameters;
    InteractionManager* interaction_manager = &simulation->interaction_manager;
    Simulable sim {
        .index = static_cast<unsigned int>(sim_parameters->q0.size()),
        .nDoF = 3,
        .parameter_index = static_cast<unsigned int>(sim_parameters->p.size()),
        .nParameters = 6,
    };

    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->mass.push_back(Triplet(i+sim.index, i+sim.index, mass));
    }

    // Resize the DoF containers
    sim_parameters->q0.conservativeResize(sim.index + sim.nDoF);
    sim_parameters->q_dot0.conservativeResize(sim.index + sim.nDoF);

    // Resize parameter containers
    sim_parameters->p.conservativeResize(sim.parameter_index + sim.nParameters);

    // Fill the parameters witht he initial conditions
    for (size_t i=0; i < sim.nDoF; i++) {
        create_parameter(sim_parameters, x0[i], sim.parameter_index + i); // initial positions
        create_parameter(sim_parameters, v0[i], sim.parameter_index + sim.nDoF + i); // initial velocities
    }

    // Set initial positions and velocities
    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->q0[sim.index + i] = x0[i];
        sim_parameters->q_dot0[sim.index + i] = v0[i];
    }

    // Initial conditions jacobian
    for (size_t i=sim.index; i < sim.index + sim.nDoF; i++) {
        sim_parameters->dq0_dp_triplets.push_back(Triplet(i,i,1));
        sim_parameters->dq_dot0_dp_triplets.push_back(Triplet(i,i,1));
    }

    simulation->simulables.push_back(sim);
    return sim;
}

#include <cassert>
#include <memory>

#include "gravity.hpp"
#include "interaction_manager.hpp"
#include "simulable.hpp"
#include "simulable_generators.hpp"
#include "mesh_boundary.hpp"
#include "simulation.hpp"
#include "simulation_parameters.hpp"
#include "spring.hpp"

Scalar distance(const std::vector<Scalar>& vertices, unsigned int i1, unsigned int i2) {
    const Vec3 p1 = Vec3(vertices[3*i1], vertices[3*i1+1], vertices[3*i1+2]);
    const Vec3 p2 = Vec3(vertices[3*i2], vertices[3*i2+1], vertices[3*i2+2]);
    return (p1 - p2).norm();
}


Simulable generate_mass_spring_pair(Simulation* simulation, Scalar node_mass, Scalar stiffness_value) {
    SimulationParameters* sim_parameters = &simulation->simulation_parameters;
    InteractionManager* interaction_manager = &simulation->interaction_manager;

    Simulable sim {
        .index = static_cast<unsigned int>(sim_parameters->q0.size()),
        .nDoF = 6,
        .parameter_index = static_cast<unsigned int>(sim_parameters->p.size()),
        .nParameters = 2
    };

    // mass matrix
    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->mass.push_back(Triplet(i+sim.index, i+sim.index, node_mass));
    }

    // Resize the DoF containers
    sim_parameters->q0.conservativeResize(sim.index + sim.nDoF);
    sim_parameters->q_dot0.conservativeResize(sim.index + sim.nDoF);

    // Resize parameter containers
    sim_parameters->p.conservativeResize(sim.parameter_index + sim.nParameters);
    Parameter stiffness = {.value=stiffness_value, .index=sim.parameter_index};
    Parameter rest_length = {.value=1, .index=sim.parameter_index+1};
    sim_parameters->p[stiffness.index] = stiffness.value;
    sim_parameters->p[rest_length.index] = rest_length.value;

    // Set initial positions and velocities
    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->q0[sim.index + i] = 0.0;
        sim_parameters->q_dot0[sim.index + i] = 0.0;
    }
    sim_parameters->q0[sim.index] = 2.0;

    interaction_manager->m_springs.push_back(Spring(sim.index + 0, sim.index + 3, stiffness, rest_length));

    return sim;
}

Simulable generate_mass_spring(Simulation* simulation,
                          const std::vector<Scalar>& vertices,
                          const std::vector<unsigned int>& indices,
                          Scalar node_mass,
                          Scalar stiffness_value,
                          Scalar bend_stiffness_value)
{
    assert(vertices.size() % 3 == 0);

    SimulationParameters* sim_parameters = &simulation->simulation_parameters;
    InteractionManager* interaction_manager = &simulation->interaction_manager;
    const Simulable sim {
        .index = static_cast<unsigned int>(sim_parameters->q0.size()),
        .nDoF = static_cast<unsigned int>(vertices.size()),
        .parameter_index = static_cast<unsigned int>(sim_parameters->p.size()),
        .nParameters = 3
    };

    // Set up mass matrix
    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->mass.push_back(Triplet(i+sim.index, i+sim.index, node_mass));
    }

    // Resize the DOF containers
    sim_parameters->q0.conservativeResize(sim.index + sim.nDoF);
    sim_parameters->q_dot0.conservativeResize(sim.index + sim.nDoF);

    // Set initial positions and velocities
    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->q0[sim.index + i] = vertices[i];
        sim_parameters->q_dot0[sim.index + i] = 0.0;
    }

    // Resize parameter containers
    sim_parameters->p.conservativeResize(sim.parameter_index + sim.nParameters);

    // Fill parameters
    Parameter stiffness = Parameter {.value=stiffness_value, .index=sim.parameter_index};
    Parameter bend_stiffness = Parameter {.value=bend_stiffness_value, .index=sim.parameter_index+1};
    Parameter rest_length = Parameter {.value=0, .index=sim.parameter_index+2};
    sim_parameters->p[stiffness.index] = stiffness.value;
    sim_parameters->p[bend_stiffness.index] = bend_stiffness.value;
    sim_parameters->p[rest_length.index] = rest_length.value; // garbage


    // Set up the springs
    std::vector<Edge> internalEdges, externalEdges;
    mesh_boundary(vertices, indices, internalEdges, externalEdges);

    unsigned int n_flex = internalEdges.size() / 2.0 + externalEdges.size();
    unsigned int n_bend = internalEdges.size() / 2.0;

    for (size_t i = 0; i < externalEdges.size(); i++) {
        Edge &e = externalEdges[i];
        rest_length.value = distance(vertices, e.a, e.b);
        interaction_manager->m_springs.push_back(
            Spring(sim.index + 3*e.a, sim.index + 3*e.b, stiffness, rest_length));
    }

    for (size_t i = 0; i < internalEdges.size(); i += 2) {
        Edge &e1 = internalEdges[i];
        Edge &e2 = internalEdges[i + 1];
        rest_length.value = distance(vertices, e1.a, e1.b);
        // Normal spring
        interaction_manager->m_springs.push_back(
            Spring(sim.index + 3*e1.a, sim.index + 3*e1.b, stiffness, rest_length));

        // Bend spring
        rest_length.value = distance(vertices, e1.opposite, e2.opposite);
        interaction_manager->m_springs.push_back(
            Spring(sim.index + 3*e1.opposite, sim.index + 3*e2.opposite, bend_stiffness, rest_length));
    }

    // Add gravity
    Vec3 gravity_vec = Vec3(0, -1, 0);
    for (size_t i=0; i < sim.nDoF; i+=3) {
        interaction_manager->m_gravity.push_back(Gravity(sim.index + i, gravity_vec));
    }

    return sim;
}

/*
 * Generates a simulable with a differentiable parameters for each spring.
 * Stiffnes values -> a value for each flex spring
 * Bend Stiffnes values -> a value for each bend spring
 */
Simulable generate_mass_spring(Simulation* simulation,
                          const std::vector<Scalar>& vertices,
                          const std::vector<unsigned int>& indices,
                          Scalar node_mass,
                          const std::vector<Scalar>& stiffness_values,
                          const std::vector<Scalar>& bend_stiffness_values)
{
    assert(vertices.size() % 3 == 0);

    SimulationParameters* sim_parameters = &simulation->simulation_parameters;
    InteractionManager* interaction_manager = &simulation->interaction_manager;
    const Simulable sim {
        .index = static_cast<unsigned int>(sim_parameters->q0.size()),
        .nDoF = static_cast<unsigned int>(vertices.size()),
        .parameter_index = static_cast<unsigned int>(sim_parameters->p.size()),
        .nParameters = static_cast<unsigned int>(stiffness_values.size() + bend_stiffness_values.size()+1)
    };

    // Set up mass matrix
    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->mass.push_back(Triplet(i+sim.index, i+sim.index, node_mass));
    }

    // Resize the DOF containers
    sim_parameters->q0.conservativeResize(sim.index + sim.nDoF);
    sim_parameters->q_dot0.conservativeResize(sim.index + sim.nDoF);

    // Set initial positions and velocities
    for (size_t i=0; i < sim.nDoF; i++) {
        sim_parameters->q0[sim.index + i] = vertices[i];
        sim_parameters->q_dot0[sim.index + i] = 0.0;
    }

    // Resize parameter containers
    sim_parameters->p.conservativeResize(sim.parameter_index + sim.nParameters);

    // Garbage rest_length
    Parameter rest_length = Parameter {.value=0, .index=sim.parameter_index + sim.nParameters-1};
    sim_parameters->p[rest_length.index] = rest_length.value; // garbage


    // Set up the springs
    std::vector<Edge> internalEdges, externalEdges;
    mesh_boundary(vertices, indices, internalEdges, externalEdges);

    unsigned int n_flex = internalEdges.size() / 2.0 + externalEdges.size();
    unsigned int n_bend = internalEdges.size() / 2.0;

    unsigned int flex_index = 0;
    unsigned int bend_index = 0;
    for (size_t i = 0; i < externalEdges.size(); i++) {
        Edge &e = externalEdges[i];
        rest_length.value = distance(vertices, e.a, e.b);
        {
            // Normal spring
            Parameter stiffness = create_parameter(sim_parameters,
                                                   stiffness_values[flex_index],
                                                   sim.parameter_index+flex_index);
            flex_index++;
            interaction_manager->m_springs.push_back(
                Spring(sim.index + 3*e.a, sim.index + 3*e.b, stiffness, rest_length));
        }
    }

    for (size_t i = 0; i < internalEdges.size(); i += 2) {
        Edge &e1 = internalEdges[i];
        Edge &e2 = internalEdges[i + 1];
        rest_length.value = distance(vertices, e1.a, e1.b);
        {
            // Normal spring
            Parameter stiffness = create_parameter(sim_parameters,
                                                   stiffness_values[flex_index],
                                                   sim.parameter_index+flex_index);
            flex_index++;
            interaction_manager->m_springs.push_back(
                Spring(sim.index + 3*e1.a, sim.index + 3*e1.b, stiffness, rest_length));
        }

        {
            // Bend spring
            Parameter bend_stiffness = create_parameter(sim_parameters,
                                                        bend_stiffness_values[bend_index],
                                                        sim.parameter_index+n_flex+bend_index);
            bend_index++;
            rest_length.value = distance(vertices, e1.opposite, e2.opposite);
            interaction_manager->m_springs.push_back(
                Spring(sim.index + 3*e1.opposite, sim.index + 3*e2.opposite, bend_stiffness, rest_length));
        }
    }

    // Add gravity
    Vec3 gravity_vec = Vec3(0, -1, 0);
    for (size_t i=0; i < sim.nDoF; i+=3) {
        interaction_manager->m_gravity.push_back(Gravity(sim.index + i, gravity_vec));
    }

    return sim;
}

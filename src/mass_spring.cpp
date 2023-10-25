#include <cassert>
#include <memory>

#include "mesh_boundary.hpp"
#include "gravity.hpp"
#include "interaction_manager.hpp"
#include "linear_algebra.hpp"
#include "simulable.hpp"
#include "simulable_generators.hpp"
#include "simulation.hpp"
#include "simulation_parameters.hpp"
#include "spring.hpp"

/* Computes euclidean distance between vertex i1 and vertex i2. */
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

    Parameter stiffness = create_parameter(sim_parameters, stiffness_value, sim.parameter_index);
    Parameter bend_stiffness = create_parameter(sim_parameters, bend_stiffness_value, sim.parameter_index+1);
    Parameter rest_length = create_parameter(sim_parameters, 0, sim.parameter_index+2); // garbage


    // Set up the springs
    std::vector<Edge> internalEdges, externalEdges;
    mesh_boundary(vertices, indices, internalEdges, externalEdges);

    unsigned int n_flex = internalEdges.size() / 2.0 + externalEdges.size();
    unsigned int n_bend = internalEdges.size() / 2.0;
    // std::cout << "n_tension = " << n_flex << std::endl;
    // std::cout << "n_bend = " << n_bend << std::endl;
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
    Vec3 gravity_vec = node_mass * Vec3(0, -1, 0);
    for (size_t i=0; i < sim.nDoF; i+=3) {
        interaction_manager->m_gravity.push_back(Gravity(sim.index + i, gravity_vec));
    }

    simulation->simulables.push_back(sim);
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
    Parameter rest_length = create_parameter(sim_parameters, 0, sim.parameter_index+sim.nParameters-1);

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
    Vec3 gravity_vec = node_mass * Vec3(0, -1, 0);
    for (size_t i=0; i < sim.nDoF; i+=3) {
        interaction_manager->m_gravity.push_back(Gravity(sim.index + i, gravity_vec));
    }

    // for (int i = 0; i < sim_parameters->p.size(); i++) {
    //     std::cout << "Value p("<< i << ") = "  << sim_parameters->p[i] << std::endl;
    // }

    simulation->simulables.push_back(sim);
    return sim;
}

void calculate_initial_state_jacobian_tilt_angle(SimulationParameters* parameters, const Simulable& sim, const Parameter& tilt_angle) {
    const Vector& x = parameters->q0;
    const Vec3 origin = Vec3(x[0], x[1], x[2]);
    for (size_t i = sim.index; i < sim.index+sim.nDoF; i+=3) {
        const Vec3 point = Vec3(x[i], x[i+1], x[i+2]);
        const Vec3 delta = point - origin;
        const Scalar radians_to_degrees_jacobian = PI/180;
        // x coordinate not affected by rotation
        // y coordinate
        parameters->dq0_dp_triplets.push_back(Triplet(i+1, tilt_angle.index, - delta.z() * radians_to_degrees_jacobian));

        // z coordinate
        parameters->dq0_dp_triplets.push_back(Triplet(i+2, tilt_angle.index, delta.y() * radians_to_degrees_jacobian));
    }
}

Parameter add_tilt_angle_parameter(SimulationParameters* parameters, const Simulable& sim, Scalar tilt_angle) {
    unsigned int index = parameters->p.size();
    parameters->p.conservativeResize(index+1);
    Parameter tilt_angle_parameter = create_parameter(parameters, tilt_angle, index);
    calculate_initial_state_jacobian_tilt_angle(parameters, sim, tilt_angle_parameter);
    return tilt_angle_parameter;
}

void add_initial_velocity_parameters(SimulationParameters* parameters, const Simulable& sim, std::vector<Parameter>* out_initial_velocities) {
    const unsigned int index = sim.index;
    const unsigned int pindex = parameters->p.size();
    const unsigned int nDoF = sim.nDoF;

    parameters->p.conservativeResize(pindex+nDoF);
    out_initial_velocities->resize(nDoF);

    for (unsigned int i=0; i < nDoF; i++) {
        // Create a paramter for each initial velocity
        (*out_initial_velocities)[i] = create_parameter(parameters, parameters->q_dot0[index + i], pindex+i);

        // Compute the initial conditions jacobian: in this case is a simple identity as the initial
        // velocities themselves are the parameters.
        //
        // std::cout << "Pindex= " << pindex << std::endl;
        parameters->dq_dot0_dp_triplets.push_back(Triplet(index+i, pindex + index+i, 1));
    }
}

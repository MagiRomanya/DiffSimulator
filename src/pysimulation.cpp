#include <raylib.h>
#include <unistd.h>
#include <vector>

#include "contact.hpp"
#include "mesh_boundary.hpp"
#include "pysimulation.hpp"
#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "raymath.h"
#include "simulable.hpp"
#include "simulable_generators.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"

PySimulation::~PySimulation() {
    if (graphics) {
        CloseWindow();
    }
}

PySimulation::PySimulation()
{
    set_up_simulation();
}

PySimulation::PySimulation(Scalar k, Scalar k_bend, bool graphics)
{
    this->graphics = graphics;
    set_up_simulation();
    reset_simulation(k, k_bend);
}

PySimulation::PySimulation(Scalar k, Scalar k_bend, Scalar tilt_angle, bool graphics) {
    this->graphics = graphics;
    set_up_simulation();
    reset_simulation(k, k_bend, tilt_angle);
}

PySimulation::PySimulation(std::vector<Scalar> k, std::vector<Scalar> k_bend, bool graphics)
{

    this->graphics = graphics;
    set_up_simulation();
    reset_simulation(k, k_bend);
}

void PySimulation::fill_containers() {
    const unsigned int nDoF = getDoF();
    const unsigned int nParameters = simulation.simulation_parameters.p.size();
    const Scalar h = simulation.simulation_parameters.TimeStep;

    // Clear the last energy derivatives object
    f = EnergyDerivatives(nDoF, nParameters);
    // Calculate the energy derivatives
    simulation.interaction_manager.calculate_energy_derivatives(&state, &f);
    std::vector<ContactData> contacts;
    simulation.contact_manager.find_contacts(state, contacts);
    simulation.contact_manager.compute_contacts_energy_derivatives(contacts, &f);

    // Construct sparse matrices & delete previous ones
    df_dx.setFromTriplets(f.df_dx_triplets.begin(), f.df_dx_triplets.end());
    df_dv.setFromTriplets(f.df_dv_triplets.begin(), f.df_dv_triplets.end());

    // Handle frozen particles
    handle_frozen_dof(simulation.simulation_parameters.frozen_dof, &df_dx);
    handle_frozen_dof(simulation.simulation_parameters.frozen_dof, &df_dv);
    handle_frozen_dof(simulation.simulation_parameters.frozen_dof, &f.force);

    // Construct the system of equations
    equation_vector = h * (f.force + h * df_dx * state.q_dot);
    equation_matrix = mass_matrix - h * df_dv - h * h * df_dx;

    handle_frozen_dof(simulation.simulation_parameters.frozen_dof, &equation_matrix);
    handle_frozen_dof(simulation.simulation_parameters.frozen_dof, &equation_vector);
}

void PySimulation::set_state(Vector xi, Vector vi) {
    state.q = xi;
    state.q_dot = vi;
}

SparseMatrix PySimulation::getEquationMatrix() { return equation_matrix; }

Vector PySimulation::getEquationVector() { return equation_vector; }

Vector PySimulation::getForce() { return f.force; }

Vector PySimulation::getPosition() { return  state.q; }

Vector PySimulation::getVelocity() { return  state.q_dot; }

Vector PySimulation::getDiffParameteres() { return simulation.simulation_parameters.p; }

SparseMatrix PySimulation::getMassMatrix() { return mass_matrix; }

Mat PySimulation::getParameterJacobian() {
    return f.df_dp;
}

SparseMatrix PySimulation::getForcePositionJacobian() { return df_dx; }

int PySimulation::getDoF() { return simulation.simulation_parameters.q0.size();}

Scalar PySimulation::getTimeStep() { return simulation.simulation_parameters.TimeStep; }

std::array<unsigned int, 2> PySimulation::getGridDimensions() { return {grid_n, grid_m}; }

std::array<unsigned int, 2> PySimulation::getNumberOfSprings() { return {n_flex, n_bend}; }

void PySimulation::render_state() {
    if (not graphics) return;
    // Cameara && inputs
    {
        UpdateCamera(&camera, CAMERA_FREE);
        if (IsKeyDown(KEY_Z)) camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    }

    // Update mesh GPU data
    Vector mesh_data = state.q.segment(mass_spring.index, mass_spring.index+mass_spring.nDoF);
    UpdateMeshBuffer(cloth_mesh, 0, mesh_data.data(), mesh_data.size()*sizeof(float), 0);
    //----------------------------------------------------------------------------------
    // Draw
    //----------------------------------------------------------------------------------
    BeginDrawing();
    {

        ClearBackground(RAYWHITE);
        DrawFPS(50, 50);
        BeginMode3D(camera);
        {
                DrawMesh(cloth_mesh, cloth_material, MatrixIdentity());
                DrawGrid(100, 1.0f);
                for (size_t i = 0; i < simulation.contact_manager.sphere_colliders.size(); i++) {
                    const Sphere& s = simulation.contact_manager.sphere_colliders[i];
                    DrawSphere(Vector3(s.center.x(), s.center.y(), s.center.z()), s.radius, GREEN);
                }
        }
        EndMode3D();
    }
    EndDrawing();
    //----------------------------------------------------------------------------------
}

SparseMatrix PySimulation::getInitialPositionJacobian() {
    const unsigned int nDoF = getDoF();
    const unsigned int nParameters = simulation.simulation_parameters.p.size();
    SparseMatrix dx0_dp = SparseMatrix(nDoF, nParameters);
    dx0_dp.setFromTriplets(simulation.simulation_parameters.dq0_dp_triplets.begin(),
                           simulation.simulation_parameters.dq0_dp_triplets.end());

    return dx0_dp;
}

SparseMatrix PySimulation::getInitialVelocityJacobian() {
    const unsigned int nDoF = getDoF();
    const unsigned int nParameters = simulation.simulation_parameters.p.size();
    SparseMatrix dv0_dp = SparseMatrix(nDoF, nParameters);
    dv0_dp.setFromTriplets(simulation.simulation_parameters.dq_dot0_dp_triplets.begin(),
                           simulation.simulation_parameters.dq_dot0_dp_triplets.end());
    return dv0_dp;
}

void PySimulation::reset_simulation(Scalar stiffness, Scalar bend_stiffness) {
    Simulation sim;
    mass_spring = generate_mass_spring(&sim, vertices, indices, node_mass, stiffness, bend_stiffness);

    // Fix the top 2 corners of the cloth
    sim.simulation_parameters.frozen_dof.push_back(0);
    sim.simulation_parameters.frozen_dof.push_back(1);
    sim.simulation_parameters.frozen_dof.push_back(2);

    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3);
    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+1);
    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+2);

    state = sim.getInitialState();

    simulation = sim;

    // Add a sphere collider
#ifdef ENABLE_CONTACT
    Sphere sphere = {Vec3(0,0,0), 2};
    simulation.contact_manager.sphere_colliders.push_back(sphere);
#endif
}

void PySimulation::reset_simulation(Scalar stiffness, Scalar bend_stiffness, Scalar tilt_angle) {
    // ANGLE IN RADIANTS!!
    angle=tilt_angle;
    create_mesh();
    Simulation sim;
    mass_spring = generate_mass_spring(&sim, vertices, indices, node_mass, stiffness, bend_stiffness);

    // Fix the top 2 corners of the cloth
    sim.simulation_parameters.frozen_dof.push_back(0);
    sim.simulation_parameters.frozen_dof.push_back(1);
    sim.simulation_parameters.frozen_dof.push_back(2);

    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3);
    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+1);
    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+2);

    state = sim.getInitialState();

    simulation = sim;

    add_tilt_angle_parameter(&simulation.simulation_parameters, mass_spring, tilt_angle);

    // Add a sphere collider
#ifdef ENABLE_CONTACT
    Sphere sphere = {Vec3(0,0,0), 2};
    simulation.contact_manager.sphere_colliders.push_back(sphere);
#endif
}

void PySimulation::reset_simulation(std::vector<Scalar> stiffness, std::vector<Scalar> bend_stiffness) {
    Simulation sim;
    mass_spring = generate_mass_spring(&sim, vertices, indices, node_mass, stiffness, bend_stiffness);

    // Fix the top 2 corners of the cloth
    sim.simulation_parameters.frozen_dof.push_back(0);
    sim.simulation_parameters.frozen_dof.push_back(1);
    sim.simulation_parameters.frozen_dof.push_back(2);

    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3);
    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+1);
    sim.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+2);

    state = sim.getInitialState();
    simulation = sim;

    // Add a sphere collider
    Sphere sphere = {Vec3(0,0,0), 2};
    simulation.contact_manager.sphere_colliders.push_back(sphere);
}

void PySimulation::set_up_simulation() {
    const int screenWidth = 800*2;
    const int screenHeight = 450*2;
    if (graphics) {
        InitWindow(screenWidth, screenHeight, "Simulator");
        camera = create_camera();
    }
    //--------------------------------------------------------------------------------------

    const Scalar stiffness = 100.0;
    const Scalar bend_stiffness = stiffness / 100;

    // Create the mesh
    //--------------------------------------------------------------------------------------
    create_mesh();
    //--------------------------------------------------------------------------------------

    // Create the simulable
    mass_spring = generate_mass_spring(&simulation, vertices, indices, node_mass, stiffness, bend_stiffness);

    // Count how many springs of each type there are
    auto nSprings = count_springs(vertices, indices);
    n_flex = nSprings[0];
    n_bend = nSprings[1];

    // Fix the top 2 corners of the cloth
    simulation.simulation_parameters.frozen_dof.push_back(0);
    simulation.simulation_parameters.frozen_dof.push_back(1);
    simulation.simulation_parameters.frozen_dof.push_back(2);

    simulation.simulation_parameters.frozen_dof.push_back((grid_n-1)*3);
    simulation.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+1);
    simulation.simulation_parameters.frozen_dof.push_back((grid_n-1)*3+2);
    //--------------------------------------------------------------------------------------
    // Create the mass matrix
    resize_containers();
    mass_matrix.setFromTriplets(simulation.simulation_parameters.mass.begin(),
                                simulation.simulation_parameters.mass.end());

    state = simulation.getInitialState();

    // Add a sphere collider
    Sphere sphere = {Vec3(0,0,0), 2};
    simulation.contact_manager.sphere_colliders.push_back(sphere);
}

void PySimulation::resize_containers() {
    const unsigned int nDoF = getDoF();
    const unsigned int nParameters = simulation.simulation_parameters.p.size();
    df_dx.resize(nDoF, nDoF);
    df_dv.resize(nDoF, nDoF);
    mass_matrix.resize(nDoF, nDoF);
    equation_matrix.resize(nDoF, nDoF);
}


void PySimulation::create_mesh() {
    const unsigned int grid_node_width = GRID_NODE_SIDE;
    grid_n = grid_node_width;
    grid_m = grid_node_width;
    const float grid_width = GRID_WIDTH_LENGTH;
    cloth_mesh = GenMeshPlaneNoGPU(grid_width, grid_width, grid_node_width-1, grid_node_width-1);
    if (graphics) {
        UploadMesh(&cloth_mesh, true);
        const Texture2D cloth_texture = LoadTexture("../resources/warning.png");
        cloth_material = LoadMaterialDefault();
        SetMaterialTexture(&cloth_material, 0, cloth_texture);
    }

    const unsigned int nDoF = cloth_mesh.vertexCount * 3;
    const unsigned int n_indices = cloth_mesh.triangleCount * 3;
    vertices = std::vector<Scalar>(cloth_mesh.vertices, cloth_mesh.vertices + nDoF);
    indices = std::vector<unsigned int>(cloth_mesh.indices, cloth_mesh.indices + n_indices);

    // Reposition the mesh in the world
    translate_vertices(vertices, Vec3(0, 0, GRID_WIDTH_LENGTH/2));
    rotate_vertices_arround_axis(vertices, Vec3(angle, 0, 0));
    translate_vertices(vertices, Vec3(0, grid_width*2, 0));
    cloth_mesh.vertices = vertices.data();
}

#include <raylib.h>
#include <vector>

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
    CloseWindow();
}

PySimulation::PySimulation()
{
    set_up_simulation();
}

PySimulation::PySimulation(Scalar k, Scalar k_bend)
{
    set_up_simulation();
    reset_simulation(k, k_bend);
}

PySimulation::PySimulation(std::vector<Scalar> k, std::vector<Scalar> k_bend)
{
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

std::vector<unsigned int> PySimulation::getSpringNodeIndices() {}

std::vector<unsigned int> PySimulation::getBendSpringNodeIndices() {}

std::array<unsigned int, 2> PySimulation::getGridDimensions() { return {grid_n, grid_m}; }

std::array<unsigned int, 2> PySimulation::getNumberOfSprings() { return {n_flex, n_bend}; }

void PySimulation::render_state() {
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
        }
        EndMode3D();
    }
    EndDrawing();
    //----------------------------------------------------------------------------------
}

SparseMatrix PySimulation::getInitialPositionJacobian() {
    const unsigned int nDoF = getDoF();
    const unsigned int nParameters = simulation.simulation_parameters.p.size();
    return SparseMatrix(nDoF, nParameters);
}

SparseMatrix PySimulation::getInitialVelocityJacobian() {
    const unsigned int nDoF = getDoF();
    const unsigned int nParameters = simulation.simulation_parameters.p.size();
    return SparseMatrix(nDoF, nParameters);
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
}

void PySimulation::set_up_simulation() {
    const int screenWidth = 800*2;
    const int screenHeight = 450*2;
    InitWindow(screenWidth, screenHeight, "Simulator");
    //--------------------------------------------------------------------------------------

    const Scalar stiffness = 100.0;
    const Scalar bend_stiffness = stiffness / 100;

    // Create the mesh
    //--------------------------------------------------------------------------------------
    const unsigned int grid_node_width = 20;
    grid_n = grid_node_width;
    grid_m = grid_node_width;
    const float grid_width = 5.0f;
    cloth_mesh = GenMeshPlane(grid_width, grid_width, grid_node_width-1, grid_node_width-1);
    const Texture2D cloth_texture = LoadTexture("../resources/warning.png");
    cloth_material = LoadMaterialDefault();
    SetMaterialTexture(&cloth_material, 0, cloth_texture);
    //--------------------------------------------------------------------------------------

    // Creating the simulable from the mesh
    //--------------------------------------------------------------------------------------
    const unsigned int nDoF = cloth_mesh.vertexCount * 3;
    const unsigned int n_indices = cloth_mesh.triangleCount * 3;
    vertices = std::vector<Scalar>(cloth_mesh.vertices, cloth_mesh.vertices + nDoF);
    indices = std::vector<unsigned int>(cloth_mesh.indices, cloth_mesh.indices + n_indices);

    // Reposition the mesh in the world
    // rotate_vertices_arround_axis(vertices, Vec3(PI/2, 0, 0));
    translate_vertices(vertices, Vec3(0, grid_width*1.3, 0));
    cloth_mesh.vertices = vertices.data();

    // Create the simulable
    mass_spring = generate_mass_spring(&simulation, vertices, indices, node_mass, stiffness, bend_stiffness);

    // Count how many springs of each type there are
    std::array<unsigned int, 2> nSprings = count_springs(vertices, indices);
    n_flex = nSprings[0];
    n_bend = nSprings[1];

    // Fix the top 2 corners of the cloth
    simulation.simulation_parameters.frozen_dof.push_back(0);
    simulation.simulation_parameters.frozen_dof.push_back(1);
    simulation.simulation_parameters.frozen_dof.push_back(2);

    simulation.simulation_parameters.frozen_dof.push_back((grid_node_width-1)*3);
    simulation.simulation_parameters.frozen_dof.push_back((grid_node_width-1)*3+1);
    simulation.simulation_parameters.frozen_dof.push_back((grid_node_width-1)*3+2);
    //--------------------------------------------------------------------------------------
    // Create the mass matrix
    resize_containers();
    mass_matrix.setFromTriplets(simulation.simulation_parameters.mass.begin(),
                                simulation.simulation_parameters.mass.end());

    state = simulation.getInitialState();
}

void PySimulation::resize_containers() {
    const unsigned int nDoF = getDoF();
    const unsigned int nParameters = simulation.simulation_parameters.p.size();
    df_dx.resize(nDoF, nDoF);
    df_dv.resize(nDoF, nDoF);
    mass_matrix.resize(nDoF, nDoF);
    equation_matrix.resize(nDoF, nDoF);
}

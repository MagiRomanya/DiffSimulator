#include <cstdio>
#include <iostream>
#include <vector>
#include <raylib.h>
#include <raymath.h>


#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "simulable.hpp"
#include "simulable_generators.hpp"
#include "simulation.hpp"

Camera3D create_camera();

int main() {
    // Window Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 800*2;
    const int screenHeight = 450*2;
    InitWindow(screenWidth, screenHeight, "Simulator");

    //--------------------------------------------------------------------------------------

    Simulation simulation;

    const Scalar mass = 1.0;
    const Scalar stiffness = 1.0;
    const Scalar bend_stiffness = 0.0;

    // Create the mesh
    const unsigned int grid_width = 2;
    Mesh cloth_mesh = GenMeshPlane(10.0f, 10.0f, grid_width-1, grid_width-1);
    const unsigned int nDoF = cloth_mesh.vertexCount * 3;
    const unsigned int n_indices = cloth_mesh.triangleCount * 3;
    std::vector<Scalar> vertices(cloth_mesh.vertices, cloth_mesh.vertices + nDoF);
    std::vector<unsigned int> indices(cloth_mesh.indices, cloth_mesh.indices + n_indices);

    Simulable mass_spring = generate_mass_spring(&simulation, vertices, indices, mass, stiffness, bend_stiffness);

    // generate_mass_spring_pair(&simulation, mass, stiffness);

    simulation.simulation_parameters.frozen_dof.push_back(0);
    simulation.simulation_parameters.frozen_dof.push_back(1);
    simulation.simulation_parameters.frozen_dof.push_back(2);

    int offset = -1;
    simulation.simulation_parameters.frozen_dof.push_back(grid_width*3+offset);
    simulation.simulation_parameters.frozen_dof.push_back(grid_width*3+offset+1);
    simulation.simulation_parameters.frozen_dof.push_back(grid_width*3+offset+2);
    PhysicsState state = simulation.getInitialState();

    Camera3D camera = create_camera();

    // for (unsigned int p=0; p < state.q.size(); p+=3) {
    //     std::cout << "Particle "<< p <<" -> "<< get_particle_position(&state, p).transpose() << std::endl;
    // }

    // Main game loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        // Cameara && inputs
        UpdateCamera(&camera, CAMERA_FREE);
        if (IsKeyDown(KEY_Z)) camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
        if (IsKeyDown(KEY_Q)) break;

        // Update physics
        std::cout << "Particle 1 -> "<< get_particle_position(&state, 0).transpose() << std::endl;
        std::cout << "Particle 2 -> "<< get_particle_position(&state, 3*grid_width+offset).transpose() << std::endl;

        simulation.step(&state);

        // Update mesh
        Vector mesh_data = state.q.segment(mass_spring.index, mass_spring.index+mass_spring.nDoF);
        UpdateMeshBuffer(cloth_mesh, 0, mesh_data.data(), mesh_data.size()*sizeof(float), 0);
        //----------------------------------------------------------------------------------
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);

                DrawMesh(cloth_mesh, LoadMaterialDefault(), MatrixIdentity());
                DrawGrid(100, 1.0f);

            EndMode3D();
            for (unsigned int i=0; i<state.q.size(); i++) {
                // const unsigned int width = 20;
                const unsigned int width = GetScreenWidth() / state.q.size();
                const Scalar percent = ((Scalar)i/(Scalar)state.q.size());
                const Color color = ColorFromHSV(360*percent, 1, 1);
                DrawRectangle(width*i, 0, width, (int) (state.q(i)*10), color);
            }
        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}

Camera3D create_camera() {
    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 2.0f, 10.0f }; // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type

    // DisableCursor();                    // Limit cursor to relative movement inside the window

    SetTargetFPS(60);                   // Set our game to run at 60 frames-per-second
    return camera;
}

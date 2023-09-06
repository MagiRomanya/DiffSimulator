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
void translate_vertices(std::vector<Scalar>& vert, const Vec3& translation);

int main() {
    // Window Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 800*2;
    const int screenHeight = 450*2;
    InitWindow(screenWidth, screenHeight, "Simulator");


    //--------------------------------------------------------------------------------------

    Simulation simulation;

    const Scalar mass = 1.0;
    const Scalar stiffness = 100.0;
    const Scalar bend_stiffness = stiffness / 100;

    // Create the mesh
    const unsigned int grid_node_width = 20;
    const float grid_width = 5.0f;
    Mesh cloth_mesh = GenMeshPlane(grid_width, grid_width, grid_node_width-1, grid_node_width-1);
    const Texture2D cloth_texture = LoadTexture("resources/warning.png");
    Material cloth_material = LoadMaterialDefault();
    SetMaterialTexture(&cloth_material, 0, cloth_texture);

    const unsigned int nDoF = cloth_mesh.vertexCount * 3;
    const unsigned int n_indices = cloth_mesh.triangleCount * 3;
    std::vector<Scalar> vertices(cloth_mesh.vertices, cloth_mesh.vertices + nDoF);
    std::vector<unsigned int> indices(cloth_mesh.indices, cloth_mesh.indices + n_indices);
    translate_vertices(vertices, Vec3(0, grid_width*1.3, 0));

    Simulable mass_spring = generate_mass_spring(&simulation, vertices, indices, mass, stiffness, bend_stiffness);

    // generate_mass_spring_pair(&simulation, mass, stiffness);

    // Fix the top 2 corners of the cloth
    simulation.simulation_parameters.frozen_dof.push_back(0);
    simulation.simulation_parameters.frozen_dof.push_back(1);
    simulation.simulation_parameters.frozen_dof.push_back(2);

    simulation.simulation_parameters.frozen_dof.push_back((grid_node_width-1)*3);
    simulation.simulation_parameters.frozen_dof.push_back((grid_node_width-1)*3+1);
    simulation.simulation_parameters.frozen_dof.push_back((grid_node_width-1)*3+2);

    PhysicsState state = simulation.getInitialState();
    Camera3D camera = create_camera();

    bool game_paused = true;

    // Main game loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        // Cameara && inputs
        {
            UpdateCamera(&camera, CAMERA_FREE);
            if (IsKeyDown(KEY_Z)) camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
            // Quit
            if (IsKeyPressed(KEY_Q)) break;
            // Pause the simulation
            if (IsKeyPressed(KEY_P)) game_paused = not game_paused;
            // Restart the simulation
            if (IsKeyPressed(KEY_R)) {
                state = simulation.getInitialState();
            }
        }

        // Update physics
        if (!game_paused) {
            simulation.step(&state);
        }

        // Update mesh
        Vector mesh_data = state.q.segment(mass_spring.index, mass_spring.index+mass_spring.nDoF);
        UpdateMeshBuffer(cloth_mesh, 0, mesh_data.data(), mesh_data.size()*sizeof(float), 0);
        //----------------------------------------------------------------------------------
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
        {

            ClearBackground(RAYWHITE);

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

    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}

Camera3D create_camera() {
    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 2.0f, 15.0f }; // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type

    // DisableCursor();                    // Limit cursor to relative movement inside the window

    SetTargetFPS(60);                   // Set our game to run at 60 frames-per-second
    return camera;
}

void translate_vertices(std::vector<Scalar>& vert, const Vec3& translation) {
    for (unsigned int i = 0; i < vert.size(); i+=3) {
        Vec3 new_pos = Vec3(vert[i], vert[i+1], vert[i+2]) + translation;
        vert[i] = new_pos.x();
        vert[i+1] = new_pos.y();
        vert[i+2] = new_pos.z();
    }
}

#include "utility_functions.hpp"

#include <raymath.h>
#include <rlgl.h>



Camera3D create_camera(unsigned int FPS) {
    // Disable backface culling
    rlDisableBackfaceCulling();
    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 2.0f, 15.0f }; // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera projection type

    // DisableCursor();                    // Limit cursor to relative movement inside the window

    SetTargetFPS(FPS);                   // Set our game to run at 60 frames-per-second
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

void rotate_vertices_arround_axis(std::vector<Scalar>& vert, const Vec3& rotation) {
    const Scalar angle = rotation.norm();
    Vector3 axis = Vector3{rotation.x() / angle,
                           rotation.y() / angle,
                           rotation.z() / angle};

    for (unsigned int i = 0; i < vert.size(); i+=3) {
        const Vector3 pos = Vector3{vert[i], vert[i+1], vert[i+2]};
        Vector3 new_pos = Vector3RotateByAxisAngle(pos, axis, angle);
        vert[i] = new_pos.x;
        vert[i+1] = new_pos.y;
        vert[i+2] = new_pos.z;
    }

}

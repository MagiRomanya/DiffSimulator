#include "contact.hpp"
#include "physics_state.hpp"

double Sphere::distance_point(const Vec3& point, bool& valid) const {
    /* Distance between a point and a sphere.
     * The result is positive when the point is outside, zero when the point is on the surface,
     * and is negative when the point is inside the sphere. */

    valid = true;
    const Vec3 dv = point - center;
    return dv.norm() - r;
}

Vec3 Sphere::outward_direction(const Vec3& point) const {
    return (point - center).normalized();
}

double InfPlane::distance_point(const Vec3& point, bool& valid) const {
    /* Distance between a point and a plane.
     * The result is positive when the point is in the region of space where the
     * normal of the plane points to the point, zero when the point is on the surface of the
     * plane, and it is negative if the normal points to the opposite direction of where the point is */
    valid = true;
    return normal.dot(point - center);
}

Vec3 InfPlane::outward_direction(const Vec3 &point) const {
    return normal;
}

double FinPlane::distance_point(const Vec3& point,  bool& valid) const {
    /* Returns the distance between a finite plane and a point.
     * The valid bool means weather or not the point is inside the finite
     * plane domain or, on the contrary, is outside of it. */

    const double d = normal.dot(point - center);
    const Vec3& intersection = point - d * normal;
    const Vec3& tangent = up.cross3(normal).normalized();
    const Vec3& bitangent = normal.cross3(tangent).normalized();
    if (abs((intersection - point).dot(tangent)) > radius) {
        valid = false;
        return 0.0f;
    }
    if (abs((intersection - point).dot(bitangent)) > radius) {
        valid = false;
        return 0.0f;
    }
    valid = true;
    return d;
}

Vec3 FinPlane::outward_direction(const Vec3 &point) const {
    return normal;
}

Contact::Contact(const Sphere* sphere) {
    geometry_type = SPHERE;
    geometry = sphere;
}

Contact::Contact(const InfPlane* plane) {
    geometry_type = INFPLANE;
    geometry = plane;
}

Contact::Contact(const FinPlane* plane) {
    geometry_type = FINPLANE;
    geometry = plane;
}

Contact::~Contact(){

}

void Contact::apply(const PhysicsState* state, EnergyDerivatives* f) {
    for (int i = 0; i < sys->get_n_particles(); i++) {
        ///////////// COLLISION DETECTION ///////////////
        unsigned int index; // Particle index
        bool valid = true;
        Vec3 point = get_particle_position(state, index);
        double dist = geometry->distance_point(point, valid);

        if (!valid) continue;
        if (dist > 0) continue;

        /////////////// COLLISION RESPONSE ////////////////
        Vec3 normal_to_surface = geometry->outward_direction(point);
        Vec3 contact_force = force(normal_to_surface, dist);
        f->force(index) = contact_force.x();
        f->force(index+1) = contact_force.y();
        f->force(index+2) = contact_force.z();

        Eigen::Matrix3d df_dx = force_derivative(state, normal_to_surface, dist);
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++) {
                f->df_dx_triplets.push_back(Triplet(index+j, index+k, df_dx(j,k)));
            }
        }
    }

}

Vec3 Contact::force(const Vec3& direction, const double dist) {
    return - contact_stiffness * dist * direction;
}

Eigen::Matrix3d Contact::force_derivative(const PhysicsState* state, const Vec3& direction, const double dist) {
    return - contact_stiffness * direction * direction.transpose()
}

#ifndef CONTACT_H_
#define CONTACT_H_

#include <raylib.h>

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct ContactData {
    Vec3   normal;              // Normal vector pointing OUT of the collision geometry
    Scalar s_distance;          // Signed distance (negative inside collider)
    Scalar stiffness = 10000;    // Contact stiffness
    unsigned int index;         // Index of the DoF it affects in the PhysicsState vectors
};

struct Sphere {
    Vec3 center;
    Scalar radius;
};

struct Plane {
    Vec3 center;
    Vec3 normal;
};

struct ContactManager {
    std::vector<Sphere> sphere_colliders;
    std::vector<Plane> plane_colliders;
    // TODO AABB, more shapes ...

    void find_contacts(const PhysicsState& state, std::vector<ContactData>& contacts);

    void compute_contacts_energy_derivatives(const std::vector<ContactData>& contacts, EnergyDerivatives* f);
};


void compute_sphere_contact_geometry(const Sphere& sphere, const Vec3& point, ContactData* out);

void compute_plane_contact_geometry(const Plane& plane, const Vec3& point, ContactData* out);

Vec3 compute_contact_force(const ContactData& contact);

Mat3 compute_contact_force_position_jacobian(const ContactData& contact);

void compute_contact_energy_derivatives(const ContactData& contact, EnergyDerivatives* f);

#endif // CONTACT_H_

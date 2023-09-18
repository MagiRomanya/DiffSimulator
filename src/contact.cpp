#include "contact.hpp"

void compute_sphere_contact_geometry(const Sphere& sphere, const Vec3& point, ContactData* out) {
    const Vec3 delta = point - sphere.center;
    const Scalar distance = delta.norm();
    out->normal = delta / distance;
    out->s_distance = distance - sphere.radius;
}

void compute_plane_contact_geometry(const Plane& plane, const Vec3& point, ContactData* out) {
    out->s_distance = plane.normal.dot(point - plane.center);
    out->normal = plane.normal;
}

Vec3 compute_contact_force(const ContactData& contact) {
    return - contact.stiffness * contact.s_distance * contact.normal;
}

Mat3 compute_contact_force_position_jacobian(const ContactData& contact) {
    return - contact.stiffness * contact.normal * contact.normal.transpose();
}

void compute_contact_energy_derivatives(const ContactData& contact, EnergyDerivatives* f) {
    // Magnitude calculation
    Vec3 force = compute_contact_force(contact);
    Mat3 df_dx = compute_contact_force_position_jacobian(contact);

    const unsigned int index = contact.index;
    // Set force
    f->force[index] = force[0];
    f->force[index+1] = force[1];
    f->force[index+2] = force[2];

    // Set force position jacobian
    for (size_t i = 0; i < 3; i++){
        for (size_t j = 0; j < 3; j++) {
            f->df_dx_triplets.push_back(Triplet(index + i, index + j, df_dx(i, j)));
        }
    }
}

void ContactManager::find_contacts(const PhysicsState &state, std::vector<ContactData> &contacts) {
    // Collision finding
    ContactData contact;
    // NOTE Here we assume that all the state is made out of particles,
    // in the future maybe this will not be true.
    for (size_t i = 0; i < state.q.size(); i+=3) {
        unsigned int index=i;
        contact.index = index;
        Vec3 point = get_particle_position(&state, index);

        // Spheres
        for (size_t s = 0; s < sphere_colliders.size(); s++) {
            compute_sphere_contact_geometry(sphere_colliders[s], point, &contact);
            if ( contact.s_distance < 0 ) contacts.push_back(contact);
        }
        // Planes
        for (size_t p = 0; p < plane_colliders.size(); p++) {
            compute_plane_contact_geometry(plane_colliders[p], point, &contact);
            if ( contact.s_distance < 0 ) contacts.push_back(contact);
        }
    }
}

void ContactManager::compute_contacts_energy_derivatives(const std::vector<ContactData>& contacts, EnergyDerivatives* f) {
    for (size_t i = 0; i < contacts.size(); i++) {
        compute_contact_energy_derivatives(contacts[i], f);
    }
}

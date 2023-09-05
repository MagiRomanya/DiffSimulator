#include "spring.hpp"
#include "linear_algebra.hpp"
#include "two_body_helper.hpp"
#include <iostream>
#include <vector>

void Spring::calculate_energy_derivatives(const PhysicsState* state, EnergyDerivatives* f) {
    get_state(state);
    two_body::add_force(p1, p2, force(), f);
    two_body::add_force_derivatives(p1, p2, force_position_derivative(), Mat3::Zero(), f);
    std::vector<unsigned int> parameters = {k.index, L0.index};
    two_body::add_force_parameter_derivatives(p1, p2, force_parameters_derivative(), f, parameters);
}

void Spring::get_state(const PhysicsState* state) {
    x1 = get_particle_position(state, p1);
    x2 = get_particle_position(state, p2);

    L = (x1 - x2).norm();
}

Vec3 Spring::force() const {
    /* Computes the spring force */
    Vec3 f = -k.value * (L - L0.value) * (x1 - x2) / L;
    return f;
}

Scalar Spring::energy() const {
    /* Computes the spring's energy */
    Scalar energy = 0.5 * k.value * (L - L0.value) * (L - L0.value);
    return energy;
}

Mat3 Spring::force_position_derivative() const {
    // u is the normalized vector between particles 1 and 2
    Vec3 u = (x1 - x2) / L;
    // Initialize the derivative matrix
    Mat3 df_dx = (L - L0.value) * Mat3::Identity();
    // The u · u.transpose() matrix
    Mat3 uut = u * u.transpose();

    // Calculate the final derivative matrix
    df_dx = - k.value / L * (df_dx + L0.value * uut);

    return df_dx; // 3x3 matrix
}

Mat3 Spring::force_velocity_derivative() const {
    return Mat3::Zero();
}

Mat Spring::force_parameters_derivative() const {
    const unsigned int nDoF        = 3; // 2 body interaction
    const unsigned int nParameters = 2; // stiffness + rest length
    Mat df_dp = Mat(nDoF, nParameters);
    Vec3 u = (x1 - x2)/L;
    Vec3 df_dk = - (L - L0.value) * u;
    Vec3 df_dL0 = k.value * u;
    df_dp << df_dk, df_dL0;
    return df_dp;
}

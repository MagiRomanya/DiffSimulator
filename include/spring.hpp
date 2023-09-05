#ifndef SPRING_H
#define SPRING_H

#include <math.h>
#include <ostream>
#include <vector>

#include "physics_state.hpp"

class Spring {
    public:
        Spring(unsigned int pi, unsigned int pj, Parameter k, Parameter L0)
            : p1(pi), p2(pj), k(k), L0(L0) {}

        unsigned int p1, p2;

        // The spring stiffness & the spring rest length
        Parameter k, L0;

        // Positions and precomputed distance
        Vec3 x1, x2;
        Scalar L;

        void calculate_energy_derivatives(const PhysicsState* state, EnergyDerivatives* f);

    private:
        void get_state(const PhysicsState* state);

        Scalar energy() const;

        Vec3 force() const;

        Mat3 force_position_derivative() const;

        Mat3 force_velocity_derivative() const;

        Mat force_parameters_derivative() const;
};

#endif

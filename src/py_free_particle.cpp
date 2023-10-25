#include <raylib.h>

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "simulable.hpp"
#include "simulable_generators.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"

class FreeParticle {
    void set_simulation(Simulation sim) {
        simulation = sim;
        state = simulation.getInitialState();
    }

    SparseMatrix getInitialPositionJacobian() {
        SparseMatrix dx0dp(getDoF(), getDoF());
        dx0dp.setFromTriplets(
            simulation.simulation_parameters.dq0_dp_triplets.begin(),
            simulation.simulation_parameters.dq0_dp_triplets.end());
        return dx0dp;
    }

    SparseMatrix getInitialVelocityJacobian() {
        SparseMatrix dv0dp(getDoF(), getDoF());
        dv0dp.setFromTriplets(
            simulation.simulation_parameters.dq_dot0_dp_triplets.begin(),
            simulation.simulation_parameters.dq_dot0_dp_triplets.end());
        return dv0dp;
    }

    SparseMatrix getMassMatrix() {
        return mass_matrix;
    }

    Scalar getTimeStep() {return simulation.simulation_parameters.TimeStep;}

    unsigned int getDoF() { return simulation.simulation_parameters.q0.size();}

    unsigned int getDiffParameters() { return simulation.simulation_parameters.p.size();}

    Vector getPosition() { return state.q;}

    Vector getVelocity() { return state.q_dot;}

    void fill_containers() {
        mass_matrix.setFromTriplets(simulation.simulation_parameters.mass.begin(),
                                    simulation.simulation_parameters.mass.end());

        f = EnergyDerivatives(getDoF(), getDiffParameters());
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

        const Scalar h = getTimeStep();
        equation_matrix = mass_matrix - h * df_dv - h * h * df_dx;
        handle_frozen_dof(simulation.simulation_parameters.frozen_dof, &equation_matrix);
    }

    void set_state(Vector xi, Vector vi) {
        state.q = xi;
        state.q_dot = vi;
    }

    SparseMatrix getEquationMatrix() { return equation_matrix; }
    Vector getForce() { return f.force; }
    Vector getDiffParameteres() { return simulation.simulation_parameters.p; }
    Mat getParameterJacobian() { return f.df_dp; }
    SparseMatrix getForcePositionJacobian() { return df_dx; }
    SparseMatrix getForceVelocityJacobian() { return df_dv; }

    Simulation get_simulation() {return simulation;}

    private:
        Simulation simulation;
        PhysicsState state;
        EnergyDerivatives f;

        SparseMatrix df_dx, df_dv, equation_matrix, mass_matrix;

        Camera3D* camera = nullptr;
};

class PyRenderer {
    PyRenderer(){
        const int screenWidth = 800*2;
        const int screenHeight = 450*2;
        // Disable raylib info and warnings logs
        SetTraceLogLevel(LOG_ERROR);
        InitWindow(screenWidth, screenHeight, "Simulator");
        // Create a camera with 300 fps
        camera = create_camera(200);
    }

    private:
        Camera3D camera;
};

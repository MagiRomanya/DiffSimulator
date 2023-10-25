#ifndef SIMULATION_H_
#define SIMULATION_H_


#include "interaction_manager.hpp"
#include "physics_state.hpp"
#include "simulable.hpp"
#include "simulation_parameters.hpp"
#include "contact.hpp"

struct Simulation {
    SimulationParameters simulation_parameters;

    InteractionManager interaction_manager;

    ContactManager contact_manager;

    std::vector<Simulable> simulables;

    PhysicsState getInitialState() {
        PhysicsState state;
        state.q = simulation_parameters.q0;
        state.q_dot = simulation_parameters.q_dot0;
        state.time = 0;
        return state;
    }

    void step(PhysicsState* state) {
        const unsigned int nDoF = simulation_parameters.q0.size();
        const unsigned int nParameters = simulation_parameters.p.size();
        EnergyDerivatives f(nDoF, nParameters);

        // Interaction force and derivatives callculation
        interaction_manager.calculate_energy_derivatives(state, &f);

        // Collision search and response
        std::vector<ContactData> contacts;
        contact_manager.find_contacts(*state, contacts);
        contact_manager.compute_contacts_energy_derivatives(contacts, &f);

        integrate(simulation_parameters, state, f);
    }
};

#endif // SIMULATION_H_

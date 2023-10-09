#ifndef PYSIMULATION_H_
#define PYSIMULATION_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "raylib.h"
#include "simulable.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"
#include <stdbool.h>
#include <vector>

#define GRID_NODE_SIDE 20
#define GRID_WIDTH_LENGTH 5.0f
#define ENABLE_CONTACT

class PySimulation {
    public:
        PySimulation();
        ~PySimulation();

        /* Two differentiable parameters. */
        PySimulation(Scalar k, Scalar k_bend, bool graphics=false);

        /* Three differentiable parameters. */
        PySimulation(Scalar k, Scalar k_bend, Scalar tilt_angle, bool graphics=false);

        /* A lot of differentiable parameters. */
        PySimulation(std::vector<Scalar> k, std::vector<Scalar> k_bend, bool graphics=false);

        void fill_containers();

        void set_state(Vector xi, Vector vi);

        SparseMatrix getEquationMatrix();

        Vector getEquationVector();

        Vector getForce();

        Vector getPosition();

        Vector getVelocity();

        Vector getDiffParameteres();

        SparseMatrix getMassMatrix();

        Mat getParameterJacobian();

        SparseMatrix getForcePositionJacobian();

        SparseMatrix getForceVelocityJacobian();

        int getDoF();

        Scalar getTimeStep();

        std::array<unsigned int, 2> getGridDimensions();

        std::array<unsigned int, 2> getNumberOfSprings();

        void render_state();

        SparseMatrix getInitialPositionJacobian();

        SparseMatrix getInitialVelocityJacobian();

        void reset_simulation(Scalar stiffness, Scalar bend_stiffness);

        void reset_simulation(Scalar stiffness, Scalar bend_stiffness, Scalar tilt_angle);

        void reset_simulation(std::vector<Scalar> stiffness, std::vector<Scalar> bend_stiffness);

    private:
        void set_up_simulation();

        void resize_containers();

        void create_mesh();

        const Scalar node_mass = 1.0;
        PhysicsState state;
        EnergyDerivatives f = EnergyDerivatives(0, 0);
        Simulation simulation;

        unsigned int grid_n, grid_m;
        unsigned int n_flex, n_bend;

        SparseMatrix df_dx, df_dv, mass_matrix, equation_matrix;
        Vector equation_vector;
        Simulable mass_spring;

        Scalar angle = 0;

        // Rendering stuff
        bool graphics = false;
        Camera3D camera;
        bool game_paused = false;
        Mesh cloth_mesh;
        std::vector<Scalar> vertices;
        std::vector<unsigned int> indices;
        Material cloth_material;
};

#endif // PYSIMULATION_H_

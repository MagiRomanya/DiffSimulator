#ifndef SIMULABLE_H_
#define SIMULABLE_H_

/*
 * A simulable object defines the degrees of freedom and the number of phyisical parameters it needs.
*/
struct Simulable {
        unsigned int index;           // Start index in the generalized coordinates in phyisics state
        unsigned int nDoF;            // Degrees of freedom that the simulable needs
        unsigned int parameter_index; // Start index in the global parameters vector in SimulationParamters
        unsigned int nParameters;     // Number of parameters that the simulable needs
};

#endif // SIMULABLE_H_

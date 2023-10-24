#include <pybind11/detail/common.h>
#include <pybind11/eigen.h> // for sparse and dense matrices
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // for std::vector as argument
#include <vector>

#include "linear_algebra.hpp"
#include "pysimulation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(symulathon, m) {
    m.doc() = "A simple simulation backend.";

    // Main simulation interface
    py::class_<PySimulation>(m, "Simulation")
        .def(py::init()) // Constructor
        .def(py::init<Scalar, Scalar, bool>(),
             py::arg("k"),
             py::arg("k_bend"),
             py::arg("graphics")=false)
        .def(py::init<Scalar, Scalar, Scalar, bool>(),
             py::arg("k"),
             py::arg("k_bend"),
             py::arg("tilt_angle"),
             py::arg("graphics")=false)
        .def(py::init<std::vector<Scalar>, std::vector<Scalar>, bool>(),
             py::arg("k"),
             py::arg("k_bend"),
             py::arg("graphics")=false)
        .def("fill_containers", &PySimulation::fill_containers)
        .def("set_state", &PySimulation::set_state)
        .def("getEquationMatrix", &PySimulation::getEquationMatrix)
        .def("getEquationVector", &PySimulation::getEquationVector)
        .def("getForce", &PySimulation::getForce)
        .def("getPosition", &PySimulation::getPosition)
        .def("getVelocity", &PySimulation::getVelocity)
        .def("getDiffParameters", &PySimulation::getDiffParameteres)
        .def("getMassMatrix", &PySimulation::getMassMatrix)
        .def("getParameterJacobian", &PySimulation::getParameterJacobian)
        .def("getForcePositionJacobian", &PySimulation::getForcePositionJacobian)
        .def("getForceVelocityJacobian", &PySimulation::getForceVelocityJacobian)
        .def("getDoF", &PySimulation::getDoF)
        .def("getTimeStep", &PySimulation::getTimeStep)
        .def("render_state", &PySimulation::render_state)
        // .def("getSpringIndices", &PySimulation::getSpringNodeIndices)
        // .def("getBendSpringIndices", &PySimulation::getBendSpringNodeIndices)
        .def("getGridDimensions", &PySimulation::getGridDimensions)
        .def("getSpringNumbers", &PySimulation::getNumberOfSprings)
        .def("getInitialPositionJacobian", &PySimulation::getInitialPositionJacobian)
        .def("getInitialVelocityJacobian", &PySimulation::getInitialVelocityJacobian)
            // Overloaded methods
        .def("reset_simulation", static_cast<void (PySimulation::*)(Scalar, Scalar)>(&PySimulation::reset_simulation))
        .def("reset_simulation", static_cast<void (PySimulation::*)(Scalar, Scalar, Scalar)>(&PySimulation::reset_simulation))
        .def("reset_simulation", static_cast<void (PySimulation::*)(std::vector<Scalar>, std::vector<Scalar>)>(&PySimulation::reset_simulation))
        .def("reset_simulation", static_cast<void (PySimulation::*)(std::vector<Scalar>)>(&PySimulation::reset_simulation))
        .def("reset_simulation", static_cast<void (PySimulation::*)(Scalar, Scalar, std::vector<Scalar>)>(&PySimulation::reset_simulation))
        .def("window_should_close", &PySimulation::window_should_close)
        .def("is_key_pressed", &PySimulation::is_key_pressed)
        ;
}

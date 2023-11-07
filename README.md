# DiffSimulator

This project aims to develop a particle system differentiable simulation.
It is designed to work as a standalone simulator binary as well as an interface that can be used through python.

## Quick start
The project is supposed to work in both Linux and Windows operating systems.
Specifically, we have tested the build in both the gcc and visual studio compilers.
This is a cmake project, and running the following commands should be enough to build the whole project.
``` sh
mkdir build
cd build
cmake ..
cmake -build .
```
The project relies on third party software that should be automatically be added as a git sub-module the first time the cmake command is run.

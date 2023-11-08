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

## Demonstration of optimizations using the differentiable simulation.

- Here we are optimizing each initial velocity of the cloth. Our goal is that the cloth reaches the sphere collider in a certain configuration.

https://github.com/MagiRomanya/DiffSimulator/assets/68277239/255870c9-ddb6-4a40-a97a-89e1596d6e5f

- Here we are optimizing 2166 elasticity paramters of the cloth. We want the final state of the cloth in frame 100 to match as best as possible a certain recorded frame.

https://github.com/MagiRomanya/DiffSimulator/assets/68277239/a4c7773d-866b-42e8-89f9-bbd1af12eb5d


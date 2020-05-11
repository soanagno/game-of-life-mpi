# acse-6-individual-assignment-acse-sa1619
Sokratis Anagnostopoulos (CID: 00986040)


# Parallel implementation of Conway's Game of Life using MPI on C++.


__Requirements:__

- For .cpp --> MPI C++ library settings.
- For .m --> Matlab R2019b


__Instructions:__

Change the following variables in the main: iter, i_grid, j_grid.
When calling the main function of the simulation, "simulate(i_grid, j_grid, iter, 1, true)", 
the last two input variables are:
the initial condition (0 random, 1 glider) and the periodic boundary condition switch (true, false), respectively.


__Description of attached video files:__

- Video A: 100 x 100 grid divided into 4 cores demonstrating periodic boundary condition and 2x2 special case with gliders.

- Video B: 100 x 100 grid divided into 13 cores demonstrating periodic boundary condition and 1x13 special case with gliders.

- Video C: 100 x 100 grid divided into 16 cores demonstrating periodic boundary condition with gliders traveling towards opposite direction.

- Video D: 100 x 100 grid divided into 16 cores demonstrating NON-periodic boundary condition with gliders.

- Video E: 100 x 100 grid divided into 18 cores demonstrating vertical configuration with gliders.

- Video F: 23 x 41 grid divided into 9 cores demonstrating horizontal configuration with random inital condition.

- Video G: 1000 x 100 grid divided into 9 cores with periodic boundary condition.

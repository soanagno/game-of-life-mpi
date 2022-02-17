// ----- GAME OF LIFE PARALLEL IMPLEMENTATION WITH MPI ----- //

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>

using namespace std;

// Definitions
int id, p, time_tot;
int i_grid, j_grid, m, n, i, j, b, i_inner, j_inner, k, iter;
int rows, columns, id_row, id_column;

// Vectors for the Datatype definitions
vector<vector<int>> block_length(16);
vector<vector<MPI_Aint>> addresses(16);
vector<vector<MPI_Datatype>> typelist(16);
vector<MPI_Datatype> comm_type(16);
vector<MPI_Aint> temp(16);

// Maping vector for Datatype order, explained in the text
vector<int> type_map = { 0, 8, 1, 11, 0, 9, 3, 10, 2 };

void closest_factors(int input, int& rows, int& columns) {

	/*  
	An efficient way to calculate closest factors of a number
	for the dimensions (rows, columns) of the processors' grid.

	---Inputs---
	input: integer number of cores
	rows: integer factor 1
	columns: integer factor 2	
	*/

	int sqrt_input = (int)sqrt(input);
	while (input % sqrt_input != 0)
		sqrt_input--;
	rows = sqrt_input; columns = input / sqrt_input;
	if (id == 0) cout << "Processor grid is divided: " << p << " into: " << rows << " by " << columns << endl;
}

void id_to_index(int id, int& id_row, int& id_column) {

	/*
	Convert core id to index.

	---Inputs---
	id: processor id 
	id_row: processor x coordinate
	id_column: processor y cooridnate
	*/

	id_column = id % columns;
	id_row = id / columns;
}

int id_from_index(int id_row, int id_column, bool periodic) {

	/*
	Get core id from index depending on the boundary condition.

	---Inputs---
	id_row: processor x coordinate
	id_column: processor y cooridnate
	periodic: deploy periodic condition (true) non-periodic condition (false)

	---Outputs---
	id of communication processor depending on the boundary condition
	*/

	if (periodic == false) {

		// Return the index neglecting indexes out of the boundaries
		// to apply non-periodic condition
		if (id_row >= rows || id_row < 0)
			return -1;
		if (id_column >= columns || id_column < 0)
			return -1;
		return id_row * columns + id_column;
	}
	else if (periodic == true) {

		// Return appropriate id efficiently without the need of extra
		// imposed conditions or additional if statements inside the communications loop. 
		// The following formula will return the neighbor ids for all internal
		// cores and will return the symmetric neighbor for a core next to the boundary.
		return (id_row + rows) % rows * columns + (id_column + columns) % columns;
	}
}

void setup_types(bool* grid, vector<vector<int>> corners) {

	/*
	Function that sets-up all the datatypes for each individual processor using only two for loops.
	Since different amount of elements needs to be considered for each Datatype, depending on if 
	a corner or a side is being sent or received, the Data types are grouped in two for loops,
	one for the corner and one for the side types. This a more general approach that can facilitate
	the application of this algorithms for N datatypes, as explained in the text.

	---Inputs---
	grid: 1D pointer of processor internal grid of booleans
	corners: 2D vector defining processors corner coordiantes (x, y)
	*/

	// Create corner types
	for (i = 0; i < 8; i++) {
		block_length[i].push_back(1);
		MPI_Get_address(&grid[corners[i][0] * n + corners[i][1]], &temp[i]);
		addresses[i].push_back(temp[i]);
		typelist[i].push_back(MPI_C_BOOL);

		MPI_Type_create_struct(block_length[i].size(), &block_length[i][0], &addresses[i][0], &typelist[i][0], &comm_type[i]);
		MPI_Type_commit(&comm_type[i]);  // Commit Datatype
	}

	// Create side types
	int jj, pixels, flag = 1; // "flag" variable changes sign each iteration. 1 is for horizontal sides and -1 is for vertical.
	for (i = 8; i < 16; i++) {
		jj = i - 8;

		// Variable that defines the size of the data to be compressed in the Datatype.
		// Switches between the size of horizontal and vertical lines using the flag variable.
		pixels = j_inner * (flag + 1)/2 - i_inner * (flag - 1)/2;  
		for (j = 1; j <= pixels; j++) {
			block_length[i].push_back(1);
			if (flag == 1) MPI_Get_address(&grid[corners[jj][0] * n + j], &temp[i]);
			else MPI_Get_address(&grid[j * n + corners[jj][1]], &temp[i]);
			addresses[i].push_back(temp[i]);
			typelist[i].push_back(MPI_C_BOOL);
		}
		flag = -flag;  // "flag" changes sign here

		MPI_Type_create_struct(block_length[i].size(), &block_length[i][0], &addresses[i][0], &typelist[i][0], &comm_type[i]);
		MPI_Type_commit(&comm_type[i]);  // Commit Datatype
	}

}

void Fill_Random(bool* grid, int i_inner, int j_inner) {

	/*
	Returns a 1D pointer array of random binary values

	---Inputs---
	grid: 1D pointer of processor internal grid of booleans
	i_inner: inner core row
	j_inner: inner core columns
	*/

	double rand_bool;
	for (i = 1; i <= i_inner + 1; i++)
		for (j = 1; j <= j_inner + 1; j++)
			grid[i * j_inner + j] = rand() % 2;
}

void Fill_Glider(bool* grid) {

	/* 
	Creates a glider for testing purposes

	---Inputs---
	grid: 1D pointer of processor internal grid of booleans
	*/
	
	int gx = 20;
	grid[(3 + gx) * n + 3] = 1;
	grid[(3 + gx) * n + 4] = 1;
	grid[(3 + gx) * n + 5] = 1;
	grid[(4 + gx) * n + 5] = 1;
	grid[(5 + gx) * n + 4] = 1;

}

void simulate(int i_grid, int j_grid, int iter, bool condition, bool boundary) {

	/*
	Sets up core grid configuration and simulates Game of Life.
	Each core produces a separate text file which is later used in post-processing.

	---Inputs---
	i_grid: rows of full computational grid
	j_grid: columns of full computational grid
	iter: number of simulation generations
	condition: initial condition (0 = random, 1 = glider)
	boundary: boundary condition (true = periodic, false = non-periodic/wall)
	*/

	closest_factors(p, rows, columns);  // Calculate optimal processor grid dimensions
	id_to_index(id, id_row, id_column); // Return index from core id


	// Efficient way to divide full grid to each core.
	// In the case where the full grid dimension cannot be divided into the requested number of cores,
	// the first few cores take 1 more line and/or column than the few last ones.
	// Make a first approximation of the core's inner allocated cells using the / (div) operator:
	i_inner = i_grid / rows; j_inner = j_grid / columns;

	// Calculate the remainder of the division for i and j
	int modi = i_grid % rows; int modj = j_grid % columns;

	// Keep initial approximated dimensions for i and j
	int i_inner0 = i_inner; int j_inner0 = j_inner;

	// Allocate 1 extra row and/or column to the first few cores until all the remainder has been allocated
	// to the first few cores. The rest of the cores keep the initial approximation ( using the div operator).
	// Note that this method avoids using any for loops.
	if (id_row + 1 <= modi) i_inner += 1;
	if (id_column + 1 <= modj) j_inner += 1;

	if (id == 0) cout << "Processor " << id << " takes a " << i_inner << " by " << j_inner << " grid" << endl;

	m = i_inner + 2;  // Rows including ghost nodes
	n = j_inner + 2;  // Columns including ghost nodes

	fstream myfile;  // File to be used for animations
	myfile.open(to_string(id) + "_pixels.txt", fstream::out);
	// First line of the file states no. of cores used, rows and columns of core grid, inner dimensions of each core and no. of iterations.
	myfile << p << " " << rows << " " << columns << " " << i_inner << " " << j_inner << " " << iter << endl;

	// Vector that calculates the corners of each INDIVIDUAL core (clockwise) which will later be used in the animations (to define each core area),
	// but mainly in the definition of the data types. Note that a simple routine could replace this vector definition in future verions
	// so that the corner points of each processor are defined even for more general cases where more than 1 core has a common interface
	// with a neighbor core. In that case the need for more Datatypes would be necessary, which could easily be implemented with this code's approach.
	vector<vector<int>> corners = { {1, 1},  {1, n - 2}, {m - 2, n - 2},  {m - 2, 1}, {0, 0},  {0, n - 1}, {m - 1, n - 1},  {m - 1, 0} };

	// Calculate corner points for grid plot (only the right bottom ones of each core are needed to define the grid of the animations)
	int l = (id_row + 1) * i_inner0 + modi;
	if (id_row + 1 <= modi && modi != 0) l = (id_row + 1) * i_inner0 + (id_row + 1);

	int r = (id_column + 1) * j_inner0 + modj;
	if (id_column + 1 <= modj && modj != 0) r = (id_column + 1) * j_inner0 + (id_column + 1);

	myfile << l << " " << r << endl;  // Second line of the file contains the l, r coordinates (x, y) of the core's corners.

	// Initialisation of array a
	bool* a = new bool[m * n]{ 0 };
	bool* a_new = new bool[m * n]{ 0 };

	// Initial Conditions
	if (condition == 0) Fill_Random(a, i_inner, j_inner);
	else if (condition == 1) Fill_Glider(a);

	// Copy a to a_new and print out first image
	for (i = 1; i <= m - 2; i++) {
		for (j = 1; j <= n - 2; j++) {
			a_new[i * n + j] = a[i * n + j];
			myfile << a[i * n + j] << " ";
		}
		myfile << endl;
	}

	// Tag number vectors for communication sequences
	vector<int> tag_num_send;
	vector<int> tag_num_receiv;
	if (rows == 1) {
		// Covers the case where the number of used cores is prime.
		// This case requires a special tag number vector since it's
		// the only case where a processor communicates with itself
		// (through the top and bottom boundaries).
		tag_num_send = { 0, 4, 1, 5, 0, 5, 3, 6, 2 };
		tag_num_receiv = { 2, 6, 3, 5, 0, 5, 1, 4, 0 };
	}
	else {
		// Covers all other cases, along with all 2 X N configurations.
		// The need for the use of a tag number vector mainly arises
		// from the 2 x N cases, where some cores communicate with their neighbors
		// more than once, using different Datatypes. This is also explained in the text.
		tag_num_send = { 0, 4, 1, 7, 0, 5, 3, 6, 2 };
		tag_num_receiv = { 2, 6, 3, 5, 0, 7, 1, 4, 0 };
	}

	setup_types(a, corners);  // Sets up required Datatypes for communication of current core with its neighbors

	MPI_Request* request = new MPI_Request[16];  // Creates vector of communication requests

	// Beginning of generations' loop
	for (k = 0; k < iter; k++) {

		// "cnt" counts the number of (sends + receives) occuring 
		// while "map" counts the total times of the accessed neighbor 
		// cores (including the current/central core), which in this case goes up to 9.

		int cnt = 0; int map = 0;  
		for (i = -1; i <= 1; i++) {
			for (j = -1; j <= 1; j++) {

				int com_i = id_row + i;
				int com_j = id_column + j;
				int com_id = id_from_index(com_i, com_j, boundary);

				// Efficiently does all required communications at once, by taking advantage of the Datatype mapping vector "type_map".
				// The Datatypes are accessed in the correct order (first corners and then sides, clockwise) using the type_map vectors.
				if (com_id >= 0 && com_id < p && map != 4) {
					MPI_Isend(MPI_BOTTOM, 1, comm_type[type_map[map]], com_id, tag_num_send[map], MPI_COMM_WORLD, &request[cnt * 2]);
					MPI_Irecv(MPI_BOTTOM, 1, comm_type[type_map[map] + 4], com_id, tag_num_receiv[map], MPI_COMM_WORLD, &request[cnt * 2 + 1]);
					cnt++;
				}
				map++;
			}
		}

		// Waits until all communication requests are completed to continue computations
		MPI_Waitall(2 * cnt, request, MPI_STATUSES_IGNORE) == MPI_SUCCESS;  


		//--- GAME OF LIFE ---//
		for (i = 1; i <= m - 2; i++) {
			for (j = 1; j <= n - 2; j++) {

				// Calculate number of neighbours, b
				b = a[i * n + j - 1] + a[i * n + j + 1] + \
					a[(i - 1) * n + j - 1] + a[(i - 1) * n + j] + a[(i - 1) * n + j + 1] + \
					a[(i + 1) * n + j - 1] + a[(i + 1) * n + j] + a[(i + 1) * n + j + 1];

				// Apply fertility rules using a minimal number of if-statements
				if (b != 2 and b != 3) a_new[i * n + j] = 0;
				else if (b == 3) a_new[i * n + j] = 1;

			}
		}

		// Replace a with a_new and print to file
		for (i = 1; i <= m - 2; i++) {
			for (j = 1; j <= n - 2; j++) {
				a[i * n + j] = a_new[i * n + j];
				myfile << a[i * n + j] << " ";
			}
			myfile << endl;
		}
	}

	myfile.close();

	// Delete pointers to avoid memory leakage
	delete[] request;
	delete[] a_new;
	delete[] a;
	
}


int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 1000);

	// Initialisations

	// fstream myfile;
	clock_t time_req;
	// int g, maxgrid = 11;
	// g = maxgrid;  // Max grid dimension that can be varied for parametric runs
	
	int max_it = 1; // Max number of simulation run for timing averaging
	
	// Definition of full grid dimensions and generations
	iter = 200;  // Max number of generations
	i_grid = 100;  // Max rows
	j_grid = 100;  // Max columns


	// for (g = 20; g <= maxgrid; g += 100) {  // used for parametric HPC runs
		for (int psi = 0; psi < max_it; psi++) {

			time_req = clock();  // Start timing

		
			// Initiate main simulation, 0 for random / 1 for glider, true/false for boundary conditions

			// simulate(g, g, iter, 0, true);  // Initiate main simulation
			simulate(i_grid, j_grid, iter, 1, true);  

			time_req = clock() - time_req;  // End timing
			time_tot += time_req;  // Total cumulative time

		}

		time_tot = time_tot / max_it; // Average timing

		/* 
		Print cpu timings to file (for HPC runs)
		if (id == 0) {
			MPI_Barrier(MPI_COMM_WORD);  // Barrier to make sure all CPUs have reached this point before time measurement
			myfile.open(to_string(p) + "_times.txt", ios_base::app);
			myfile << g << " " << p << " " << (double)time_req / CLOCKS_PER_SEC << endl;
			myfile.close();
		}
		*/
	// }
	
	MPI_Finalize();
	
}

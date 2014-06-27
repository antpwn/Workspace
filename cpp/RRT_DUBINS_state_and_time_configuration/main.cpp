#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
using namespace std;

#include <GL/glut.h>
#include <math.h>
#include <kdtree.h>

#include <ctype.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
/////////////////////////////////  begin definitions for RRT algorithm ///////////////////////////////////////////

//default_random_engine generator;
//uniform_real_distribution<double> distribution(0.0, 1.0);
float** vertices_final;
int** edges_final;
float** traj_final;
float** solutions_states_relaxed_Dubins_final;
float** solutions_marks_relaxed_Dubins_final;
float** solutions_states_full_Dubins_final;
float** solutions_marks_full_Dubins_final;
float** solutions_TTH_Relaxed_Dubins_matrix_final;
float** solutions_TTH_full_Dubins_matrix_final;
vector<int> sizes_values_final;
vector<float> monte_carlo_avg_Relaxed_Dubins_final;
vector<float> monte_carlo_avg_full_Dubins_final;
vector<float> convergence_values_final;

float** color_convergence_matrix;

int hit_index;
int hit_index_RD;

const float PI = 3.141592653589793238463;

#define d2r  PI/180
#define r2d  180/PI

#define states_number 4
#define dubins_mark_size 12
#define plain_col  2
#define obst_col  3
//#define obst_col  2

#define edges_row  2   
#define traj_row  13  // traj containes a conbination of cost function value for the corisponding edge in 0 and path marks in 1-12
//					cost value																	--- 0
// path_marks		R, c, L,																	--- 1-3 
//					initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 4-7
//					v_norm,																		--- 8
//					finalCircle_centre[0],14finalCircle_centre[1],finalArcAngle,phi_contact,		--- 9-12

const int run_mode = 1; //	1 - RRT   2- RRT*
const int steer_length_mode = 1; //	1 - length to sample point   2- constant ds length   3- random length (0,ds]
int plot_mode; //  1 - only tree  2- only best solution  3 - full tree and solution
int Dubins_sol_mode; //  1 - Relaxed Dubins  2- Full Dubins
const int goal_mode = 1; //	1- point goal, 2- area goal, 3- moving target goal
const int solution_mode = 2; //	1- single solution for the simulation, 2- multiple solutions for the simulation (monte carlo)
const int kdtree_mode = 2; //	1- old search for nearest and near, 2- new log time kd-tree search for nearest and near

const int steps = 1000000;
const int monte_carlo_number = 50;
int convergence_size = 21;

const int angle_state_vector_size = 1;
const int angle_state_vector[angle_state_vector_size] = { 3 };

int edges_col = steps;
int traj_col = steps;

const float d_radius = 1;

const float ds = 1.0;
//const float ds = 5;
const float ds_obs = 0.05; // size of segments to check for obstacles collision
const float ds_plot = 0.1; // size of segments to plot
const float V_max = 10; // in this simulation we will always travel at V_max
const float Sim_R_max = 10;
const float Sim_R_min = 1;

const int dim = states_number;
const float gama = 4;

const float max_dis_av = 30000; // maximum possible distance between any two vertices

const float hit_dis = 0.5;
const float hit_dis_RD = 0.5; // value for relaxed Dubins

const float plain[states_number][plain_col] = { { 0, 20 }, //x1 - defines edges of the potential search space
		{ 0, 20 }, //x2 - plain row i is made of two numbers that state the min and max
		{ -PI, PI }, //x3 = heading angle, has the circling feature
		{ 0, 20 } }; //x4 = time - allowed values of the i state

const float space_volume = 20 * 20 * 2 * PI * 20;

// define obstacles
const int obst_num = 0;
const float obstacles[1 * (states_number - 1)][obst_col] = { { 9, 11, -2 }, // same logic as plain + velocity projection to that axis
		{ 0, 18, 3 }, //, each 'states_number' rows are an obstacle.
		{ -PI, PI, 0 } };
//const float obstacles[1 * (states_number)][obst_col] = { 	{ 9, 11},  // same logic as plain + velocity projection to that axis
//																{ 0, 10},	//, each 'states_number' rows are an obstacle.
//																{-PI,PI},
//																{0,20	}};

const float x1_0 = 1;
const float x2_0 = 1;
const float x3_0 = 45 * d2r;
const float x4_0 = 0;

const float goal_point_x1 = 15;
const float goal_point_x2 = 15;
const float goal_point_x3 = 135 * d2r;
const float goal_point_x4 = 40;

float euclidean_dis = pow(
		pow(goal_point_x2 - x2_0, 2) + pow(goal_point_x1 - x1_0, 2), 0.5);

int firts_point = 0;

/////////////////////////////////  begin definitions for OpenGL drawing ///////////////////////////////////////////

#define window_pos_left  100
#define window_pos_hight  50

#define window_width  1100
#define window_height 700
float red = 1.0f, blue = 0.0f, green = 1.0f;

// angle of rotation for the camera direction
float angle = 0.0;
// actual vector representing the camera's direction
float lx = 0.0f, lz = -1.0f;
// XZ position of the camera
float x = 0.0f, z = 20.0f;

// center position
float x_center = plain[0][0] + (plain[0][1] - plain[0][0]) / 2, y_center =
		plain[1][0] + (plain[1][1] - plain[1][0]) / 2, z_center = plain[2][0]
		+ (plain[2][1] - plain[2][0]) / 2;

// camera position
float x_cam = x_center + 20.0, y_cam = y_center + 45.0, z_cam = z_center + 40.0;

// up position
float x_up = 0.0f, y_up = 0.0f, z_up = 1.0f;

float deltaAngle = 0.0f;
float deltaMOveX_usr = 0.0f;
float deltaMOveY_usr = 0.0f;
float deltaMOveZ_usr = 0.0f;
float deltaMove = 0;
int xOrigin_L = -1;
int yOrigin_L = -1;
int xOrigin_R = -1;
int yOrigin_R = -1;
int w_up = 0;
int w_down = 0;
float delta_center_cam = 0;
float delta_center_camXY = 0;
float delta_x = 0;
float delta_y = 0;
float delta_z = 0;

vector<float> comm_mov;

float sa = 0;
float ca = 0;
float sb = 0;
float cb = 0;

const float colormap[3][64] = { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625,
		0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9375,
		0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625, 0.5 }, { 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375,
		0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		0.9375, 0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625, 0.5, 0.4375, 0.375,
		0.3125, 0.25, 0.1875, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0 }, { 0.5625, 0.6250, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 0.9375, 0.875, 0.8125, 0.75, 0.6875, 0.625, 0.5625, 0.5,
		0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125, 0.0625, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };

//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// RRT functions //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

////////
vector<float> vector_turn_x(vector<float> vec_in, float theta) {
	vector<float> vec_out;
	vec_out.push_back(vec_in[0]);
	vec_out.push_back(vec_in[1] * cos(theta) - vec_in[2] * sin(theta));
	vec_out.push_back(vec_in[1] * sin(theta) + vec_in[2] * cos(theta));
	return vec_out;
}
////////
vector<float> vector_turn_y(vector<float> vec_in, float phi) {
	vector<float> vec_out;
	vec_out.push_back(vec_in[0] * cos(phi) + vec_in[2] * sin(phi));
	vec_out.push_back(vec_in[1]);
	vec_out.push_back(-vec_in[0] * sin(phi) + vec_in[2] * cos(phi));
	return vec_out;
}
////////
vector<float> vector_turn_z(vector<float> vec_in, float psi) {
	vector<float> vec_out;
	vec_out.push_back(vec_in[0] * cos(psi) - vec_in[1] * sin(psi));
	vec_out.push_back(vec_in[0] * sin(psi) + vec_in[1] * cos(psi));
	vec_out.push_back(vec_in[2]);
	return vec_out;
}

////////
float angle_to_pi_interval(float angle) {
	while (angle >= PI) {
		angle = angle - 2 * PI;
	}
	while (angle < -PI) {
		angle = angle + 2 * PI;
	}
	return angle;
}
////////
float angle_to_req_angle_interval(float angle, float req_angle) {
	while (angle >= req_angle + PI) {
		angle = angle - 2 * PI;
	}
	while (angle < req_angle - PI) {
		angle = angle + 2 * PI;
	}
	return angle;
}

////////
float** matrix_creat_float(int row_num, int column_num) {
	float** matrix = new float*[row_num]; //allocate space for row pointers
	for (int i = 0; i < row_num; i++) {
		matrix[i] = new float[column_num];
	} //allocate space for each row
	return matrix;
}
////////
int** matrix_creat_int(int row_num, int column_num) {
	int** matrix = new int*[row_num]; //allocate space for row pointers
	for (int i = 0; i < row_num; i++) {
		matrix[i] = new int[column_num];
	} //allocate space for each row
	return matrix;
}

////////
void matrix_delete_float(float** mat, int row_num) {
	for (int i = 0; i < row_num; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}
////////
void matrix_delete_int(int** mat, int row_num) {
	for (int i = 0; i < row_num; i++) {
		delete[] mat[i];
	}
	delete[] mat;
}

////////
void matrix_print_float(int row_num, int column_num, float** mat) {
	for (int r = 0; r < row_num; r++) {
		for (int c = 0; c < column_num; c++)
			cout << mat[r][c] << "\t";
		cout << endl;
	}
	cout << endl;
}
////////
void matrix_print_int(int row_num, int column_num, int** mat) {
	for (int r = 0; r < row_num; r++) {
		for (int c = 0; c < column_num; c++)
			cout << mat[r][c] << " ";
		cout << endl;
	}
	cout << endl;
}

////////
void matrix_print_float_to_file(int row_num, int column_num, float** mat,
		ofstream& results_file) {
	for (int r = 0; r < row_num; r++) {
		for (int c = 0; c < column_num; c++)
			results_file << mat[r][c] << "\t";
		results_file << endl;
	}
	results_file << endl;
}
////////
void matrix_print_int_to_file(int row_num, int column_num, int** mat,
		ofstream& results_file) {
	for (int r = 0; r < row_num; r++) {
		for (int c = 0; c < column_num; c++)
			results_file << mat[r][c] << " ";
		results_file << endl;
	}
	results_file << endl;
}

////////
void vector_print_float(int vector_size, vector<float> vector) {
	for (int c = 0; c < vector_size; c++) {
		cout << vector[c] << "\t";
	}
	cout << endl;

}
////////
void vector_print_int(int vector_size, vector<int> vector) {
	for (int c = 0; c < vector_size; c++) {
		cout << vector[c] << "\t";
	}
	cout << endl;

}

////////
void vector_print_float_to_file(int vector_size, vector<float> vector,
		ofstream& results_file) {
	for (int c = 0; c < vector_size; c++) {
		results_file << vector[c] << "\t";
	}
	results_file << endl;

}
////////
void vector_print_int_to_file(int vector_size, vector<int> vector,
		ofstream& results_file) {
	for (int c = 0; c < vector_size; c++) {
		results_file << vector[c] << "\t";
	}
	results_file << endl;

}

////////
void matrix_print_static(int row_num, int column_num, float mat[][plain_col]) {
	for (int r = 0; r < row_num; r++) {
		for (int c = 0; c < column_num; c++)
			cout << mat[r][c] << " ";
		cout << endl;
	}
	cout << endl;
}

////////
void insert_vector_to_matrix_float(float** matrix, int matrix_row_num,
		vector<float> vec, int wanted_column) {
	for (int r = 0; r < matrix_row_num; r++) {
		matrix[r][wanted_column] = vec[r];
	}
}
////////
void insert_vector_to_matrix_int(int** matrix, int matrix_row_num,
		vector<int> vec, int wanted_column) {
	for (int r = 0; r < matrix_row_num; r++) {
		matrix[r][wanted_column] = vec[r];
	}
}

////////
vector<float> export_vector_from_matrix(float** matrix, int matrix_row_num,
		int wanted_column) {
	vector<float> vec;
	for (int r = 0; r < matrix_row_num; r++) {
		vec.push_back(matrix[r][wanted_column]);
	}
	return vec;
}
////////

float* invert_vector_to_array_float(vector<float> vec_old, int vec_size) {
	float* vec_new = new float[vec_size]; //allocate space for row pointers

	for (int r = 0; r < vec_size; r++) {
		vec_new[r] = vec_old[r];
	}
	return vec_new;
}

////////
void matrix_copy_float(int row_num, int column_num, float** a, float** b) {
	for (int r = 0; r < row_num; r++) {
		for (int c = 0; c < column_num; c++) {
			a[r][c] = b[r][c];
		}
	}
}
////////
void matrix_copy_int(int row_num, int column_num, int** a, int** b) {
	for (int r = 0; r < row_num; r++) {
		for (int c = 0; c < column_num; c++) {
			a[r][c] = b[r][c];
		}
	}
}
////////
vector<float> vector_copy_float(int vector_size, vector<float> a) {
	vector<float> b;
	for (int i = 0; i < vector_size; i++) {
		b.push_back(a[i]);
	}
	return b;
}
////////
vector<int> vector_copy_int(int vector_size, vector<int> a) {
	vector<int> b;
	for (int i = 0; i < vector_size; i++) {
		b.push_back(a[i]);
	}
	return b;
}
////////

////////
float** matrix_add_matrix(float** matrix_first, int matrix_first_row,
		int matrix_first_column, float** matrix_second, int matrix_second_row,
		int matrix_second_column, int matrix_column, float** matrix_joint) {
	matrix_joint = new float*[matrix_first_row + matrix_second_row]; //allocate space for row pointers
	for (int i = 0; i < matrix_first_row + matrix_second_row; i++) {
		matrix_joint[i] = new float[matrix_column];
	} //allocate space for each row

	for (int i = 0; i < matrix_first_row + matrix_second_row; i++) {
		for (int j = 0; j < matrix_column; j++) {
			if ((i < matrix_first_row) && (j < matrix_first_column))
				matrix_joint[i][j] = matrix_first[i][j];
			else if ((i >= matrix_first_row) && (j < matrix_second_column))
				matrix_joint[i][j] = matrix_second[i - matrix_first_row][j];
			else
				matrix_joint[i][j] = -999;
		}
	}

	return matrix_joint;

}
////////

vector<int> read_solution_file_sizes() {
	string line;
	int size;
	vector<int> vect;

	ifstream results_file("results_file_1000000_50_1.txt");
	if (!results_file.is_open())
		results_file.open("results_file_1000000_50_1.txt");
	getline(results_file, line); // soluitons sizes values are:
//		  cout << line << '\n';
	getline(results_file, line); // soluitons matrix for relaxed dubins size is:
//		  cout << line << '\n';
	results_file >> size;
	getline(results_file, line); // size
	//	      cout << size << '\n';
	vect.push_back(size);

	getline(results_file, line); // Dubins markers matrix for relaxed dubins size is:
//		  cout << line << '\n';
	results_file >> size;
	getline(results_file, line); // size
	//	      cout << size << '\n';
	vect.push_back(size);

	getline(results_file, line); // soluitons matrix for full dubins size is:
//		  cout << line << '\n';
	results_file >> size;
	getline(results_file, line); // size
	//	      cout << size << '\n';
	vect.push_back(size);

	getline(results_file, line); // Dubins markers matrix for full dubins size is:
//		  cout << line << '\n';
	results_file >> size;
	getline(results_file, line); // size
	//	      cout << size << '\n';
	vect.push_back(size);

	results_file.close();

	return vect;
}

////////
void read_solution_file(float** solutions_states_relaxed_Dubins,
		float** solutions_marks_relaxed_Dubins,
		float** solutions_states_full_Dubins,
		float** solutions_marks_full_Dubins,
		float** solutions_TTH_Relaxed_Dubins_matrix,
		float** solutions_TTH_full_Dubins_matrix, float** vectors_matrix_temp,
		vector<int> sizes_values) {
	string line;
	int size;
	float val;

	ifstream results_file("results_file_1000000_50_1.txt");
	if (!results_file.is_open())
		results_file.open("results_file_1000000_50_1.txt");

	// !!! read solutions_states_relaxed_Dubins matrix !!! //
	getline(results_file, line); // soluitons sizes values are:
	getline(results_file, line); // soluitons matrix for relaxed dubins size is:
	getline(results_file, line);
	getline(results_file, line); // Dubins markers matrix for relaxed dubins size is:
	getline(results_file, line);
	getline(results_file, line); // soluitons matrix for full dubins size is:
	getline(results_file, line);
	getline(results_file, line); // Dubins markers matrix for full dubins size is:
	//		  cout << line << '\n';
	getline(results_file, line);
	getline(results_file, line); //
	getline(results_file, line); //

	getline(results_file, line); // soluitons matrix for relaxed dubins is:
	for (int i = 0; i < monte_carlo_number * states_number; i++) {
		for (int j = 0; j < sizes_values[0]; j++) {
			results_file >> val;
			solutions_states_relaxed_Dubins[i][j] = val;
		}
		getline(results_file, line);
	}
//	      	  	  	  	  matrix_print_float(monte_carlo_number*states_number, size,solutions_states_relaxed_Dubins);
	getline(results_file, line);
	getline(results_file, line);

	// !!! read solutions_marks_relaxed_Dubins matrix !!! //
	getline(results_file, line); // Dubins markers matrix for relaxed dubins is:
	for (int i = 0; i < monte_carlo_number; i++) {
		for (int j = 0; j < sizes_values[1]; j++) {
			results_file >> val;
			solutions_marks_relaxed_Dubins[i][j] = val;
		}
		getline(results_file, line);
	}
	getline(results_file, line);
	getline(results_file, line);

	// !!! read solutions_states_full_Dubins matrix !!! //
	getline(results_file, line); // soluitons matrix for full dubins is:
	for (int i = 0; i < monte_carlo_number * states_number; i++) {
		for (int j = 0; j < sizes_values[2]; j++) {
			results_file >> val;
			solutions_states_full_Dubins[i][j] = val;
		}
		getline(results_file, line);
	}
//	      	  	  	  	  matrix_print_float(states_number, size,mat);
	getline(results_file, line);
	getline(results_file, line);

	// !!! read solutions_marks_full_Dubins matrix !!! //
	getline(results_file, line); // Dubins markers matrix for full dubins is:
	for (int i = 0; i < monte_carlo_number; i++) {
		for (int j = 0; j < sizes_values[3]; j++) {
			results_file >> val;
			solutions_marks_full_Dubins[i][j] = val;
		}
		getline(results_file, line);
	}

	getline(results_file, line); //
	getline(results_file, line); //
	getline(results_file, line); //

	vector<float> monte_carlo_avg_Relaxed_Dubins;
	vector<float> monte_carlo_avg_full_Dubins;
	vector<float> convergence_values;

	getline(results_file, line); // convergence matrix for Relaxed Dubins is:
	for (int i = 0; i < monte_carlo_number; i++) {
		for (int j = 0; j < convergence_size; j++) {
			results_file >> val;
			solutions_TTH_Relaxed_Dubins_matrix[i][j] = val;
		}
		getline(results_file, line);
	}

	getline(results_file, line); //
	getline(results_file, line); // number of vertices in the tree
	for (int j = 0; j < convergence_size; j++) {
		results_file >> val;
		convergence_values.push_back(val);
	}
	getline(results_file, line);

	getline(results_file, line); // avrages results
	for (int j = 0; j < convergence_size; j++) {
		results_file >> val;
		monte_carlo_avg_Relaxed_Dubins.push_back(val);
	}
	getline(results_file, line);

	getline(results_file, line); //
	getline(results_file, line); //

	getline(results_file, line); // convergence matrix for full Dubins is:
	for (int i = 0; i < monte_carlo_number; i++) {
		for (int j = 0; j < convergence_size; j++) {
			results_file >> val;
			solutions_TTH_full_Dubins_matrix[i][j] = val;
		}
		getline(results_file, line);
	}

	getline(results_file, line); //
	getline(results_file, line); // number of vertices in the tree
	for (int j = 0; j < convergence_size; j++) {
		results_file >> val;
//		  convergence_values.push_back(val); // assuming same as before
	}
	getline(results_file, line);

	getline(results_file, line); // avrages results
	for (int j = 0; j < convergence_size; j++) {
		results_file >> val;
		monte_carlo_avg_full_Dubins.push_back(val);
	}
	results_file.close();

	for (int i = 0; i < convergence_size; i++) {
		vectors_matrix_temp[0][i] = monte_carlo_avg_Relaxed_Dubins[i];
		vectors_matrix_temp[1][i] = monte_carlo_avg_full_Dubins[i];
		vectors_matrix_temp[2][i] = convergence_values[i];
	}
}

////////
float dis2(vector<float> a, vector<float> b) {
	float sumpower2 = 0;
	//	vector_print_float(states_number, a);
	//	vector_print_float(states_number, b);

	for (int i = 0; i < 2; i++) {
		sumpower2 = sumpower2 + pow((a[i] - b[i]), 2);
	}
	float dist2 = sqrt(sumpower2);
	return dist2;
}

////////
float dis3(vector<float> a, vector<float> b) {

	int angle_state_mone = 0;
	float sumpower2 = 0;
	for (int i = 0; i < 3; i++) {
		if ((angle_state_mone < angle_state_vector_size)
				&& (i == angle_state_vector[angle_state_mone] - 1)) {
			float dis_ab = a[i] - b[i];
			if (fabs(dis_ab) > PI)
				dis_ab = 2 * PI - fabs(dis_ab);
			sumpower2 = sumpower2 + pow(dis_ab, 2);
			angle_state_mone++;
		} else
			sumpower2 = sumpower2 + pow((a[i] - b[i]), 2);
	}
	float dist3 = sqrt(sumpower2);
	return dist3;
}

////////
float dis_state_number(vector<float> a, vector<float> b) {

	int angle_state_mone = 0;
	float sumpower2 = 0;
	for (int i = 0; i < states_number; i++) {
		if ((angle_state_mone < angle_state_vector_size)
				&& (i == angle_state_vector[angle_state_mone] - 1)) {
			float dis_ab = a[i] - b[i];
			if (fabs(dis_ab) > PI)
				dis_ab = 2 * PI - fabs(dis_ab);
			sumpower2 = sumpower2 + pow(dis_ab, 2);
			angle_state_mone++;
		} else
			sumpower2 = sumpower2 + pow((a[i] - b[i]), 2);
	}
	float dis_state_num = sqrt(sumpower2);
	return dis_state_num;
}

////////
float cost_func(vector<float> path_marks) {

	// quadratic acceleration value cost function

	// path_marks_and_vector= {	R, c, L,																	--- 0-2
	//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
	//							v_norm,																		--- 7
	//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
	//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

	float cost = 0;

	if (path_marks[0] == -999)
		return -999;

	float R = path_marks[0];
	float control_value = pow(V_max, 2) / R;

	float L_tot = path_marks[2];
	float initialArcAngle = path_marks[5];
	float L_first_circle = initialArcAngle * R;

	if (L_tot <= L_first_circle) {
		cost = pow(control_value, 2) * L_tot / V_max;
		return cost;
	}

	float L_straight = path_marks[7];

	if ((L_tot > L_first_circle) && (L_tot < L_first_circle + L_straight)) {
		cost = pow(control_value, 2) * L_first_circle / V_max;
		return cost;
	}

	float finalArcAngle = path_marks[10];
	float L_last_circle = finalArcAngle * R;

	if ((L_tot > L_first_circle + L_straight)
			&& (L_tot < L_first_circle + L_straight + L_last_circle)) {
		cost = pow(control_value, 2) * (L_tot - L_straight) / V_max;
		return cost;
	} else
		cost = pow(control_value, 2) * (L_first_circle + L_last_circle) / V_max;

	return cost;
}

////////
float vertex_cost_value(int vertex_index, int** edges, float** traj) {
	float cost;
	cost = 0;

	int temp_index = vertex_index;
	int prev_index, end_index;

	while (temp_index != 0) {
		prev_index = edges[0][temp_index - 1];
		//		cout << traj[0][temp_index - 1] << endl;
		cost = cost + traj[0][temp_index - 1];
		temp_index = prev_index;
	}

	return cost;
}

////////
int path_size_find(float** vertices, int** edges, int solution_index) {
	int path_size = 1;
	int prev_index, temp_index;
	temp_index = solution_index;

	while (temp_index != 0) {
		prev_index = edges[0][temp_index - 1];
		temp_index = prev_index;
		path_size++;
	}

	return path_size;
}

////////
int path_extraction(float** soluiton_matrix, float** vertices, int** edges,
		int solution_index, int path_size) {
	int prev_index, temp_index;
	temp_index = solution_index;

	vector<float> x1;
	vector<float> x2;
	vector<float> x3;
	vector<float> x4;

	x1.push_back(vertices[0][solution_index]);
	x2.push_back(vertices[1][solution_index]);
	x3.push_back(vertices[2][solution_index]);
	x4.push_back(vertices[3][solution_index]);

	while (temp_index != 0) {
		prev_index = edges[0][temp_index - 1];

		x1.insert(x1.begin(), vertices[0][prev_index]);
		x2.insert(x2.begin(), vertices[1][prev_index]);
		x3.insert(x3.begin(), vertices[2][prev_index]);
		x4.insert(x4.begin(), vertices[3][prev_index]);
//		if  (solution_index==5549)
//		{
//			cout << "prev_index:" << prev_index << endl;
//			cout << vertices[0][prev_index]<< "," << vertices[1][prev_index]<< "," << vertices[2][prev_index]<< "," << vertices[3][prev_index] << endl;
//		}
		temp_index = prev_index;
	}

	for (int i = 0; i < path_size; i++) {
		soluiton_matrix[0][i] = x1[i];
		soluiton_matrix[1][i] = x2[i];
		soluiton_matrix[2][i] = x3[i];
		soluiton_matrix[3][i] = x4[i];
	}
//					matrix_print_float(states_number, path_size, soluiton_matrix);

	return path_size;
}

////////
int dubins_mark_extraction(float** soluiton_matrix, float** vertices,
		int** edges, float** traj, int solution_index, int path_size) {
	int prev_index, temp_index;
	temp_index = solution_index;

	vector<int> indexes;

	indexes.push_back(solution_index);
//				cout << "step 3.1 done" << endl;
//				cout << "temp_index is:" << temp_index << endl;
//	int sign=0;
//	int sign_mone=0;
//	if (temp_index==166939)
//		sign=1;
	prev_index = edges[0][temp_index - 1];

	while (prev_index != 0) {
//		cout << "step 3.1.1 done" << endl;
//		if 	(sign==1)
//		{
//			sign_mone++;
//			if (sign_mone>27)
//				////////
//			{
//				matrix_print_int_to_file(edges_row, edges_col, edges,results_file);
//				results_file.close();
//			}
//			cout << "temp_index is:" << temp_index << endl;
//			cout << "temp_index is:" << edges[0][temp_index - 1] << endl;
//		}

//		cout << "step 3.1.2 done" << endl;
		indexes.insert(indexes.begin(), prev_index);
//		cout << "step 3.1.3 done" << endl;
		temp_index = prev_index;
//		cout << "step 3.1.4 done" << endl;
		prev_index = edges[0][temp_index - 1];
	}
//				cout << "step 3.2 done" << endl;

//	if (solution_index==5549) {
//		vector_print_int(path_size - 1,indexes);
//		}
//				cout << "step 3.3 done" << endl;
//				cout << "outer loop size:" << path_size - 1 << endl;
//				cout << "inner loop size:" << dubins_mark_size << endl;
//			vector_print_int(path_size - 1,indexes);

	for (int i = 0; i < path_size - 1; i++) {
		for (int j = 0; j < dubins_mark_size; j++) {
			soluiton_matrix[0][i * dubins_mark_size + j] =
					traj[1 + j][indexes[i] - 1];
//			cout << "traj point:" << traj[1 + j][indexes[i] - 1] << endl;
		}
	}
//				cout << "step 3.4 done" << endl;

//					matrix_print_float(1, (path_size-1)*dubins_mark_size, soluiton_matrix);

	return path_size;
}

////////
void target_sim(float** goal_path, int steps, vector<float> x_0, float V,
		float dt, float heading) {
	vector<float> path_next(states_number);
	insert_vector_to_matrix_float(goal_path, states_number, x_0, 0);
	//matrix_print_float(states_number, 2, goal_path);

	for (int i = 1; i < steps; i++) {
		path_next[0] = goal_path[0][i - 1] + V * dt * cos(heading);
		path_next[1] = goal_path[1][i - 1] + V * dt * sin(heading);
		path_next[2] = goal_path[2][i - 1] + dt;

		insert_vector_to_matrix_float(goal_path, states_number, path_next, i);
	}
}

////////
void RRT_function(int steps, float** vertices, int** edges, float** traj);
////////
void RRT_star_function(int steps, float** vertices, int** edges, float** traj);
////////
vector<float> Dubins_state_time_steer(vector<float> vector_start,
		vector<float> vector_end);
////////
float dubins(vector<float> pinit, float thinit, vector<float> pf, float thf,
		float R);
vector<float> dubins_LSL(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R);
vector<float> dubins_LSR(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R);
vector<float> dubins_RSL(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R);
vector<float> dubins_RSR(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R);

////////
int goal_point_search_Dubins(int steps, float** vertices, int** edges,
		vector<float> goal_point) {
	float time_best = plain[states_number - 1][1];
	int i = 0;
	int index = -999;
	int *index_p;
	float time_new, distance_2d_from_target;
	//		vector_print_float(states_number, vector_new);
	if (kdtree_mode == 2) {
		struct kdtree *kd_vertices_goal_search = kd_create(states_number - 1);
		float time_new;
		for (int mone = 0; mone < steps + 1; mone++) {

			float vector_temp[states_number - 1] = { vertices[0][mone],
					vertices[1][mone], vertices[2][mone] };
			float *vector_temp_p;
			vector_temp_p = vector_temp;
			//			cout << "vector temp: " << vector_temp[0] << ' ' << vector_temp[1] << ' ' << vector_temp[2] << endl;
			kd_insertf(kd_vertices_goal_search, vector_temp_p,
					&edges[1][mone - 1]);
			vector_temp_p[2] = vector_temp_p[2] + 2 * PI;
			kd_insertf(kd_vertices_goal_search, vector_temp_p,
					&edges[1][mone - 1]);
			vector_temp_p[2] = vector_temp_p[2] - 4 * PI;
			kd_insertf(kd_vertices_goal_search, vector_temp_p,
					&edges[1][mone - 1]);

		}
		struct kdres *hit_points_set;
		float vector_new[states_number - 1];
		int *temp_index;
		//		cout << "goal point: " << goal_point[0] << ' ' << goal_point[1] << endl;
		float goal_area[states_number - 1] = { goal_point[0], goal_point[1],
				goal_point[2] };
//		float *goal_area_p;
//		goal_area_p = goal_area;
		hit_points_set = kd_nearest_rangef(kd_vertices_goal_search, goal_area,
				hit_dis);
		while (!kd_res_end(hit_points_set)) {
			temp_index = (int*) kd_res_itemf(hit_points_set, vector_new);
			//			 cout << "temp index: " << temp_index << endl;
			//			 cout << "vector new: " << vector_new[0] << ' ' << vector_new[1] << endl;
			time_new = vertices[3][*temp_index];
			//			 cout << "time new: " << time_new << endl;
			if (time_new < time_best) {
				time_best = time_new;
				//					 cout << "time best: " << time_best << endl;
				*index_p = *temp_index;
				//					 cout << "temp index: " << temp_index << endl;
				//					 cout << "index p: " << index_p << endl;
				index = *index_p;
				//					 cout << "index: " << index << endl;

			}
			//			 cout << index << endl;
			kd_res_next(hit_points_set);

		}
		kd_res_free(hit_points_set);
		kd_free(kd_vertices_goal_search);

		return index;
	} else {
		do {
			vector<float> vector_new = export_vector_from_matrix(vertices,
					states_number, i);
			//vector_print_float(states_number, vector_new);
			//vector_print_float(states_number, goal_point);

			time_new = vector_new[3];
			distance_2d_from_target = dis3(vector_new, goal_point);

			if ((time_new < time_best) && (distance_2d_from_target < hit_dis)) {
				time_best = time_new;
				index = i;
			}
			i++;
		} while (i < steps + 1);

		return index;
	}
	delete[] index_p;
}

////////
int goal_point_search_Relaxed_Dubins(int steps, float** vertices, int** edges,
		vector<float> goal_point) {
	float time_best = plain[states_number - 1][1];
	int i = 0;
	int index = -999;
	int *index_p;
	float time_new, distance_2d_from_target;
	//		vector_print_float(states_number, vector_new);
	if (kdtree_mode == 2) {
		struct kdtree *kd_vertices_goal_search = kd_create(states_number - 2);
		float time_new;
		for (int mone = 0; mone < steps + 1; mone++) {

			float vector_temp[states_number - 2] = { vertices[0][mone],
					vertices[1][mone] };
			float *vector_temp_p;
			vector_temp_p = vector_temp;
			//			cout << "vector temp: " << vector_temp[0] << ' ' << vector_temp[1] << endl;
			kd_insertf(kd_vertices_goal_search, vector_temp_p,
					&edges[1][mone - 1]);
		}
		struct kdres *hit_points_set;
		float vector_new[states_number - 2];
		int *temp_index;
		//		cout << "goal point: " << goal_point[0] << ' ' << goal_point[1] << endl;
		float goal_area[states_number - 2] = { goal_point[0], goal_point[1] };
		float *goal_area_p;
		goal_area_p = goal_area;
		hit_points_set = kd_nearest_rangef(kd_vertices_goal_search, goal_area,
				hit_dis_RD);
		while (!kd_res_end(hit_points_set)) {
			temp_index = (int*) kd_res_itemf(hit_points_set, vector_new);
			//			 cout << "temp index: " << temp_index << endl;
			//			 cout << "vector new: " << vector_new[0] << ' ' << vector_new[1] << endl;
			time_new = vertices[3][*temp_index];
			//			 cout << "time new: " << time_new << endl;
			if (time_new < time_best) {
				time_best = time_new;
				//					 cout << "time best: " << time_best << endl;
				*index_p = *temp_index;
				//					 cout << "temp index: " << temp_index << endl;
				//					 cout << "index p: " << index_p << endl;
				index = *index_p;
				//					 cout << "index: " << index << endl;

			}
			//			 cout << index << endl;
			kd_res_next(hit_points_set);
		}
		kd_res_free(hit_points_set);
		kd_free(kd_vertices_goal_search);

		return index;
	} else {
		do {
			vector<float> vector_new = export_vector_from_matrix(vertices,
					states_number, i);
			//vector_print_float(states_number, vector_new);
			//vector_print_float(states_number, goal_point);

			time_new = vector_new[3];
			distance_2d_from_target = dis2(vector_new, goal_point);

			if ((time_new < time_best)
					&& (distance_2d_from_target < hit_dis_RD)) {
				time_best = time_new;
				index = i;
			}
			i++;
		} while (i < steps + 1);

		return index;
	}
	delete[] index_p;
}

////////
float vector_avrage(vector<float> vector, int vector_size) {
	int div_size = vector_size;
	float accumulator = 0;
	for (int i = 0; i < vector_size; i++) {
		if (vector[i] == -999)
			div_size--;
		else
			accumulator = accumulator + vector[i];
	}
	if (div_size == 0) {
		return -999;
	}
	return accumulator / div_size;
}

////////
vector<float> matrix_column_avrages(float** matrix, int row_num, int col_num) {
	vector<float> avrages_vector(col_num);
	for (int j = 0; j < col_num; j++) {
		int div_size = row_num;
		float accumulator = 0;
		for (int i = 0; i < row_num; i++) {
			if (matrix[i][j] == -999)
				div_size--;
			else
				accumulator = accumulator + matrix[i][j];
		}
		if (div_size == 0) {
			avrages_vector[j] = -999;
		}
		avrages_vector[j] = accumulator / div_size;
	}
	return avrages_vector;
}

float** matrix_indexes_to_values(int** matrix_indexes, int matrix_indexes_row,
		int matrix_indexes_column, float** vertices) {
	float** matrix_values = matrix_creat_float(matrix_indexes_row,
			matrix_indexes_column);
	int temp;
	for (int i = 0; i < matrix_indexes_row; i++) {
		for (int j = 0; j < matrix_indexes_column; j++) {
			if (matrix_indexes[i][j] == -999)
				matrix_values[i][j] = -999;
			else
				temp = matrix_indexes[i][j];
			matrix_values[i][j] = vertices[2][temp];
		}
	}
	return matrix_values;
}

///////////////////////////////////////////////////////////////////
int main_RRT_sim(float** vertices_final, int** edges_final,
		float** traj_final) {

	hit_index = -999;
	hit_index_RD = -999;
	//const int steps = 500;
	vector<int> convergence_values;

//	convergence_values.push_back(1000);
//	convergence_values.push_back(1259);
//	convergence_values.push_back(1585);
//	convergence_values.push_back(1995);
//	convergence_values.push_back(2512);
//	convergence_values.push_back(3162);
//	convergence_values.push_back(3981);
//	convergence_values.push_back(5012);
//	convergence_values.push_back(6310);
//	convergence_values.push_back(7943);
	convergence_values.push_back(10000);
	convergence_values.push_back(12589);
	convergence_values.push_back(15849);
	convergence_values.push_back(19953);
	convergence_values.push_back(25119);
	convergence_values.push_back(31623);
	convergence_values.push_back(39811);
	convergence_values.push_back(50119);
	convergence_values.push_back(63096);
	convergence_values.push_back(79433);
	convergence_values.push_back(100000);
	convergence_values.push_back(125893);
	convergence_values.push_back(158489);
	convergence_values.push_back(199526);
	convergence_values.push_back(251189);
	convergence_values.push_back(316228);
	convergence_values.push_back(398107);
	convergence_values.push_back(501187);
	convergence_values.push_back(630957);
	convergence_values.push_back(794328);
	convergence_values.push_back(1000000);

	//	cout << convergence_values[1] << '\t' << convergence_values[2] << endl;
	//	vector_print_int(`, convergence_values);

	int monte_carlo_steps;
	if (solution_mode == 1) {
		monte_carlo_steps = 1;
		convergence_size = 1;
	} else
		monte_carlo_steps = monte_carlo_number;

	int** solutions_TTH_Dubins_matrix_indexes = matrix_creat_int(
			monte_carlo_steps, convergence_size); // TTH - Time To Hit
	float** solutions_TTH_Dubins_matrix = matrix_creat_float(monte_carlo_steps,
			convergence_size); // TTH - Time To Hit
	int** solutions_TTH_Relaxed_Dubins_matrix_indexes = matrix_creat_int(
			monte_carlo_steps, convergence_size); // TTH - Time To Hit
	float** solutions_TTH_Relaxed_Dubins_matrix = matrix_creat_float(
			monte_carlo_steps, convergence_size); // TTH - Time To Hit
	int path_size;
	//for (int steps_mone = 0; steps_mone < steps; steps_mone++)
	//{
	//vector<float> monte_carlo_TTH(monte_carlo_steps);

	// defining the initial vertex

	vector<float> x_0;
	x_0.push_back(x1_0);
	x_0.push_back(x2_0);
	x_0.push_back(x3_0);
	x_0.push_back(x4_0);

	const int steps_target = 400;

	float x1_0_target = 15;
	float x2_0_target = 0;
	float x3_0_target = 90 * d2r;
	float x4_0_target = 0;

	vector<float> goal_point;
	goal_point.push_back(goal_point_x1);
	goal_point.push_back(goal_point_x2);
	goal_point.push_back(goal_point_x3);
	goal_point.push_back(goal_point_x4); // fourth number for time has no meaning has no meaning

	int temp;
	//	float space_volume=1;
	//	for (int dim_mone=0;dim_mone<states_number;dim_mone++)
	//	{
	//		space_volume =space_volume*(plain[dim_mone][1]-plain[dim_mone][0]);
	//	}

	int vertices_col = steps + 1;
	int edges_col = steps;
	int traj_col = steps;
	float** vertices; // = matrix_creat_float(states_number, vertices_col);
	int** edges; // = matrix_creat_int(edges_row, edges_col);
	float** traj;
	float** solutions_states_full_Dubins;
	float** solutions_marks_full_Dubins;
	float** solutions_states_relaxed_Dubins;
	float** solutions_marks_relaxed_Dubins;
	float** solutions_states_full_Dubins_new;
	float** solutions_marks_full_Dubins_new;
	float** solutions_states_relaxed_Dubins_new;
	float** solutions_marks_relaxed_Dubins_new;

	float** soluiton_states_matrix;
	float** soluiton_marks_matrix;

	float** goal_target_path = matrix_creat_float(states_number,
			steps_target + 1);

	vertices = matrix_creat_float(states_number, vertices_col);
	edges = matrix_creat_int(edges_row, edges_col);
	traj = matrix_creat_float(traj_row, traj_col);

	int path_size_full_dubins, path_size_max_full_dubins;
	int path_size_relaxed_dubins, path_size_max_relaxed_dubins;
	int matrix_column;

	for (int i_m_c = 0; i_m_c < monte_carlo_steps; i_m_c++) {
		cout << "Now starting solution number: " << i_m_c + 1 << endl;
		cout << "Path assinment method is: " << steer_length_mode << endl;

		//initializing vertices matrix

		if (solution_mode == 2) {
			vertices = matrix_creat_float(states_number, vertices_col);
			edges = matrix_creat_int(edges_row, edges_col);
			traj = matrix_creat_float(traj_row, traj_col);
		}
		insert_vector_to_matrix_float(vertices, states_number, x_0, 0);
		//cout<<"This is part 1.";

		if (run_mode == 1) {
			RRT_function(steps, vertices, edges, traj);
		}
		if (run_mode == 2) {
			RRT_star_function(steps, vertices, edges, traj);
		}
//		matrix_print_float(states_number, vertices_col, vertices);
//		matrix_print_int(edges_row, edges_col, edges);
//		matrix_print_float(traj_row, traj_col, traj);

//		vector<float> print_vertices = export_vector_from_matrix(vertices, states_number, 6282);
//		vector_print_float(states_number, print_vertices);
//		cout << edges[0][6281] << ' '<< edges[1][6281] << endl;
//		vector<float> print_traj = export_vector_from_matrix(traj, traj_row, 6282);
//		vector_print_float(traj_row, print_traj);

		if (goal_mode == 2) {
			float goal_area[states_number][plain_col] = { { 18, 20 },
					{ 18, 20 }, { 30, 35 } };
		}

		if (goal_mode == 3) {

			vector<float> x_0_target;
			x_0_target.push_back(x1_0_target);
			x_0_target.push_back(x2_0_target);
			x_0_target.push_back(x3_0_target);

			float target_heading = 90 * d2r;
			float t_0_target = x3_0_target;
			const float dt_target = 0.1;
			const float V_target = 0.5;

			target_sim(goal_target_path, steps_target + 1, x_0_target, V_target,
					dt_target, target_heading);
			//matrix_print_float(states_number, steps_target + 1, goal_target_path);
		}

		if ((goal_mode == 1) && (solution_mode == 1)) {

			goal_point.push_back(goal_point_x1);
			goal_point.push_back(goal_point_x2);
			goal_point.push_back(goal_point_x3);
			goal_point.push_back(goal_point_x4); // in this simulation the forth number has no meaning
			hit_index = goal_point_search_Dubins(steps, vertices, edges,
					goal_point);
			hit_index_RD = goal_point_search_Relaxed_Dubins(steps, vertices,
					edges, goal_point);

			if (hit_index_RD != -999) {
				cout << "Solution for Relaxed Dubins" << endl;
				cout << "The hit index is:" << hit_index_RD << endl;
				cout << "The vertex is:" << endl;
				vector<float> ver_temp;
				ver_temp = export_vector_from_matrix(vertices, states_number,
						hit_index_RD);
				vector_print_float(states_number, ver_temp);
			}
			if (hit_index != -999) {
				cout << "Solution for Dubins" << endl;
				cout << "The hit index is:" << hit_index << endl;
				cout << "The vertex is:" << endl;
				vector<float> ver_temp;
				ver_temp = export_vector_from_matrix(vertices, states_number,
						hit_index);
				vector_print_float(states_number, ver_temp);

			}
			if ((hit_index_RD == -999) && (hit_index == -999)) {
				cout << "No solution was found!" << endl;
				//monte_carlo_TTH[i_m_c] = -999;
			}
		}
		if ((goal_mode == 1) && (solution_mode == 2)) {

			goal_point[0] = goal_point_x1;
			goal_point[1] = goal_point_x2;
			goal_point[2] = goal_point_x3;
			goal_point[3] = goal_point_x4; // in this simulation the forth number has no meaning

			for (int i_conv = 0; i_conv < convergence_size; i_conv++) {

				//					vector_print_int(convergence_size, convergence_values);
				temp = goal_point_search_Dubins(convergence_values[i_conv],
						vertices, edges, goal_point);
				solutions_TTH_Dubins_matrix_indexes[i_m_c][i_conv] = temp;
				if (temp == -999)
					solutions_TTH_Dubins_matrix[i_m_c][i_conv] = -999;
				else
					solutions_TTH_Dubins_matrix[i_m_c][i_conv] =
							vertices[3][temp];
			}

			for (int i_conv = 0; i_conv < convergence_size; i_conv++) {

				//					vector_print_int(convergence_size, convergence_values);
				temp = goal_point_search_Relaxed_Dubins(
						convergence_values[i_conv], vertices, edges,
						goal_point);
				solutions_TTH_Relaxed_Dubins_matrix_indexes[i_m_c][i_conv] =
						temp;
				if (temp == -999)
					solutions_TTH_Relaxed_Dubins_matrix[i_m_c][i_conv] = -999;
				else
					solutions_TTH_Relaxed_Dubins_matrix[i_m_c][i_conv] =
							vertices[states_number - 1][temp];
			}

			temp = goal_point_search_Relaxed_Dubins(steps, vertices, edges,
					goal_point);
			path_size = path_size_find(vertices, edges, temp);
			soluiton_states_matrix = matrix_creat_float(states_number,
					path_size);
			soluiton_marks_matrix = matrix_creat_float(1,
					(path_size - 1) * dubins_mark_size);
//							cout << "step 3 done" << endl;

			path_size_relaxed_dubins = path_extraction(soluiton_states_matrix,
					vertices, edges, temp, path_size);
//							matrix_print_float(states_number, path_size_relaxed_dubins, soluiton_states_matrix);

			path_size_relaxed_dubins = dubins_mark_extraction(
					soluiton_marks_matrix, vertices, edges, traj, temp,
					path_size);
//							matrix_print_float(1, (path_size_relaxed_dubins - 1) * dubins_mark_size, soluiton_marks_matrix);

//							cout << "step 4 done" << endl;

//							cout << "path size is:" << path_size_relaxed_dubins << endl;

			if (i_m_c == 0) {
				path_size_max_relaxed_dubins = path_size_relaxed_dubins;
				solutions_states_relaxed_Dubins = matrix_creat_float(states_number, path_size_max_relaxed_dubins);
//							matrix_print_float(states_number, path_size_max_relaxed_dubins, solutions_states_relaxed_Dubins);
//							matrix_print_float(states_number, path_size_max_relaxed_dubins, soluiton_states_matrix);
				matrix_copy_float(states_number, path_size_max_relaxed_dubins,
						solutions_states_relaxed_Dubins,
						soluiton_states_matrix);
//							matrix_print_float(states_number, path_size_max_relaxed_dubins, solutions_states_relaxed_Dubins);

				solutions_marks_relaxed_Dubins = matrix_creat_float(1,((path_size_max_relaxed_dubins - 1) * dubins_mark_size));
				matrix_copy_float(1,
						(path_size_max_relaxed_dubins - 1) * dubins_mark_size,
						solutions_marks_relaxed_Dubins, soluiton_marks_matrix);
			} else {

				matrix_column = path_size_relaxed_dubins;
				if (path_size_max_relaxed_dubins > path_size_relaxed_dubins)
					matrix_column = path_size_max_relaxed_dubins;
				solutions_states_relaxed_Dubins_new = matrix_add_matrix(
						solutions_states_relaxed_Dubins, states_number * i_m_c,
						path_size_max_relaxed_dubins, soluiton_states_matrix,
						states_number, path_size_relaxed_dubins, matrix_column,
						solutions_states_relaxed_Dubins_new);
				matrix_delete_float(solutions_states_relaxed_Dubins,
						states_number * i_m_c);

				solutions_states_relaxed_Dubins = matrix_creat_float(
						states_number * i_m_c + states_number, matrix_column);

				matrix_copy_float(states_number * i_m_c + states_number,
						matrix_column, solutions_states_relaxed_Dubins,
						solutions_states_relaxed_Dubins_new);
				matrix_delete_float(solutions_states_relaxed_Dubins_new,
						states_number * i_m_c + states_number);

				matrix_column = (path_size_relaxed_dubins - 1)
						* dubins_mark_size;
				if (path_size_max_relaxed_dubins > path_size_relaxed_dubins)
					matrix_column = (path_size_max_relaxed_dubins - 1)
							* dubins_mark_size;
				solutions_marks_relaxed_Dubins_new = matrix_add_matrix(
						solutions_marks_relaxed_Dubins,
						i_m_c,(path_size_max_relaxed_dubins - 1) * dubins_mark_size,
						soluiton_marks_matrix,
						1, (path_size_relaxed_dubins - 1) * dubins_mark_size,
						matrix_column, solutions_marks_relaxed_Dubins_new);

				matrix_delete_float(solutions_marks_relaxed_Dubins, i_m_c);

				solutions_marks_relaxed_Dubins = matrix_creat_float(i_m_c + 1,
						matrix_column);
				matrix_copy_float(i_m_c + 1, matrix_column,
						solutions_marks_relaxed_Dubins,
						solutions_marks_relaxed_Dubins_new);
				matrix_delete_float(solutions_marks_relaxed_Dubins_new,
						i_m_c + 1);

				if (path_size_max_relaxed_dubins < path_size_relaxed_dubins)
					path_size_max_relaxed_dubins = path_size_relaxed_dubins;
//							matrix_print_float(states_number*(i_m_c+1), path_size_max_full_dubins, solutions_states_full_Dubins);

			}
			matrix_delete_float(soluiton_states_matrix, states_number);
			matrix_delete_float(soluiton_marks_matrix, 1);

			// Dubins_marks= {	R, c, L,																	--- 0-2
			//				initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
			//				v_norm,																		--- 7
			//				finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11

			// find and add solutions for full dubins problem
			temp = goal_point_search_Dubins(steps, vertices, edges, goal_point);
			path_size = path_size_find(vertices, edges, temp);
			soluiton_states_matrix = matrix_creat_float(states_number,
					path_size);
			soluiton_marks_matrix = matrix_creat_float(1,
					(path_size - 1) * dubins_mark_size);

			path_size_full_dubins = path_extraction(soluiton_states_matrix,
					vertices, edges, temp, path_size);
			path_size_full_dubins = dubins_mark_extraction(
					soluiton_marks_matrix, vertices, edges, traj, temp,
					path_size);
//							matrix_print_float(states_number, path_size_full_dubins, soluiton_states_matrix);
//							cout << "solutions size is: "<< path_size << endl;

			if (i_m_c == 0) {
				path_size_max_full_dubins = path_size_full_dubins;
				solutions_states_full_Dubins = matrix_creat_float(states_number,
						path_size_max_full_dubins);
				solutions_marks_full_Dubins = matrix_creat_float(1,
						(path_size_max_full_dubins - 1) * dubins_mark_size);
//							matrix_print_float(states_number, path_size_max_full_dubins, solutions_states_full_Dubins);
//							matrix_print_float(states_number, path_size_max_full_dubins, soluiton_states_matrix);
				matrix_copy_float(states_number, path_size_max_full_dubins,
						solutions_states_full_Dubins, soluiton_states_matrix);
				matrix_copy_float(1,
						(path_size_max_full_dubins - 1) * dubins_mark_size,
						solutions_marks_full_Dubins, soluiton_marks_matrix);
			} else {

				matrix_column = path_size_full_dubins;
				if (path_size_max_full_dubins > path_size_full_dubins)
					matrix_column = path_size_max_full_dubins;
				solutions_states_full_Dubins_new = matrix_add_matrix(
						solutions_states_full_Dubins, states_number * i_m_c,
						path_size_max_full_dubins, soluiton_states_matrix,
						states_number, path_size_full_dubins, matrix_column,
						solutions_states_full_Dubins_new);
				matrix_delete_float(solutions_states_full_Dubins,
						states_number * i_m_c);

				solutions_states_full_Dubins = matrix_creat_float(
						states_number * i_m_c + states_number, matrix_column);
				matrix_copy_float(states_number * i_m_c + states_number,
						matrix_column, solutions_states_full_Dubins,
						solutions_states_full_Dubins_new);
				matrix_delete_float(solutions_states_full_Dubins_new,
						states_number * i_m_c + states_number);

				matrix_column = (path_size_full_dubins - 1) * dubins_mark_size;
				if (path_size_max_full_dubins > path_size_full_dubins)
					matrix_column = (path_size_max_full_dubins - 1)
							* dubins_mark_size;
				solutions_marks_full_Dubins_new = matrix_add_matrix(
						solutions_marks_full_Dubins, i_m_c,
						(path_size_max_full_dubins - 1) * dubins_mark_size,
						soluiton_marks_matrix, 1,
						(path_size_full_dubins - 1) * dubins_mark_size,
						matrix_column, solutions_marks_full_Dubins_new);
				matrix_delete_float(solutions_marks_full_Dubins, i_m_c);

				solutions_marks_full_Dubins = matrix_creat_float(i_m_c + 1,
						matrix_column);
				matrix_copy_float(i_m_c + 1, matrix_column,
						solutions_marks_full_Dubins,
						solutions_marks_full_Dubins_new);
				matrix_delete_float(solutions_marks_full_Dubins_new, i_m_c + 1);

				if (path_size_max_full_dubins < path_size_full_dubins)
					path_size_max_full_dubins = path_size_full_dubins;
//							matrix_print_float(states_number*(i_m_c+1), path_size_max_full_dubins, solutions_states_full_Dubins);

			}
			matrix_delete_float(soluiton_states_matrix, states_number);
			matrix_delete_float(soluiton_marks_matrix, 1);

		}
		if (i_m_c != monte_carlo_steps - 1) {
			matrix_delete_float(vertices, states_number);
			matrix_delete_int(edges, edges_row);
			matrix_delete_float(traj, traj_row);
		}
	}
	matrix_delete_float(goal_target_path, states_number);

	if ((goal_mode == 1) && (solution_mode == 2)) {
		vector<float> monte_carlo_avg;

		ofstream results_file;
		results_file.open("results_file_1000000_50_1.txt");

		results_file << "soluitons sizes values are:" << endl;
		results_file << "soluitons matrix for relaxed dubins size is:" << endl;
		results_file << path_size_max_relaxed_dubins << endl;
		results_file << "Dubins markers matrix for relaxed dubins size is:"
				<< endl;
		results_file << (path_size_max_relaxed_dubins - 1) * dubins_mark_size
		<< endl;
		results_file << "soluitons matrix for full dubins size is:" << endl;
		results_file << path_size_max_full_dubins << endl;
		results_file << "Dubins markers matrix for full dubins size is:"
				<< endl;
		results_file << (path_size_max_full_dubins - 1) * dubins_mark_size
		<< endl;
		results_file << endl;
		results_file << endl;

		results_file << "soluitons matrix for relaxed dubins is:" << endl;
		matrix_print_float_to_file(monte_carlo_steps * states_number,
				path_size_max_relaxed_dubins, solutions_states_relaxed_Dubins,
				results_file);
		results_file << endl;

		results_file << "Dubins markers matrix for relaxed dubins is:" << endl;
		matrix_print_float_to_file(monte_carlo_steps,
				(path_size_max_relaxed_dubins - 1) * dubins_mark_size,
				solutions_marks_relaxed_Dubins, results_file);
		results_file << endl;

		results_file << "soluitons matrix for full dubins is:" << endl;
		matrix_print_float_to_file(monte_carlo_steps * states_number,
				path_size_max_full_dubins, solutions_states_full_Dubins,
				results_file);
		results_file << endl;

		results_file << "Dubins markers matrix for full dubins is:" << endl;
		matrix_print_float_to_file(monte_carlo_steps,
				(path_size_max_full_dubins - 1) * dubins_mark_size,
				solutions_marks_full_Dubins, results_file);
		results_file << endl;

		cout << endl;
		cout << "convergence matrix for Relaxed Dubins is:" << endl;
		matrix_print_float(monte_carlo_steps, convergence_size,
				solutions_TTH_Relaxed_Dubins_matrix);
		monte_carlo_avg = matrix_column_avrages(
				solutions_TTH_Relaxed_Dubins_matrix, monte_carlo_steps,
				convergence_size);
		cout << "number of vertices in the tree" << endl;
		vector_print_int(convergence_size, convergence_values);
		cout << "avrages results" << endl;
		vector_print_float(convergence_size, monte_carlo_avg);

		cout << endl;
		cout << endl;
		cout << "convergence matrix for full Dubins is:" << endl;
		matrix_print_float(monte_carlo_steps, convergence_size,
				solutions_TTH_Dubins_matrix);
		monte_carlo_avg = matrix_column_avrages(solutions_TTH_Dubins_matrix,
				monte_carlo_steps, convergence_size);
		cout << "number of vertices in the tree" << endl;
		vector_print_int(convergence_size, convergence_values);
		cout << "avrages results" << endl;
		vector_print_float(convergence_size, monte_carlo_avg);

		results_file << endl;
		results_file << "convergence matrix for Relaxed Dubins is:" << endl;
		matrix_print_float_to_file(monte_carlo_steps, convergence_size,
				solutions_TTH_Relaxed_Dubins_matrix, results_file);
		monte_carlo_avg = matrix_column_avrages(
				solutions_TTH_Relaxed_Dubins_matrix, monte_carlo_steps,
				convergence_size);
		results_file << "number of vertices in the tree" << endl;
		vector_print_int_to_file(convergence_size, convergence_values,
				results_file);
		results_file << "avrages results" << endl;
		vector_print_float_to_file(convergence_size, monte_carlo_avg,
				results_file);

		results_file << endl;
		results_file << endl;
		results_file << "convergence matrix for full Dubins is:" << endl;
		matrix_print_float_to_file(monte_carlo_steps, convergence_size,
				solutions_TTH_Dubins_matrix, results_file);
		monte_carlo_avg = matrix_column_avrages(solutions_TTH_Dubins_matrix,
				monte_carlo_steps, convergence_size);
		results_file << "number of vertices in the tree" << endl;
		vector_print_int_to_file(convergence_size, convergence_values,
				results_file);
		results_file << "avrages results" << endl;
		vector_print_float_to_file(convergence_size, monte_carlo_avg,
				results_file);

		results_file.close();
	}

	cout << "end of sim" << endl;

	//matrix_print_float(states_number, vertices_col, vertices);
	//matrix_print_int(edges_row, edges_col, edges);
//	matrix_print_float(traj_row, traj_col, traj);

	matrix_copy_float(states_number, vertices_col, vertices_final, vertices);
	matrix_copy_int(edges_row, edges_col, edges_final, edges);
	matrix_copy_float(traj_row, traj_col, traj_final, traj);

	return 0;

}

////////
double uniform(double a, double b) {
	return rand() / (RAND_MAX + 1.0) * (b - a) + a;
}

////////
vector<float> sample() {
	vector<float> point;

	for (int i = 0; i < states_number; i++) {
		//uniform_real_distribution<float> distribution(plain[i][0], plain[i][1]);
		//point[i] = distribution(generator);
		point.push_back(uniform(plain[i][0], plain[i][1]));
	}
	return point;
}

////////
int factorial(int n) {
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

////////
float n_sphere_unit_vol() {
	float V;
	int k;
	if (dim % 2 == 0) {
		k = dim / 2;
		V = pow(PI, k) / factorial(k);
	} else {
		k = (dim - 1) / 2;
		V = 2 * factorial(k) * pow(4 * PI, k) / factorial(dim);
	}
	return V;
}

////////
float RRT_star_min_dis(int n) {

	//	return 5;

	float n_f = n;
	float dim_f = dim;
	float eta = ds;
	float V_3d = n_sphere_unit_vol();
	float dim_radius = pow(
			1.5 * space_volume / V_3d * (log(n_f) / n_f) * exp(1 + 1 / dim_f),
			1.0 / dim_f);
	//	cout << "dim_radius:" << dim_radius << ' ' << "eta:" << eta << endl;

	if (eta < dim_radius)
		return eta;

	return dim_radius;
}

////////
int nearest(vector<float> rand_point, float** vertices, int last_step) {

	int vertex_index = 0;
	vector<float> temp_vector;
	temp_vector = export_vector_from_matrix(vertices, states_number,
			vertex_index);
	float distance = dis_state_number(temp_vector, rand_point);
	float new_distance;
	for (int i = 1; i < last_step; i++) {
		temp_vector = export_vector_from_matrix(vertices, states_number, i);
		new_distance = dis_state_number(temp_vector, rand_point);
		if ((new_distance < distance) && (temp_vector[3] < rand_point[3])) {
			distance = new_distance;
			vertex_index = i;
		}
	}
	return vertex_index;
}

////////
vector<float> linear_interpulation(vector<float> start_point,
		vector<float> end_point, float wanted_distance) {
	float dis_full = dis3(start_point, end_point);
	//	if (dis_full < wanted_distance)
	//	{
	//		return end_point;
	//	}

	vector<float> new_point;

	for (int i = 0; i < states_number; i++) {
		new_point.push_back(
				start_point[i]
						+ (wanted_distance / dis_full)
								* (end_point[i] - start_point[i]));
	}
	return new_point;
}
////////
vector<float> linear_interpulation_2d(vector<float> start_point,
		vector<float> end_point, float wanted_distance) {
	float dis_full = dis2(start_point, end_point);
	//	if (dis_full < wanted_distance)
	//	{
	//		return end_point;
	//	}

	vector<float> new_point;

	for (int i = 0; i < 2; i++) {
		new_point.push_back(
				start_point[i]
						+ (wanted_distance / dis_full)
								* (end_point[i] - start_point[i]));
	}
	return new_point;
}

////////
float S_length_from_path_marks(vector<float> path_marks) {
	// path_marks_and_vector= {	R, c, L,																	--- 0-2
	//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
	//							v_norm,																		--- 7
	//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
	//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

	float s;

	float R = path_marks[0];

	float factor_S = sqrt(1 + 1 / pow(V_max, 2));
	float factor_C = sqrt(1 + 1 / pow(V_max, 2) + 1 / pow(R, 2));

	float L_tot = path_marks[2];
	float initialCircle_centre1 = path_marks[3];
	float initialCircle_centre2 = path_marks[4];
	float initialArcAngle = path_marks[5];
	float phi_pinit = path_marks[6];
	float L_first_circle = initialArcAngle * R;

	if (L_tot <= L_first_circle) {
		s = L_tot * factor_C;
		return s;
	}

	float L_straight = path_marks[7];

	if ((L_tot > L_first_circle) && (L_tot < L_first_circle + L_straight)) {
		s = L_first_circle * factor_C + (L_tot - L_first_circle) * factor_S;
		return s;
	}

	float finalCircle_centre1 = path_marks[8];
	float finalCircle_centre2 = path_marks[9];
	float finalArcAngle = path_marks[10];
	float phi_contact = path_marks[11];
	float L_last_circle = finalArcAngle * R;

	if ((L_tot > L_first_circle) && (L_tot > L_first_circle + L_straight)
			&& (L_tot < L_first_circle + L_straight + L_last_circle)) {
		s = L_first_circle * factor_C + L_straight * factor_S
				+ (L_tot - L_first_circle - L_straight) * factor_C;
		return s;
	} else
		s = L_first_circle * factor_C + L_straight * factor_S
				+ L_last_circle * factor_C;

	return s;
}

////////
vector<float> path_Dubins_extend(vector<float> start_point,
		vector<float> end_point, vector<float> path_marks, float ds_new) {

	// copy last known path_marks_and_vector
	vector<float> path_marks_and_vector;

	for (int i = 0; i < dubins_mark_size; i++) {
		path_marks_and_vector.push_back(path_marks[i]);
	}
	for (int i = 0; i < states_number; i++) {
		path_marks_and_vector.push_back(end_point[i]);
	}

	// path_marks_and_vector= {	R, c, L,																	--- 0-2
	//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
	//							v_norm,																		--- 7
	//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
	//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

	// calculate at confirm extension size
	float ds_last = S_length_from_path_marks(path_marks);
	if (ds_last >= ds_new)
		return path_marks_and_vector;

	float s_factor = ds_new / ds_last;

	// Parameter notations
	float R = path_marks[0];
	float c = path_marks[1];
	float L_tot = path_marks[2];
	float initialCircle_centre1 = path_marks[3];
	float initialCircle_centre2 = path_marks[4];
	float initialArcAngle = path_marks[5];
	float phi_pinit = path_marks[6];
	float L_first_circle = initialArcAngle * R;
	float L_straight = path_marks[7];
	float finalCircle_centre1 = path_marks[8];
	float finalCircle_centre2 = path_marks[9];
	float finalArcAngle = path_marks[10];
	float phi_contact = path_marks[11];
	float L_last_circle = finalArcAngle * R;

	// increase unnecessary values
	L_tot = L_tot * s_factor;
	initialArcAngle = initialArcAngle * s_factor;
	L_first_circle = L_first_circle * s_factor;
	L_straight = L_straight * s_factor;
	finalArcAngle = finalArcAngle * s_factor;
	L_last_circle = L_last_circle * s_factor;

	// calculate new straight line end point
	vector<float> line_start, line_end;

	if ((c == 1) || (c == 3) || (c == 11) || (c == 13) || (c == 21)
			|| (c == 23)) // start circle is L
			{
		line_start.push_back(
				initialCircle_centre1 + R * cos(initialArcAngle + phi_pinit));
		line_start.push_back(
				initialCircle_centre2 + R * sin(initialArcAngle + phi_pinit));
		line_start.push_back(phi_pinit + PI / 2 + initialArcAngle);
	} else {
		line_start.push_back(
				initialCircle_centre1 + R * cos(initialArcAngle - phi_pinit));
		line_start.push_back(
				initialCircle_centre2 - R * sin(initialArcAngle - phi_pinit));
		line_start.push_back(phi_pinit - PI / 2 - initialArcAngle);
	}

	line_end.push_back(line_start[0] + L_straight * cos(line_start[2]));
	line_end.push_back(line_start[1] + L_straight * sin(line_start[2]));
	line_end.push_back(line_start[2]);

	// reconstruct final circle
	if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
			|| (c == 24)) // end circle is L
			{
		phi_contact = line_end[2] - PI / 2;
		finalCircle_centre1 = line_end[0] - R * cos(phi_contact);
		finalCircle_centre2 = line_end[1] - R * sin(phi_contact);
	} else {
		phi_contact = line_end[2] + PI / 2;
		finalCircle_centre1 = line_end[0] - R * cos(phi_contact);
		finalCircle_centre2 = line_end[1] + R * sin(phi_contact);
	}

	// replace path_marks values
	path_marks_and_vector[2] = L_tot;
	path_marks_and_vector[5] = initialArcAngle;
	path_marks_and_vector[7] = L_straight;
	path_marks_and_vector[8] = finalCircle_centre1;
	path_marks_and_vector[9] = finalCircle_centre2;
	path_marks_and_vector[10] = finalArcAngle;
	path_marks_and_vector[11] = phi_contact;

	// calculate new path end point
	if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
			|| (c == 24)) // path is either LSL or RSL
			{
		path_marks_and_vector[12] = finalCircle_centre1
				+ R * cos(L_last_circle / R + phi_contact);
		path_marks_and_vector[13] = finalCircle_centre2
				+ R * sin(L_last_circle / R + phi_contact);
		path_marks_and_vector[14] = phi_contact + PI / 2 + L_last_circle / R;
		path_marks_and_vector[15] = start_point[3] + L_tot / V_max;
	} else {
		path_marks_and_vector[12] = finalCircle_centre1
				+ R * cos(L_last_circle / R - phi_contact);
		path_marks_and_vector[13] = finalCircle_centre2
				- R * sin(L_last_circle / R - phi_contact);
		path_marks_and_vector[14] = phi_contact - PI / 2 - L_last_circle / R;
		path_marks_and_vector[15] = start_point[3] + L_tot / V_max;
	}

	return path_marks_and_vector;
}

////////
vector<float> path_Dubins_cut(vector<float> start_point,
		vector<float> end_point, vector<float> path_mark_vector, float ds_new) {
	vector<float> path_marks_and_vector;

	// path_marks_and_vector= {	R, c, L,																	--- 0-2
	//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
	//							v_norm,																		--- 7
	//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
	//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

	float R = path_mark_vector[0];
	float c = path_mark_vector[1];
	float L_tot = path_mark_vector[2];
	float initialCircle_centre1 = path_mark_vector[3];
	float initialCircle_centre2 = path_mark_vector[4];
	float initialArcAngle = path_mark_vector[5];
	float phi_pinit = path_mark_vector[6];
	float L_first_circle = initialArcAngle * R;
	float L_straight = path_mark_vector[7];
	float finalCircle_centre1 = path_mark_vector[8];
	float finalCircle_centre2 = path_mark_vector[9];
	float finalArcAngle = path_mark_vector[10];
	float phi_contact = path_mark_vector[11];
	float L_last_circle = finalArcAngle * R;

	float factor_S = sqrt(1 + 1 / pow(V_max, 2));
	float factor_C = sqrt(1 + 1 / pow(V_max, 2) + 1 / pow(R, 2));

	float s_first_circle = L_first_circle * factor_C;
	float s_straight = L_straight * factor_S;
	float s_last_circle = L_last_circle * factor_C;

//	check= s_first_circle+s_straight+s_last_circle;

	float dL_new;
	if (ds_new <= s_first_circle)
		dL_new = ds_new / factor_C;
	if ((ds_new > s_first_circle) && (ds_new < s_first_circle + s_straight))
		dL_new = s_first_circle / factor_C
				+ (ds_new - s_first_circle) / factor_S;
	if ((ds_new > s_first_circle + s_straight)
			&& (ds_new < s_first_circle + s_straight + s_last_circle))
		dL_new = s_first_circle / factor_C + s_straight / factor_S
				+ (ds_new - s_first_circle - s_straight) / factor_C;
	if (ds_new > s_first_circle + s_straight + s_last_circle)
		dL_new = s_first_circle / factor_C + s_straight / factor_S
				+ s_last_circle / factor_C;

	path_marks_and_vector.push_back(path_mark_vector[0]);
	path_marks_and_vector.push_back(path_mark_vector[1]);
	path_marks_and_vector.push_back(dL_new);
	path_marks_and_vector.push_back(path_mark_vector[3]);
	path_marks_and_vector.push_back(path_mark_vector[4]);
	path_marks_and_vector.push_back(path_mark_vector[5]);
	path_marks_and_vector.push_back(path_mark_vector[6]);

	///  path cut is on initial circle
	if (L_first_circle >= dL_new) {
		for (int j = 7; j < 12; j++) {
			path_marks_and_vector.push_back(-999);
		}
		if ((c == 1) || (c == 3) || (c == 11) || (c == 13) || (c == 21)
				|| (c == 23)) // path is either LSL or LSR
				{
			path_marks_and_vector.push_back(
					initialCircle_centre1 + R * cos(dL_new / R + phi_pinit));
			path_marks_and_vector.push_back(
					initialCircle_centre2 + R * sin(dL_new / R + phi_pinit));
			path_marks_and_vector.push_back(phi_pinit + PI / 2 + dL_new / R);
			path_marks_and_vector.push_back(start_point[3] + dL_new / V_max);
		} else // path is either RSR or RSL
		{
			path_marks_and_vector.push_back(
					initialCircle_centre1 + R * cos(dL_new / R - phi_pinit));
			path_marks_and_vector.push_back(
					initialCircle_centre2 - R * sin(dL_new / R - phi_pinit));
			path_marks_and_vector.push_back(phi_pinit - PI / 2 - dL_new / R);
			path_marks_and_vector.push_back(start_point[3] + dL_new / V_max);
		}

		return path_marks_and_vector;
	}

	path_marks_and_vector.push_back(path_mark_vector[7]);
	path_marks_and_vector.push_back(path_mark_vector[8]);
	path_marks_and_vector.push_back(path_mark_vector[9]);
	path_marks_and_vector.push_back(path_mark_vector[10]);
	path_marks_and_vector.push_back(path_mark_vector[11]);

	///  path cut is on the straight line
	if (L_first_circle + L_straight >= dL_new) {
		path_marks_and_vector[10] = -999;

		vector<float> line_start;
		vector<float> line_end;

		if ((c == 1) || (c == 3) || (c == 11) || (c == 13) || (c == 21)
				|| (c == 23)) // start circle is L
				{
			line_start.push_back(
					initialCircle_centre1
							+ R * cos(initialArcAngle + phi_pinit));
			line_start.push_back(
					initialCircle_centre2
							+ R * sin(initialArcAngle + phi_pinit));
			line_start.push_back(phi_pinit + PI / 2 + initialArcAngle);
		} else {
			line_start.push_back(
					initialCircle_centre1
							+ R * cos(initialArcAngle - phi_pinit));
			line_start.push_back(
					initialCircle_centre2
							- R * sin(initialArcAngle - phi_pinit));
			line_start.push_back(phi_pinit - PI / 2 - initialArcAngle);
		}

		if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
				|| (c == 24)) // end circle is L
				{
			line_end.push_back(finalCircle_centre1 + R * cos(phi_contact));
			line_end.push_back(finalCircle_centre2 + R * sin(phi_contact));
			line_end.push_back(phi_contact + PI / 2);
		} else {
			line_end.push_back(finalCircle_centre1 + R * cos(-phi_contact));
			line_end.push_back(finalCircle_centre2 - R * sin(-phi_contact));
			line_end.push_back(phi_contact - PI / 2);
		}

		vector<float> line_cut = linear_interpulation_2d(line_start, line_end,
				dL_new - L_first_circle);

		path_marks_and_vector.push_back(line_cut[0]);
		path_marks_and_vector.push_back(line_cut[1]);
		path_marks_and_vector.push_back(line_start[2]);
		path_marks_and_vector.push_back(start_point[3] + dL_new / V_max);
		return path_marks_and_vector;
	}

	///  path cut is on final circle

	if (L_first_circle + L_straight + L_last_circle >= dL_new) {
		float dL_last_circle = dL_new - L_first_circle - L_straight;
		if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
				|| (c == 24)) // path is either LSL or RSL
				{
			path_marks_and_vector.push_back(
					finalCircle_centre1
							+ R * cos(dL_last_circle / R + phi_contact));
			path_marks_and_vector.push_back(
					finalCircle_centre2
							+ R * sin(dL_last_circle / R + phi_contact));
			path_marks_and_vector.push_back(
					phi_contact + PI / 2 + dL_last_circle / R);
			path_marks_and_vector.push_back(start_point[3] + dL_new / V_max);
		} else {
			path_marks_and_vector.push_back(
					finalCircle_centre1
							+ R * cos(dL_last_circle / R - phi_contact));
			path_marks_and_vector.push_back(
					finalCircle_centre2
							- R * sin(dL_last_circle / R - phi_contact));
			path_marks_and_vector.push_back(
					phi_contact - PI / 2 - dL_last_circle / R);
			path_marks_and_vector.push_back(start_point[3] + dL_new / V_max);
		}
		return path_marks_and_vector;
	}

	///  path extention

	if (L_tot < dL_new) {
		float dL_addition = dL_new - L_tot;

		if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
				|| (c == 24)) // path is either LSL or RSL
				{
			path_marks_and_vector.push_back(
					finalCircle_centre1
							+ R * cos(dL_addition / R + phi_contact));
			path_marks_and_vector.push_back(
					finalCircle_centre2
							+ R * sin(dL_addition / R + phi_contact));
			path_marks_and_vector.push_back(
					phi_contact + PI / 2 + dL_addition / R);
			path_marks_and_vector.push_back(start_point[3] + dL_new / V_max);
		} else {
			path_marks_and_vector.push_back(
					finalCircle_centre1
							+ R * cos(dL_addition / R - phi_contact));
			path_marks_and_vector.push_back(
					finalCircle_centre2
							- R * sin(dL_addition / R - phi_contact));
			path_marks_and_vector.push_back(
					phi_contact - PI / 2 - dL_addition / R);
			path_marks_and_vector.push_back(start_point[3] + dL_new / V_max);
		}
		return path_marks_and_vector;
	}
}

////////
bool is_obstacle(vector<float> point) {
	if (obst_num == 0) {
		return false;
	}

	bool bool_var = false;
	bool obst_check = true;

	float time = point[3];
	for (int i = 1; i < obst_num + 1; i++) {
		for (int j = 0; j < states_number - 1; j++) {
			if ((point[j]
					< obstacles[(i - 1) * (states_number - 1) + j][0]
							+ time
									* obstacles[(i - 1) * (states_number - 1)
											+ j][2])
					|| (point[j]
							> obstacles[(i - 1) * (states_number - 1) + j][1]
									+ time
											* obstacles[(i - 1)
													* (states_number - 1) + j][2]))
//			if ((point[j]<obstacles[(i - 1)*(states_number-1) + j][0]) ||
//				(point[j]>obstacles[(i - 1)*(states_number-1) + j][1]))
					{
				obst_check = false;
			}
		}
		if (obst_check) {
			bool_var = true;
		}
		obst_check = true;
	}
	return bool_var;
}

////////
bool in_plain(vector<float> point) {

	bool bool_var = true;

	for (int j = 0; j < states_number; j++) {
		if ((point[j] < plain[j][0]) || (point[j] > plain[j][1])) {
			bool_var = false;
		}
	}

	return bool_var;
}

////////
bool is_vectors_equal(vector<float> a, vector<float> b) {
	bool bool_var = true;

	for (int i = 0; i < states_number; i++) {
		if (a[i] != b[i]) {
			bool_var = false;
		}
	}

	return bool_var;
}

////////
bool path_plain_obstacle_check(vector<float> path_marks,
		vector<float> start_point) {
	// path_marks_and_vector= {	R, c, L,																	--- 0-2
	//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
	//							v_norm,																		--- 7
	//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
	//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

	vector<float> new_vector;
	float R = path_marks[0];
	float factor_S = sqrt(1 + 1 / pow(V_max, 2));
	float factor_C = sqrt(1 + 1 / pow(V_max, 2) + 1 / pow(R, 2));

	float dL_obs_C = ds_obs / factor_C;
	float dL_obs_S = ds_obs / factor_S;
	float dL_check = dL_obs_C;
	float c = path_marks[1];
	float path_full_length = path_marks[2];
	float initialCircle_centre1 = path_marks[3];
	float initialCircle_centre2 = path_marks[4];
	float initialArcAngle = path_marks[5];
	float phi_pinit = path_marks[6];
	float L_first_circle = initialArcAngle * R;

	new_vector.push_back(-999);
	new_vector.push_back(-999);
	new_vector.push_back(-999);
	new_vector.push_back(-999);

	//check initial circle
	while ((dL_check < path_full_length) && (dL_check < L_first_circle)) {

		if ((c == 1) || (c == 3) || (c == 11) || (c == 13) || (c == 21)
				|| (c == 23)) // path is either LSL or LSR
				{
			new_vector[0] = initialCircle_centre1
					+ R * cos(dL_check / R + phi_pinit);
			new_vector[1] = initialCircle_centre2
					+ R * sin(dL_check / R + phi_pinit);
			new_vector[2] = phi_pinit + PI / 2 + dL_check / R;
			new_vector[3] = start_point[3] + dL_check / V_max;
		} else // path is either RSR or RSL
		{
			new_vector[0] = initialCircle_centre1
					+ R * cos(dL_check / R - phi_pinit);
			new_vector[1] = initialCircle_centre2
					- R * sin(dL_check / R - phi_pinit);
			new_vector[2] = phi_pinit - PI / 2 - dL_check / R;
			new_vector[3] = start_point[3] + dL_check / V_max;
		}

		new_vector[2] = angle_to_pi_interval(new_vector[2]);
		if (is_obstacle(new_vector) || !in_plain(new_vector)) {
			return false;
		}
		dL_check = dL_check + dL_obs_C;
	}

	if (path_marks[7] == -999)
		return true;

	//check straight line
	float L_straight = path_marks[7];
	float finalCircle_centre1 = path_marks[8];
	float finalCircle_centre2 = path_marks[9];
	float finalArcAngle = path_marks[10];
	float phi_contact = path_marks[11];

	vector<float> line_start;
	vector<float> line_end;
	if ((c == 1) || (c == 3) || (c == 11) || (c == 13) || (c == 21)
			|| (c == 23)) // start circle is L
			{
		line_start.push_back(
				initialCircle_centre1 + R * cos(initialArcAngle + phi_pinit));
		line_start.push_back(
				initialCircle_centre2 + R * sin(initialArcAngle + phi_pinit));
		line_start.push_back(phi_pinit + PI / 2 + initialArcAngle);
	} else {
		line_start.push_back(
				initialCircle_centre1 + R * cos(initialArcAngle - phi_pinit));
		line_start.push_back(
				initialCircle_centre2 - R * sin(initialArcAngle - phi_pinit));
		line_start.push_back(phi_pinit - PI / 2 - initialArcAngle);
	}

	if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
			|| (c == 24)) // end circle is L
			{
		line_end.push_back(finalCircle_centre1 + R * cos(phi_contact));
		line_end.push_back(finalCircle_centre2 + R * sin(phi_contact));
		line_end.push_back(phi_contact + PI / 2);
	} else {
		line_end.push_back(finalCircle_centre1 + R * cos(-phi_contact));
		line_end.push_back(finalCircle_centre2 - R * sin(-phi_contact));
		line_end.push_back(phi_contact - PI / 2);
	}

	while ((dL_check < path_full_length)
			&& (dL_check < L_first_circle + L_straight)) {

		vector<float> line_cut = linear_interpulation_2d(line_start, line_end,
				dL_check - L_first_circle);

		new_vector[0] = line_cut[0];
		new_vector[1] = line_cut[1];
		new_vector[2] = line_start[2];
		new_vector[3] = start_point[3] + dL_check / V_max;

		new_vector[2] = angle_to_pi_interval(new_vector[2]);
		if (is_obstacle(new_vector) || !in_plain(new_vector)) {
			return false;
		}
		dL_check = dL_check + dL_obs_S;
	}

	if (path_marks[10] == -999)
		return true;

	//check final circle
	float L_last_circle = finalArcAngle * R;

	while ((dL_check < path_full_length)
			&& (dL_check < L_first_circle + L_straight + L_last_circle)) {

		float dL_last_circle = dL_check - L_first_circle - L_straight;
		if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
				|| (c == 24)) // path is either LSL or RSL
				{
			new_vector[0] = finalCircle_centre1
					+ R * cos(dL_last_circle / R + phi_contact);
			new_vector[1] = finalCircle_centre2
					+ R * sin(dL_last_circle / R + phi_contact);
			new_vector[2] = phi_contact + PI / 2 + dL_last_circle / R;
			new_vector[3] = start_point[3] + dL_check / V_max;
		} else {
			new_vector[0] = finalCircle_centre1
					+ R * cos(dL_last_circle / R - phi_contact);
			new_vector[1] = finalCircle_centre2
					- R * sin(dL_last_circle / R - phi_contact);
			new_vector[2] = phi_contact - PI / 2 - dL_last_circle / R;
			new_vector[3] = start_point[3] + dL_check / V_max;
		}

		new_vector[2] = angle_to_pi_interval(new_vector[2]);
		if (is_obstacle(new_vector) || !in_plain(new_vector)) {
			return false;
		}
		dL_check = dL_check + dL_obs_C;
	}

	return true;
}

////////
vector<float> steer(vector<float> vertex_start, vector<float> rand_point) {
	// path_marks_and_vector= {	R, c, L,																	--- 0-2
	//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
	//							v_norm,																		--- 7
	//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
	//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

	vector<float> vertex_new_and_path_marks;
	float path_t = rand_point[3] - vertex_start[3];
	float path_L_min = dis2(vertex_start, rand_point);
	float path_V_min = path_L_min / path_t;

	if (path_t <= 0) {
		vertex_new_and_path_marks.push_back(-999);
		return vertex_new_and_path_marks;
	}

	vector<float> path_marks = Dubins_state_time_steer(vertex_start,
			rand_point);
//					cout << path_marks[0] << endl;
	if (path_marks[0] == -999) {
		vertex_new_and_path_marks.push_back(-999);
		return vertex_new_and_path_marks;
	}

	float s_path = S_length_from_path_marks(path_marks);
//	cout << s_path << endl;

	float ds_new;
	if (steer_length_mode == 1) {
		if (s_path <= ds) {
			if (is_obstacle(rand_point) || !in_plain(rand_point)) {
				vertex_new_and_path_marks.push_back(-999);
				return vertex_new_and_path_marks;
			}

			for (int i = 0; i < path_marks.size(); i++) {
				vertex_new_and_path_marks.push_back(path_marks[i]);
			}
			for (int i = 0; i < rand_point.size(); i++) {
				vertex_new_and_path_marks.push_back(rand_point[i]);
			}
		} else {
			vertex_new_and_path_marks = path_Dubins_cut(vertex_start,
					rand_point, path_marks, ds);
		}
	}

	if (steer_length_mode == 3) {
		ds_new = uniform(0, ds);

		if (ds_new < s_path) {
			vertex_new_and_path_marks = path_Dubins_cut(vertex_start,
					rand_point, path_marks, ds_new);
		} else {
			vertex_new_and_path_marks = path_Dubins_extend(vertex_start,
					rand_point, path_marks, ds_new);
		}
	}

	if (steer_length_mode == 2) {
		ds_new = ds;
		if (ds_new < s_path) {
			vertex_new_and_path_marks = path_Dubins_cut(vertex_start,
					rand_point, path_marks, ds_new);
		} else {
			vertex_new_and_path_marks = path_Dubins_extend(vertex_start,
					rand_point, path_marks, ds_new);
		}
	}

	vector<float> vertex_new;
	for (int i = 0; i < states_number; i++) {
		vertex_new.push_back(vertex_new_and_path_marks[dubins_mark_size + i]);
	}
	//changing angle to reside between -PI to PI
	vertex_new[2] = angle_to_pi_interval(vertex_new[2]);
	vertex_new_and_path_marks[dubins_mark_size + 2] = angle_to_pi_interval(
			vertex_new[2]);

//	vector_print_float(states_number,vertex_new);

	//	cout << is_obstacle(vertex_new)  << ' ' << !in_plain(vertex_new) << endl;
	if (is_obstacle(vertex_new) || !in_plain(vertex_new)) {

		for (int i = 0; i < states_number; i++) {
			vertex_new[i] = -999;
		}
		return vertex_new;
	}

	for (int i = 0; i < dubins_mark_size; i++) {
		path_marks[i] = vertex_new_and_path_marks[i];
	}

	if (!path_plain_obstacle_check(path_marks, vertex_new)) {
		vertex_new_and_path_marks[0] = -999;
		return vertex_new_and_path_marks;
	}

	return vertex_new_and_path_marks;
}

////////
vector<float> forced_steer(vector<float> vertex_start,
		vector<float> vertex_end) {
	// path_marks_and_vector= {	R, c, L,																	--- 0-2
	//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
	//							v_norm,																		--- 7
	//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
	//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

	vector<float> path_marks;
	float path_t = vertex_end[3] - vertex_start[3];
	float path_L_min = dis2(vertex_start, vertex_end);
	float path_V_min = path_L_min / path_t;

	if ((path_t <= 0) || (path_V_min > V_max)) {
		path_marks.push_back(-999);
		return path_marks;
	}

	path_marks = Dubins_state_time_steer(vertex_start, vertex_end);
	if (path_marks[0] == -999) {
		return path_marks;
	}

	//	cout << is_obstacle(vertex_new)  << ' ' << !in_plain(vertex_new) << endl;

	if (!path_plain_obstacle_check(path_marks, vertex_start)) {
		path_marks[0] = -999;
		return path_marks;
	}

	return path_marks;
}

////////
void RRT_function(int steps, float** vertices, int** edges, float** traj) {
	int i = 0;
	//// kdtree code ////
	struct kdtree *kd_vertices = kd_create(states_number);

	vector<float> rand_point, vertex_new, vertex_new_and_path_marks,
			dubins_mark;
	for (int vertex_mone = 0; vertex_mone < states_number; vertex_mone++) {
		vertex_new.push_back(-999);
	}
	for (int vertex_mone = 0; vertex_mone < dubins_mark_size; vertex_mone++) {
		dubins_mark.push_back(-999);
	}
	vector<int> new_edge;
	for (int edge_mone = 0; edge_mone < edges_row; edge_mone++) {
		new_edge.push_back(-999);
	}
	vector<float> new_traj;
	for (int traj_mone = 0; traj_mone < traj_row; traj_mone++)
		new_traj.push_back(-999);
	int vertex_index;

	// kdtree code ////
	////////
	float* vertex_temp = invert_vector_to_array_float(
			export_vector_from_matrix(vertices, states_number, 0),
			states_number);
	//	cout << vertex_temp[0] << ' ' << vertex_temp[1] << ' ' << vertex_temp[2] << ' ' << vertex_temp[3] << endl;
	//	cout << &vertex_temp << endl;
	kd_insertf(kd_vertices, vertex_temp, &firts_point);
	vertex_temp[2] = vertex_temp[2] + 2 * PI;
	kd_insertf(kd_vertices, vertex_temp, &firts_point);
	vertex_temp[2] = vertex_temp[2] - 4 * PI;
	kd_insertf(kd_vertices, vertex_temp, &firts_point);
	delete vertex_temp;
	double pos_n[states_number];
	int *vertex_index_pointer;
	//	cout << vertex_temp[0] << ' ' << vertex_temp[1] << ' ' << vertex_temp[2] << ' ' << endl;
	//	cout << &vertex_temp << endl;

	while (i < steps) {

		i = i + 1;
		rand_point = sample();
//		rand_point[0] =1.5;
//		rand_point[1] =1.7;
//		rand_point[2] =60*d2r;
//		rand_point[3] =0.08625;

		//// kdtree code ////

//		if (i==15)
//			cout << i << endl;
//		cout << "random sample: ";
//		cout << rand_point[0] << ' ' << rand_point[1] << ' ' << rand_point[2] << ' ' << rand_point[3] << endl;
		if (kdtree_mode == 2) {

			struct kdres *kd_result;
			float* vertex_temp = invert_vector_to_array_float(rand_point,
					states_number);
			kd_result = kd_nearestf(kd_vertices, vertex_temp);
			delete vertex_temp;
			vertex_index_pointer = (int*) kd_res_item(kd_result, pos_n);
			vertex_index = *vertex_index_pointer;
			kd_res_free(kd_result);

//						int vertex_index2 = nearest(rand_point, vertices, i);
//						cout << "vertex index from nearest is: " << vertex_index << endl;
//						cout << "pos is: " << pos_n[0] << ' '  << pos_n[1] << ' ' << pos_n[2] << ' ' << pos_n[3] <<  endl;

		} else {
			vertex_index = nearest(rand_point, vertices, i);
		}

		vector<float> vertex_last = export_vector_from_matrix(vertices,
				states_number, vertex_index);
//			if ((i == 15) && (fabs(vertex_last[0] - 4.44384) < 0.0001)) {
//				cout << "random sample: ";
//				cout << rand_point[0] << ' ' << rand_point[1] << ' '
//						<< rand_point[2] << ' ' << rand_point[3] << endl;
//
//				cout << "vertex last is: ";
//				cout << vertex_last[0] << ' ' << vertex_last[1] << ' '
//						<< vertex_last[2] << ' ' << vertex_last[3] << endl;
//			}
		vertex_new_and_path_marks = steer(vertex_last, rand_point);

		if (vertex_new_and_path_marks[0] == -999)
			i = i - 1;
		else {
			for (int j = 0; j < states_number; j++) {
				vertex_new[j] = vertex_new_and_path_marks[dubins_mark_size + j];
			}
			for (int j = 0; j < dubins_mark_size; j++) {
				dubins_mark[j] = vertex_new_and_path_marks[j];
			}
//				if ((i == 15) && (fabs(vertex_last[0] - 4.44384) < 0.0001)) {
//					cout << "vertex new is: ";
//					cout << vertex_new[0] << ' ' << vertex_new[1] << ' '
//							<< vertex_new[2] << ' ' << vertex_new[3] << endl;
//					vector_print_float(dubins_mark_size, dubins_mark);
//				}
			new_edge[0] = vertex_index;
			new_edge[1] = i;

			new_traj[0] = cost_func(dubins_mark);
			for (int j = 0; j < dubins_mark_size; j++) {
				new_traj[j + 1] = dubins_mark[j];
			}

			insert_vector_to_matrix_float(vertices, states_number, vertex_new,
					i);
			insert_vector_to_matrix_int(edges, edges_row, new_edge, i - 1);
			insert_vector_to_matrix_float(traj, traj_row, new_traj, i - 1);

			//// kdtree code ////
			if (kdtree_mode == 2) {
				float* vertex_temp = invert_vector_to_array_float(vertex_new,
						states_number);
				kd_insertf(kd_vertices, vertex_temp, &edges[1][i - 1]);
				vertex_temp[2] = vertex_temp[2] + 2 * PI;
				kd_insertf(kd_vertices, vertex_temp, &edges[1][i - 1]);
				vertex_temp[2] = vertex_temp[2] - 4 * PI;
				kd_insertf(kd_vertices, vertex_temp, &edges[1][i - 1]);
				delete vertex_temp;
			}

			if ((i % 10000) == 0) {
				cout << i << endl;
			}
		}

		//		matrix_print_int(edges_row, i, edges);

	}
//		matrix_print_float(states_number, steps+1, vertices);
//		matrix_print_int(edges_row, steps, edges);
//		matrix_print_float(traj_row, steps, traj);

	//// kdtree code ////
	kd_free(kd_vertices);
}

////////
void RRT_star_function(int steps, float** vertices, int** edges, float** traj) {
	int i = 0;

	struct kdtree *kd_vertices = kd_create(4);
	struct kdres *vertex_near_set;

	vector<float> rand_point, vertex_new, vertex_new_and_path_marks,
			new_dubins_mark, best_dubins_mark, temp_dubins_mark;
	vector<float> vertex_last;
	int vertex_index;
	float** path;
	vector<int> new_edge;
	for (int edge_mone = 0; edge_mone < edges_row; edge_mone++)
		new_edge.push_back(-999);
	vector<float> new_traj;
	for (int traj_mone = 0; traj_mone < traj_row; traj_mone++)
		new_traj.push_back(-999);

	for (int j = 0; j < states_number; j++) {
		vertex_new.push_back(-999);
	}
	for (int j = 0; j < dubins_mark_size; j++) {
		new_dubins_mark.push_back(-999);
		best_dubins_mark.push_back(-999);
	}
	//// kdtree code ////

	vector<float> initial_vector = export_vector_from_matrix(vertices,
			states_number, 0);
	float* initial_vector2 = invert_vector_to_array_float(initial_vector,
			states_number);
	kd_insertf(kd_vertices, initial_vector2, &firts_point);
	initial_vector2[2] = initial_vector2[2] + 2 * PI;
	kd_insertf(kd_vertices, initial_vector2, &firts_point);
	initial_vector2[2] = initial_vector2[2] - 4 * PI;
	kd_insertf(kd_vertices, initial_vector2, &firts_point);
	delete initial_vector2;

	double pos_n[states_number];
	int *vertex_index_pointer;
	float new_cost, near_radius, best_new_vertex_cost, new_vertex_near_cost;
	int* temp_vertex_near_index_p;

	while (i < steps) {
		i = i + 1;
		rand_point = sample();
//				cout << "random sample: ";
//				cout << rand_point[0] << ' ' << rand_point[1] << ' ' << rand_point[2] << ' ' << endl;

		//// kdtree code ////

//		if (i==13)
//			cout << i << endl;
//						cout << "random sample: ";
//						cout << rand_point[0] << ' ' << rand_point[1] << ' ' << rand_point[2] << ' '  << rand_point[3] << endl;

		if (kdtree_mode == 2) {
			struct kdres *kd_result;
			float* vertex_temp = invert_vector_to_array_float(rand_point,
					states_number);
			kd_result = kd_nearestf(kd_vertices, vertex_temp);
			delete vertex_temp;
			vertex_index_pointer = (int*) kd_res_item(kd_result, pos_n);
			vertex_index = *vertex_index_pointer;
			kd_res_free(kd_result);

			//			int vertex_index2 = nearest(rand_point, vertices, i);
//						cout << "vertex index from nearest is: " << vertex_index << endl;
		} else {
			vertex_index = nearest(rand_point, vertices, i);
		}

		vertex_last = export_vector_from_matrix(vertices, states_number,
				vertex_index);
//		cout << "vertex last is: ";
//		cout << vertex_last[0] << ' ' << vertex_last[1] << ' ' << vertex_last[2] << ' '  << vertex_last[3] << endl;

		vertex_new_and_path_marks = steer(vertex_last, rand_point);
//				cout << "vertex new is: ";
//				cout << vertex_new_and_path_marks[12] << ' ' << vertex_new_and_path_marks[13] << ' ' << vertex_new_and_path_marks[14] << ' ' << vertex_new_and_path_marks[15] << endl;

		if (vertex_new_and_path_marks[0] != -999) // vertex_new[0]==-999 means that the new path treverses an obstacle or steps outside the plain
				{
			int temp_vertex_near_index, best_vertex_near_index;
			vector<float> temp_vertex_near;
			float pos[states_number], dist, temp_vertex_near_cost,
					temp_to_new_cost, new_to_temp_cost, best_to_new_cost;

			for (int j = 0; j < states_number; j++) {
				vertex_new[j] = vertex_new_and_path_marks[dubins_mark_size + j];
			}
			for (int j = 0; j < dubins_mark_size; j++) {
				new_dubins_mark[j] = vertex_new_and_path_marks[j];
			}

//			cout << "vertex new is: ";
//			cout << vertex_new[0] << ' ' << vertex_new[1] << ' ' << vertex_new[2] << ' ' << vertex_new[3] << endl;
//			vector_print_float(dubins_mark_size, new_dubins_mark);

			new_cost = cost_func(new_dubins_mark);

//						cout << "new cost: " << new_cost << endl;
			new_edge[0] = vertex_index;
			new_edge[1] = i;
			new_traj[0] = new_cost;
			for (int j = 0; j < dubins_mark_size; j++) {
				new_traj[j + 1] = new_dubins_mark[j];
			}

			insert_vector_to_matrix_float(vertices, states_number, vertex_new,
					i);
			insert_vector_to_matrix_int(edges, edges_row, new_edge, i - 1);
			insert_vector_to_matrix_float(traj, traj_row, new_traj, i - 1);
//						cout << traj[0][0] << ' ' << traj[0][1] << ' ' << traj[0][2] << ' ' << endl;

			best_new_vertex_cost = vertex_cost_value(i, edges, traj);
//						cout << "best_new_vertex_cost: " << best_new_vertex_cost << endl;
			best_vertex_near_index = vertex_index;
//						cout << "best_vertex_near_index: " << best_vertex_near_index << endl;
			best_to_new_cost = new_cost;
			for (int j = 0; j < dubins_mark_size; j++) {
				best_dubins_mark[j] = new_dubins_mark[j];
			}
			near_radius = RRT_star_min_dis(i);
//						cout << "near_radius: " << near_radius << endl;

			if (kdtree_mode == 2) {
				//// kdtree code ////
				float* vertex_temp2 = invert_vector_to_array_float(vertex_new,
						states_number);
				vertex_near_set = kd_nearest_rangef(kd_vertices, vertex_temp2,
						near_radius);
				delete vertex_temp2;
				while (!kd_res_end(vertex_near_set)) {
					temp_vertex_near_index_p = (int*) kd_res_itemf(
							vertex_near_set, pos);
					temp_vertex_near_index = *temp_vertex_near_index_p;
					temp_vertex_near = export_vector_from_matrix(vertices,
							states_number, temp_vertex_near_index);
					temp_vertex_near_cost = vertex_cost_value(
							temp_vertex_near_index, edges, traj);
					temp_dubins_mark = forced_steer(temp_vertex_near,
							vertex_new);

					temp_to_new_cost = cost_func(temp_dubins_mark);
					if ((temp_dubins_mark[0] != -999)
							&& (temp_vertex_near_cost
									+ cost_func(temp_dubins_mark)
									< best_new_vertex_cost)) {
						best_new_vertex_cost = temp_vertex_near_cost
								+ temp_to_new_cost;
						best_vertex_near_index = temp_vertex_near_index;
						best_to_new_cost = temp_to_new_cost;
						for (int j = 0; j < dubins_mark_size; j++) {
							best_dubins_mark[j] = temp_dubins_mark[j];
						}
					}
					kd_res_next(vertex_near_set);
				}

//								cout << edges[0][i-1] << endl;
				edges[0][i - 1] = best_vertex_near_index;
//								cout << edges[0][i-1] << endl;
//								cout << traj[0][i-1] << endl;
				traj[0][i - 1] = best_to_new_cost;
				for (int j = 0; j < dubins_mark_size; j++) {
					traj[j + 1][i - 1] = best_dubins_mark[j];
				}

//								cout << traj[0][i-1] << endl;
				kd_res_free(vertex_near_set);

				/// reconnection part
				float* vertex_temp3 = invert_vector_to_array_float(vertex_new,
						states_number);
				vertex_near_set = kd_nearest_rangef(kd_vertices, vertex_temp3,
						near_radius);
				delete vertex_temp3;
				while (!kd_res_end(vertex_near_set)) {
					temp_vertex_near_index_p = (int*) kd_res_itemf(
							vertex_near_set, pos);
					temp_vertex_near_index = *temp_vertex_near_index_p;
					temp_vertex_near = export_vector_from_matrix(vertices,
							states_number, temp_vertex_near_index);
					temp_vertex_near_cost = vertex_cost_value(
							temp_vertex_near_index, edges, traj);
					temp_dubins_mark = forced_steer(vertex_new,
							temp_vertex_near);
					new_to_temp_cost = cost_func(temp_dubins_mark);

					if ((temp_dubins_mark[0] != -999)
							&& (temp_vertex_near_cost
									> best_new_vertex_cost + new_to_temp_cost)) {
//												cout << edges[0][temp_vertex_near_index-1] << endl;
						edges[0][temp_vertex_near_index - 1] = i;
//												cout << edges[0][temp_vertex_near_index-1] << endl;
//												cout << traj[0][temp_vertex_near_index-1] << endl;
						traj[0][temp_vertex_near_index - 1] = new_to_temp_cost;
//												cout << traj[0][temp_vertex_near_index-1] << endl;
						for (int j = 0; j < dubins_mark_size; j++) {
							traj[j + 1][temp_vertex_near_index - 1] =
									temp_dubins_mark[j];
						}
					}
					kd_res_next(vertex_near_set);
				}
				kd_res_free(vertex_near_set);
			} else {
				// code for non-kdtree prog.
			}

			float* vertex_temp4 = invert_vector_to_array_float(vertex_new,
					states_number);
			kd_insertf(kd_vertices, vertex_temp4, &edges[1][i - 1]);
			vertex_temp4[2] = vertex_temp4[2] + 2 * PI;
			kd_insertf(kd_vertices, vertex_temp4, &edges[1][i - 1]);
			vertex_temp4[2] = vertex_temp4[2] - 4 * PI;
			kd_insertf(kd_vertices, vertex_temp4, &edges[1][i - 1]);
			delete vertex_temp4;

			if ((i % 10000) == 0) {
				cout << i << endl;
			}
		} else {
			i = i - 1;
		}
		//		matrix_print_int(edges_row, i, edges);
	}
//		matrix_print_float(states_number, steps+1, vertices);
//		matrix_print_int(edges_row, steps, edges);
//	delete vertex_index_pointer;
//	delete temp_vertex_near_index_p;
	kd_free(kd_vertices);
}

////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// OpenGL functions//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void x_usr_movement_calc(int deltaMOveX_usr) {

	delta_center_cam = sqrt(
			pow(x_center - x_cam, 2) + pow(y_center - y_cam, 2));

	delta_x = deltaMOveX_usr * (y_center - y_cam) / delta_center_cam;
	delta_y = deltaMOveX_usr * (x_center - x_cam) / delta_center_cam;

	x_center = x_center + delta_x;
	y_center = y_center + delta_y;
	x_cam = x_cam + delta_x;
	y_cam = y_cam + delta_y;

}

void mouseButton(int button, int state, int x, int y) {

	// only start motion if the left button is pressed
	if (button == GLUT_LEFT_BUTTON) {

		// when the button is released
		if (state == GLUT_UP) {
			xOrigin_L = -1;
			yOrigin_L = -1;

		} else { // state = GLUT_DOWN
			xOrigin_L = x;
			yOrigin_L = y;
		}
	}
	if (button == GLUT_RIGHT_BUTTON) {

		// when the button is released
		if (state == GLUT_UP) {
			xOrigin_R = -1;
			yOrigin_R = -1;

		} else { // state = GLUT_DOWN
			xOrigin_R = x;
			yOrigin_R = y;
		}
	}
	if (button == 3) { // Wheel reports as button 3(scroll up) and button 4(scroll down)

		deltaMOveZ_usr = -1;

		comm_mov[0] = 0;
		comm_mov[1] = 0;
		comm_mov[2] = deltaMOveZ_usr;

//		float theta2=acos(sqrt(pow(x_up,2)+pow(y_up,2))/x_up);
//		comm_mov=vector_turn_x(comm_mov,theta2);

//		float phi2=atan2(x_up,z_up);
//		comm_mov=vector_turn_y(comm_mov,phi2);

		float theta = atan2(z_cam - z_center,
				sqrt(pow(x_cam - x_center, 2) + pow(y_cam - y_center, 2)));
		comm_mov = vector_turn_x(comm_mov, -theta);

		float psi = atan2(x_cam - x_center, y_cam - y_center);
		comm_mov = vector_turn_z(comm_mov, -psi);

		delta_x = comm_mov[0];
		delta_y = comm_mov[1];
		delta_z = comm_mov[2];

		x_center = x_center + delta_x;
		y_center = y_center + delta_y;
		z_center = z_center + delta_z;
		x_cam = x_cam + delta_x;
		y_cam = y_cam + delta_y;
		z_cam = z_cam + delta_z;

		w_up = 0;
	}

	if (button == 4) { // Wheel reports as button 3(scroll up) and button 4(scroll down)

		deltaMOveZ_usr = 1;

		comm_mov[0] = 0;
		comm_mov[1] = 0;
		comm_mov[2] = deltaMOveZ_usr;

//		float theta2=acos(sqrt(pow(x_up,2)+pow(y_up,2))/x_up);
//		comm_mov=vector_turn_x(comm_mov,theta2);

//		float phi2=atan2(x_up,z_up);
//		comm_mov=vector_turn_y(comm_mov,phi2);

		float theta = atan2(z_cam - z_center,
				sqrt(pow(x_cam - x_center, 2) + pow(y_cam - y_center, 2)));
		comm_mov = vector_turn_x(comm_mov, -theta);

		float psi = atan2(x_cam - x_center, y_cam - y_center);
		comm_mov = vector_turn_z(comm_mov, -psi);

		delta_x = comm_mov[0];
		delta_y = comm_mov[1];
		delta_z = comm_mov[2];

		x_center = x_center + delta_x;
		y_center = y_center + delta_y;
		z_center = z_center + delta_z;
		x_cam = x_cam + delta_x;
		y_cam = y_cam + delta_y;
		z_cam = z_cam + delta_z;

		w_up = 0;
	}
}

void mouseMove(int x, int y) {

	// this will only be true when the left button is down

	if (xOrigin_L >= 0) {

		deltaMOveX_usr = (x - xOrigin_L) * 0.001f;
		deltaMOveY_usr = -(y - yOrigin_L) * 0.001f;

		comm_mov[0] = deltaMOveX_usr;
		comm_mov[1] = deltaMOveY_usr;
		comm_mov[2] = 0;

//		float theta2=acos(sqrt(pow(x_up,2)+pow(y_up,2))/x_up);
//		comm_mov=vector_turn_x(comm_mov,theta2);

//		float phi2=atan2(x_up,z_up);
//		comm_mov=vector_turn_y(comm_mov,phi2);

		float theta = atan2(z_cam - z_center,
				sqrt(pow(x_cam - x_center, 2) + pow(y_cam - y_center, 2)));
		comm_mov = vector_turn_x(comm_mov, -theta);

		float psi = atan2(x_cam - x_center, y_cam - y_center);
		comm_mov = vector_turn_z(comm_mov, -psi);

		delta_x = comm_mov[0];
		delta_y = comm_mov[1];
		delta_z = comm_mov[2];

		x_center = x_center + delta_x;
		y_center = y_center + delta_y;
		z_center = z_center + delta_z;
		x_cam = x_cam + delta_x;
		y_cam = y_cam + delta_y;
		z_cam = z_cam + delta_z;

	}
	if (xOrigin_R >= 0) {

		deltaMOveX_usr = (x - xOrigin_R) * 0.001f;
		deltaMOveY_usr = -(y - yOrigin_R) * 0.001f;

		comm_mov[0] = deltaMOveX_usr;
		comm_mov[1] = deltaMOveY_usr;
		comm_mov[2] = 0;

		float theta = atan2(z_cam - z_center,
				sqrt(pow(x_cam - x_center, 2) + pow(y_cam - y_center, 2)));
		comm_mov = vector_turn_x(comm_mov, -theta);

		float psi = atan2(x_cam - x_center, y_cam - y_center);
		comm_mov = vector_turn_z(comm_mov, -psi);

		delta_x = comm_mov[0];
		delta_y = comm_mov[1];
		delta_z = comm_mov[2];

		x_cam = x_cam + delta_x;
		y_cam = y_cam + delta_y;
		z_cam = z_cam + delta_z;
	}
}

void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;
	float ratio = 1.0 * w / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45, ratio, 1, 100);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

void draw_line(vector<float> vector_a, vector<float> vector_b, float r, float g,
		float b, float width) {
	glLineWidth(width);
	glColor3f(r, g, b);
	glBegin(GL_LINES);
	glVertex3f(vector_a[0], vector_a[1], vector_a[2]);
	glVertex3f(vector_b[0], vector_b[1], vector_b[2]);
	glEnd();
}

void draw_path(vector<float> vector_a, vector<float> vector_b,
		vector<float> path_marks, float r, float g, float b, float width) {
	// traj containes a conbination of cost function value for the corisponding edge in 0 and path marks in 1-12
	//					cost value																	--- 0
	// path_marks		R, c, L,																	--- 1-3
	//					initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 4-7
	//					v_norm,																		--- 8
	//					finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 9-12

//	cout << "vector_a=[" << vector_a[0] << "," << vector_a[1] << "," << vector_a[2] << "," << vector_a[3] << "]" << endl;
//	cout << "vector_b=[" << vector_b[0] << "," << vector_b[1] << "," << vector_b[2] << "," << vector_b[3] << "]" << endl;

//		vector_print_float(states_number, vector_a);
//		vector_print_float(states_number, vector_b);
//		vector_print_float(dubins_mark_size + 1, path_marks);

	vector<float> last_vector;
	vector<float> new_vector;
	float ds_check = ds_plot;
	float dL_plot = sqrt(pow(ds_plot, 2) * (1 + 1 / pow(V_max, 2)));
	float dL_check = dL_plot;
	float R = path_marks[1];
	float c = path_marks[2];
	float path_full_length = path_marks[3];
	float initialCircle_centre1 = path_marks[4];
	float initialCircle_centre2 = path_marks[5];
	float initialArcAngle = path_marks[6];
	float phi_pinit = path_marks[7];
	float L_first_circle = initialArcAngle * R;

	last_vector.push_back(vector_a[0]);
	last_vector.push_back(vector_a[1]);
	last_vector.push_back(vector_a[3]);
	new_vector.push_back(-999);
	new_vector.push_back(-999);
	new_vector.push_back(-999);

	//check initial circle
	while ((dL_check < path_full_length) && (dL_check < L_first_circle)) {

		if ((c == 1) || (c == 3) || (c == 11) || (c == 13) || (c == 21)
				|| (c == 23)) // path is either LSL or LSR
				{
			new_vector[0] = initialCircle_centre1
					+ R * cos(dL_check / R + phi_pinit);
			new_vector[1] = initialCircle_centre2
					+ R * sin(dL_check / R + phi_pinit);
			new_vector[2] = vector_a[3] + dL_check / V_max;
		} else // path is either RSR or RSL
		{
			new_vector[0] = initialCircle_centre1
					+ R * cos(dL_check / R - phi_pinit);
			new_vector[1] = initialCircle_centre2
					- R * sin(dL_check / R - phi_pinit);
			new_vector[2] = vector_a[3] + dL_check / V_max;
		}

//			cout << "last_vector=[" << last_vector[0] << "," << last_vector[1]
//					<< "," << last_vector[2] << "]" << endl;
//			cout << "new_vector=[" << new_vector[0] << "," << new_vector[1]
//					<< "," << new_vector[2] << "]" << endl;
		draw_line(last_vector, new_vector, r, g, b, width);
		last_vector = vector_copy_float(3, new_vector);
		dL_check = dL_check + dL_plot;
	}

	if (path_marks[8] == -999) {
		new_vector[0] = vector_b[0];
		new_vector[1] = vector_b[1];
		new_vector[2] = vector_b[3];
		draw_line(last_vector, new_vector, r, g, b, width);
		return;
	}

	//check straight line
	float L_straight = path_marks[8];
	float finalCircle_centre1 = path_marks[9];
	float finalCircle_centre2 = path_marks[10];
	float finalArcAngle = path_marks[11];
	float phi_contact = path_marks[12];

	vector<float> line_start;
	vector<float> line_end;
	if ((c == 1) || (c == 3) || (c == 11) || (c == 13) || (c == 21)
			|| (c == 23)) // start circle is L
			{
		line_start.push_back(
				initialCircle_centre1 + R * cos(initialArcAngle + phi_pinit));
		line_start.push_back(
				initialCircle_centre2 + R * sin(initialArcAngle + phi_pinit));
	} else {
		line_start.push_back(
				initialCircle_centre1 + R * cos(initialArcAngle - phi_pinit));
		line_start.push_back(
				initialCircle_centre2 - R * sin(initialArcAngle - phi_pinit));
	}

	if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
			|| (c == 24)) // end circle is L
			{
		line_end.push_back(finalCircle_centre1 + R * cos(phi_contact));
		line_end.push_back(finalCircle_centre2 + R * sin(phi_contact));
	} else {
		line_end.push_back(finalCircle_centre1 + R * cos(-phi_contact));
		line_end.push_back(finalCircle_centre2 - R * sin(-phi_contact));
	}

	while ((dL_check < path_full_length)
			&& (dL_check < L_first_circle + L_straight)) {

		vector<float> line_cut = linear_interpulation_2d(line_start, line_end,
				dL_check - L_first_circle);

		new_vector[0] = line_cut[0];
		new_vector[1] = line_cut[1];
		new_vector[2] = vector_a[3] + dL_check / V_max;

//			cout << "last_vector=[" << last_vector[0] << "," << last_vector[1]
//					<< "," << last_vector[2] << "]" << endl;
//			cout << "new_vector=[" << new_vector[0] << "," << new_vector[1]
//					<< "," << new_vector[2] << "]" << endl;

		draw_line(last_vector, new_vector, r, g, b, width);
		last_vector = vector_copy_float(3, new_vector);
		dL_check = dL_check + dL_plot;
	}

	if (path_marks[11] == -999) {
		new_vector[0] = vector_b[0];
		new_vector[1] = vector_b[1];
		new_vector[2] = vector_b[3];
		draw_line(last_vector, new_vector, r, g, b, width);
		return;
	}

	//check final circle
	float L_last_circle = finalArcAngle * R;

	while ((dL_check < path_full_length)
			&& (dL_check < L_first_circle + L_straight + L_last_circle)) {

		float dL_last_circle = dL_check - L_first_circle - L_straight;
		if ((c == 1) || (c == 4) || (c == 11) || (c == 14) || (c == 21)
				|| (c == 24)) // path is either LSL or RSL
				{
			new_vector[0] = finalCircle_centre1
					+ R * cos(dL_last_circle / R + phi_contact);
			new_vector[1] = finalCircle_centre2
					+ R * sin(dL_last_circle / R + phi_contact);
			new_vector[2] = vector_a[3] + dL_check / V_max;
		} else {
			new_vector[0] = finalCircle_centre1
					+ R * cos(dL_last_circle / R - phi_contact);
			new_vector[1] = finalCircle_centre2
					- R * sin(dL_last_circle / R - phi_contact);
			new_vector[2] = vector_a[3] + dL_check / V_max;
		}

//			cout << "last_vector=[" << last_vector[0] << "," << last_vector[1]
//					<< "," << last_vector[2] << "]" << endl;
//			cout << "new_vector=[" << new_vector[0] << "," << new_vector[1]
//					<< "," << new_vector[2] << "]" << endl;

		draw_line(last_vector, new_vector, r, g, b, width);
		last_vector = vector_copy_float(3, new_vector);
		dL_check = dL_check + dL_plot;
	}

	new_vector[0] = vector_b[0];
	new_vector[1] = vector_b[1];
	new_vector[2] = path_full_length / V_max;
//		cout << "last_vector=[" << last_vector[0] << "," << last_vector[1]
//				<< "," << last_vector[2] << "]" << endl;
//		cout << "new_vector=[" << new_vector[0] << "," << new_vector[1] << ","
//				<< new_vector[2] << "]" << endl;
	if (new_vector[2] > last_vector[2])
		draw_line(last_vector, new_vector, r, g, b, width);
	return;
}

void draw_RRT_graph(float** vertices, int** edges, float** traj) {

	int start_index, end_index;
	vector<float> vector_start(states_number);
	vector<float> vector_end(states_number);
	vector<float> current_traj(traj_row);

	float max_value = plain[3][1];
	float min_value = plain[3][0];

	float color_section = (max_value - min_value) / 64;
	int color_index;

	for (int i = 0; i < edges_col; i++) {
		start_index = edges[0][i];
		end_index = edges[1][i];
		vector_start = export_vector_from_matrix(vertices, states_number,
				start_index);
		vector_end = export_vector_from_matrix(vertices, states_number,
				end_index);
		float time_avg = (vector_end[3] + vector_start[3]) / 2.0;
		color_index = (time_avg - min_value) / color_section;
		current_traj = export_vector_from_matrix(traj, traj_row, i);
		draw_path(vector_start, vector_end, current_traj,
				colormap[0][color_index], colormap[1][color_index],
				colormap[2][color_index], 1.5);
	}
	//	vector_start[0]=5;
	//	vector_start[1]=5;
	//	draw_line_for_tree(vector_start,vector_end);

}

void draw_solution_path(float** vertices, int** edges, float** traj) {

	int temp_index;
	if (Dubins_sol_mode == 1)
		temp_index = hit_index_RD;
	else
		temp_index = hit_index;
	int prev_index, end_index;
	vector<float> vector_prev(states_number);
	vector<float> vector_temp(states_number);
	vector<float> current_traj(traj_row);

	while (temp_index != 0) {
		prev_index = edges[0][temp_index - 1];
		vector_prev = export_vector_from_matrix(vertices, states_number,
				prev_index);
		vector_temp = export_vector_from_matrix(vertices, states_number,
				temp_index);
		current_traj = export_vector_from_matrix(traj, traj_row,
				temp_index - 1);
		draw_path(vector_prev, vector_temp, current_traj, 0.0, 0.0, 0.0, 4.0);
		temp_index = prev_index;
	}
	//	vector_start[0]=5;
	//	vector_start[1]=5;
	//	draw_line_for_tree(vector_start,vector_end);

}

void draw_monte_carlo_solutions(float** solutions_states_relaxed_Dubins,
		float** solutions_marks_relaxed_Dubins,
		float** solutions_states_full_Dubins,
		float** solutions_marks_full_Dubins,
		float** solutions_TTH_Relaxed_Dubins_matrix,
		float** solutions_TTH_full_Dubins_matrix,
		vector<float> monte_carlo_avg_Relaxed_Dubins,
		vector<float> monte_carlo_avg_full_Dubins,
		vector<float> convergence_values) {

	//	plot_mode //  1 - all relaxed dubins solutions  2 - only best relaxed dubins solution
	//                3 - all full dubins solutions 	4 - only best full dubins solution
	//				  5 - best euclidean, relaxed and full solutions with optimal solutions

	int soluitons_matrix_relaxed_dubins_size = sizes_values_final[0];
	int marks_matrix_relaxed_dubins_size = sizes_values_final[1];
	int soluitons_matrix_full_dubins_size = sizes_values_final[2];
	int marks_matrix_full_dubins_size = sizes_values_final[3];

	vector<float> vector_prev, vector_next;
	for (int i = 0; i < states_number; i++) {
		vector_prev.push_back(-999);
		vector_next.push_back(-999);
	}

	vector<float> current_traj;
	for (int i = 0; i < dubins_mark_size + 1; i++) {
		current_traj.push_back(-999);
	}

	double r, g, b;

	if (plot_mode == 1) {
		for (int i = 0; i < monte_carlo_number; i++) {
			for (int mone = 1; mone < soluitons_matrix_relaxed_dubins_size;
					mone++) {
				if (solutions_states_relaxed_Dubins[i * states_number][mone]
						!= -999) {
					for (int j = 0; j < states_number; j++) {
						vector_prev[j] = solutions_states_relaxed_Dubins[i
								* states_number + j][mone - 1];
						vector_next[j] = solutions_states_relaxed_Dubins[i
								* states_number + j][mone];
					}
					for (int j = 1; j < dubins_mark_size + 1; j++) {
						current_traj[j] =
								solutions_marks_relaxed_Dubins[i][(mone - 1)
										* dubins_mark_size + j - 1];
					}
//						if (mone == 13) {
//							cout << "13" << endl;
//							vector_print_float(states_number, vector_prev);
//							vector_print_float(states_number, vector_next);
//							vector_print_float(dubins_mark_size + 1,
//									current_traj);
//						}
					r = color_convergence_matrix[i][0];
					g = color_convergence_matrix[i][1];
					b = color_convergence_matrix[i][2];
					draw_path(vector_prev, vector_next, current_traj, r, g, b,
							4.0);
				}
			}
		}
	}

	if (plot_mode == 2) {
		int i = 0;
		float best_time =
				solutions_TTH_Relaxed_Dubins_matrix[0][convergence_size - 1];

		for (int m_c_number = 1; m_c_number < monte_carlo_number;
				m_c_number++) {
			if (solutions_TTH_Relaxed_Dubins_matrix[m_c_number][convergence_size
					- 1] < best_time) {
				best_time =
						solutions_TTH_Relaxed_Dubins_matrix[m_c_number][convergence_size
								- 1];
				i = m_c_number;
			}
		}

		for (int mone = 1; mone < soluitons_matrix_relaxed_dubins_size;
				mone++) {
			if (solutions_states_relaxed_Dubins[i * states_number][mone]
					!= -999) {
				for (int j = 0; j < states_number; j++) {
					vector_prev[j] = solutions_states_relaxed_Dubins[i
							* states_number + j][mone - 1];
					vector_next[j] = solutions_states_relaxed_Dubins[i
							* states_number + j][mone];
				}
				for (int j = 1; j < dubins_mark_size + 1; j++) {
					current_traj[j] = solutions_marks_relaxed_Dubins[i][(mone
							- 1) * dubins_mark_size + j - 1];
				}
//						if (mone == 13) {
//							cout << "13" << endl;
//							vector_print_float(states_number, vector_prev);
//							vector_print_float(states_number, vector_next);
//							vector_print_float(dubins_mark_size + 1,
//									current_traj);
//						}
				r = color_convergence_matrix[i][0];
				g = color_convergence_matrix[i][1];
				b = color_convergence_matrix[i][2];
				draw_path(vector_prev, vector_next, current_traj, r, g, b, 4.0);
			}
		}
	}

	if (plot_mode == 3) {
		for (int i = 0; i < monte_carlo_number; i++) {
			for (int mone = 1; mone < soluitons_matrix_full_dubins_size;
					mone++) {
				if (solutions_states_full_Dubins[i * states_number][mone]
						!= -999) {
					for (int j = 0; j < states_number; j++) {
						vector_prev[j] = solutions_states_full_Dubins[i
								* states_number + j][mone - 1];
						vector_next[j] = solutions_states_full_Dubins[i
								* states_number + j][mone];
					}
					for (int j = 1; j < dubins_mark_size + 1; j++) {
						current_traj[j] = solutions_marks_full_Dubins[i][(mone
								- 1) * dubins_mark_size + j - 1];
					}
//						if (mone == 13) {
//							cout << "13" << endl;
//							vector_print_float(states_number, vector_prev);
//							vector_print_float(states_number, vector_next);
//							vector_print_float(dubins_mark_size + 1,
//									current_traj);
//						}
					r = color_convergence_matrix[i][0];
					g = color_convergence_matrix[i][1];
					b = color_convergence_matrix[i][2];
					draw_path(vector_prev, vector_next, current_traj, r, g, b,
							4.0);
				}
			}
		}
	}

	if (plot_mode == 4) {
		int i = 0;
		float best_time = solutions_TTH_full_Dubins_matrix[0][convergence_size
				- 1];

		for (int m_c_number = 1; m_c_number < monte_carlo_number;
				m_c_number++) {
			if (solutions_TTH_full_Dubins_matrix[m_c_number][convergence_size
					- 1] < best_time) {
				best_time =
						solutions_TTH_full_Dubins_matrix[m_c_number][convergence_size
								- 1];
				i = m_c_number;
			}
		}

		for (int mone = 1; mone < soluitons_matrix_full_dubins_size; mone++) {
			if (solutions_states_full_Dubins[i * states_number][mone] != -999) {
				for (int j = 0; j < states_number; j++) {
					vector_prev[j] = solutions_states_full_Dubins[i
							* states_number + j][mone - 1];
					vector_next[j] = solutions_states_full_Dubins[i
							* states_number + j][mone];
				}
				for (int j = 1; j < dubins_mark_size + 1; j++) {
					current_traj[j] = solutions_marks_full_Dubins[i][(mone - 1)
							* dubins_mark_size + j - 1];
				}
//						if (mone == 13) {
//							cout << "13" << endl;
//							vector_print_float(states_number, vector_prev);
//							vector_print_float(states_number, vector_next);
//							vector_print_float(dubins_mark_size + 1,
//									current_traj);
//						}
				r = color_convergence_matrix[i][0];
				g = color_convergence_matrix[i][1];
				b = color_convergence_matrix[i][2];
				draw_path(vector_prev, vector_next, current_traj, r, g, b, 4.0);
			}
		}
	}

	if (plot_mode == 5) {
		// plot best relaxed dubins solution
		int i = 0;
		float best_time =
				solutions_TTH_Relaxed_Dubins_matrix[0][convergence_size - 1];

		for (int m_c_number = 1; m_c_number < monte_carlo_number;
				m_c_number++) {
			if (solutions_TTH_Relaxed_Dubins_matrix[m_c_number][convergence_size
					- 1] < best_time) {
				best_time =
						solutions_TTH_Relaxed_Dubins_matrix[m_c_number][convergence_size
								- 1];
				i = m_c_number;
			}
		}

		for (int mone = 1; mone < soluitons_matrix_relaxed_dubins_size;
				mone++) {
			if (solutions_states_relaxed_Dubins[i * states_number][mone]
					!= -999) {
				for (int j = 0; j < states_number; j++) {
					vector_prev[j] = solutions_states_relaxed_Dubins[i
							* states_number + j][mone - 1];
					vector_next[j] = solutions_states_relaxed_Dubins[i
							* states_number + j][mone];
				}
				for (int j = 1; j < dubins_mark_size + 1; j++) {
					current_traj[j] = solutions_marks_relaxed_Dubins[i][(mone
							- 1) * dubins_mark_size + j - 1];
				}
//						if (mone == 13) {
//							cout << "13" << endl;
//							vector_print_float(states_number, vector_prev);
//							vector_print_float(states_number, vector_next);
//							vector_print_float(dubins_mark_size + 1,
//									current_traj);
//						}
				draw_path(vector_prev, vector_next, current_traj, 1.0, 0.0, 0.0,
						4.0);
			}
		}

		// plot best full dubins solution

		i = 0;
		best_time = solutions_TTH_full_Dubins_matrix[0][convergence_size - 1];

		for (int m_c_number = 1; m_c_number < monte_carlo_number;
				m_c_number++) {
			if (solutions_TTH_full_Dubins_matrix[m_c_number][convergence_size
					- 1] < best_time) {
				best_time =
						solutions_TTH_full_Dubins_matrix[m_c_number][convergence_size
								- 1];
				i = m_c_number;
			}
		}

		for (int mone = 1; mone < soluitons_matrix_full_dubins_size; mone++) {
			if (solutions_states_full_Dubins[i * states_number][mone] != -999) {
				for (int j = 0; j < states_number; j++) {
					vector_prev[j] = solutions_states_full_Dubins[i
							* states_number + j][mone - 1];
					vector_next[j] = solutions_states_full_Dubins[i
							* states_number + j][mone];
				}
				for (int j = 1; j < dubins_mark_size + 1; j++) {
					current_traj[j] = solutions_marks_full_Dubins[i][(mone - 1)
							* dubins_mark_size + j - 1];
				}
//						if (mone == 13) {
//							cout << "13" << endl;
//							vector_print_float(states_number, vector_prev);
//							vector_print_float(states_number, vector_next);
//							vector_print_float(dubins_mark_size + 1,
//									current_traj);
//						}
				r = color_convergence_matrix[i][0];
				g = color_convergence_matrix[i][1];
				b = color_convergence_matrix[i][2];
				draw_path(vector_prev, vector_next, current_traj, 0.0, 1.0, 0.0,
						4.0);
			}
		}

		// plot optimal euclidean solution

		vector<float> start_vector, end_vector;

		start_vector.push_back(x1_0);
		start_vector.push_back(x2_0);
		start_vector.push_back(x4_0);

		end_vector.push_back(goal_point_x1);
		end_vector.push_back(goal_point_x2);
		end_vector.push_back(euclidean_dis / V_max);

		draw_line(start_vector, end_vector, 0.0, 0.0, 1.0, 4.0);

		// plot optimal full dubins solution

		vector<float> pinit, pf;
		float thinit, thf;
		pinit.push_back(x1_0);
		pinit.push_back(x2_0);
		thinit = x3_0;
		pf.push_back(goal_point_x1);
		pf.push_back(goal_point_x2);
		thf = goal_point_x3;

		float type = dubins(pinit, thinit, pf, thf, Sim_R_min);
//								cout << type << endl;

		//RSR_length - 1
		//RLR_length - 2
		//RSL_length - 3
		//LSL_length - 4
		//LRL_length - 5
		//LSR_length - 6

		current_traj[1] = Sim_R_min;
		vector<float> full_dubins_marks;
		if (type == 1) {
			full_dubins_marks = dubins_RSR(pinit, thinit, pf, thf, Sim_R_min);
			current_traj[2] = 2;
		}
		if (type == 3) {
			full_dubins_marks = dubins_RSL(pinit, thinit, pf, thf, Sim_R_min);
			current_traj[2] = 4;
		}
		if (type == 4) {
			full_dubins_marks = dubins_LSL(pinit, thinit, pf, thf, Sim_R_min);
			current_traj[2] = 1;
		}
		if (type == 5) {
			full_dubins_marks = dubins_LSR(pinit, thinit, pf, thf, Sim_R_min);
			current_traj[2] = 3;
		}
		// path_marks_and_vector= {	R, c, L,																	--- 0-2
		//							initialCircle_centre[0],initialCircle_centre[1],initialArcAngle,phi_pinit,	--- 3-6
		//							v_norm,																		--- 7
		//							finalCircle_centre[0],finalCircle_centre[1],finalArcAngle,phi_contact,		--- 8-11
		//							new_vertex[0],new_vertex[1],new_vertex[2],new_vertex[3]}					--- 12-15

//								vector_print_float(dubins_mark_size-2, full_dubins_marks);

		for (int i = 3; i < dubins_mark_size + 1; i++) {
			current_traj[i] = full_dubins_marks[i - 3];
		}

		start_vector[2] = x3_0;
		start_vector.push_back(x4_0);
		end_vector[2] = goal_point_x3;
		end_vector.push_back(goal_point_x4);

//						vector_print_float(states_number, start_vector);
//						vector_print_float(states_number, end_vector);
//						vector_print_float(dubins_mark_size+1, current_traj);

		draw_path(start_vector, end_vector, current_traj, 0.0, 0.7, 0.0, 4.0);

	}

	//	vector_start[0]=5;
	//	vector_start[1]=5;
	//	draw_line_for_tree(vector_start,vector_end);

}

void draw_floor() {
	glLineWidth(1);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (GLfloat i = -20.0; i <= 40.0; i += 2.0) {
		glVertex3f(i, 40.0, 0);
		glVertex3f(i, -20.0, 0);
		glVertex3f(40.0, i, 0);
		glVertex3f(-20.0, i, 0);
	}
}

void draw_plain() {
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][0], plain[1][0], plain[3][0]);
	glVertex3f(plain[0][1], plain[1][0], plain[3][0]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][1], plain[1][0], plain[3][0]);
	glVertex3f(plain[0][1], plain[1][1], plain[3][0]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][1], plain[1][1], plain[3][0]);
	glVertex3f(plain[0][0], plain[1][1], plain[3][0]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][0], plain[1][1], plain[3][0]);
	glVertex3f(plain[0][0], plain[1][0], plain[3][0]);
	glEnd();

	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][0], plain[1][0], plain[3][0]);
	glVertex3f(plain[0][0], plain[1][0], plain[3][1]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][1], plain[1][0], plain[3][0]);
	glVertex3f(plain[0][1], plain[1][0], plain[3][1]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][1], plain[1][1], plain[3][0]);
	glVertex3f(plain[0][1], plain[1][1], plain[3][1]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][0], plain[1][1], plain[3][0]);
	glVertex3f(plain[0][0], plain[1][1], plain[3][1]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][0], plain[1][0], plain[3][1]);
	glVertex3f(plain[0][1], plain[1][0], plain[3][1]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][1], plain[1][0], plain[3][1]);
	glVertex3f(plain[0][1], plain[1][1], plain[3][1]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][1], plain[1][1], plain[3][1]);
	glVertex3f(plain[0][0], plain[1][1], plain[3][1]);
	glEnd();
	glLineWidth(2.5);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(plain[0][0], plain[1][1], plain[3][1]);
	glVertex3f(plain[0][0], plain[1][0], plain[3][1]);
	glEnd();
}

void draw_obstacles() {

	vector<float> obst_color;
	obst_color.push_back(1);
	obst_color.push_back(0);
	obst_color.push_back(1);

	for (int mone = 1; mone < obst_num + 1; mone++) {

		float time_max = plain[3][1];
		float time_min = plain[3][0];
		float d_time = time_max - time_min;
		float d_time_temp = 100000;
		float d_time_temp2 = 100000;

		// define 4 low points

		float lowminmin_x = obstacles[(mone - 1) * (states_number - 1) + 0][0];
		float lowminmin_y = obstacles[(mone - 1) * (states_number - 1) + 1][0];
		float lowminmin_z = time_min;

		float lowminmax_x = obstacles[(mone - 1) * (states_number - 1) + 0][0];
		float lowminmax_y = obstacles[(mone - 1) * (states_number - 1) + 1][1];
		float lowminmax_z = time_min;

		float lowmaxmin_x = obstacles[(mone - 1) * (states_number - 1) + 0][1];
		float lowmaxmin_y = obstacles[(mone - 1) * (states_number - 1) + 1][0];
		float lowmaxmin_z = time_min;

		float lowmaxmax_x = obstacles[(mone - 1) * (states_number - 1) + 0][1];
		float lowmaxmax_y = obstacles[(mone - 1) * (states_number - 1) + 1][1];
		float lowmaxmax_z = time_min;

		// define 4 high points

		float highminmin_x = obstacles[(mone - 1) * (states_number - 1) + 0][0]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 0][2];
		float highminmin_y = obstacles[(mone - 1) * (states_number - 1) + 1][0]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 1][2];
		float highminmin_z = time_max;
		if (highminmin_x < plain[0][0])
			d_time_temp = (plain[0][0] - lowminmin_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highminmin_x > plain[0][1])
			d_time_temp = (plain[0][1] - lowminmin_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highminmin_y < plain[1][0])
			d_time_temp2 = (plain[1][0] - lowminmin_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (highminmin_y > plain[1][1])
			d_time_temp2 = (plain[1][1] - lowminmin_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (d_time_temp2 < d_time_temp)
			d_time_temp = d_time_temp2;
		if (d_time_temp < d_time) {
			highminmin_x =
					obstacles[(mone - 1) * (states_number - 1) + 0][0]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 0][2];
			highminmin_y =
					obstacles[(mone - 1) * (states_number - 1) + 1][0]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 1][2];
			highminmin_z = time_min + d_time_temp;
			d_time_temp = 100000;
			d_time_temp2 = 100000;
		}

		float highminmax_x = obstacles[(mone - 1) * (states_number - 1) + 0][0]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 0][2];
		float highminmax_y = obstacles[(mone - 1) * (states_number - 1) + 1][1]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 1][2];
		float highminmax_z = time_max;
		if (highminmax_x < plain[0][0])
			d_time_temp = (plain[0][0] - lowminmax_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highminmax_x > plain[0][1])
			d_time_temp = (plain[0][1] - lowminmax_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highminmax_y < plain[1][0])
			d_time_temp2 = (plain[1][0] - lowminmax_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (highminmax_y > plain[1][1])
			d_time_temp2 = (plain[1][1] - lowminmax_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (d_time_temp2 < d_time_temp)
			d_time_temp = d_time_temp2;
		if (d_time_temp < d_time) {
			highminmax_x =
					obstacles[(mone - 1) * (states_number - 1) + 0][0]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 0][2];
			highminmax_y =
					obstacles[(mone - 1) * (states_number - 1) + 1][1]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 1][2];
			highminmax_z = time_min + d_time_temp;
			d_time_temp = 100000;
			d_time_temp2 = 100000;
		}

		float highmaxmin_x = obstacles[(mone - 1) * (states_number - 1) + 0][1]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 0][2];
		float highmaxmin_y = obstacles[(mone - 1) * (states_number - 1) + 1][0]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 1][2];
		float highmaxmin_z = time_max;
		if (highmaxmin_x < plain[0][0])
			d_time_temp = (plain[0][0] - lowmaxmin_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highmaxmin_x > plain[0][1])
			d_time_temp = (plain[0][1] - lowmaxmin_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highmaxmin_y < plain[1][0])
			d_time_temp2 = (plain[1][0] - lowmaxmin_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (highmaxmin_y > plain[1][1])
			d_time_temp2 = (plain[1][1] - lowmaxmin_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (d_time_temp2 < d_time_temp)
			d_time_temp = d_time_temp2;
		if (d_time_temp < d_time) {
			highmaxmin_x =
					obstacles[(mone - 1) * (states_number - 1) + 0][1]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 0][2];
			highmaxmin_y =
					obstacles[(mone - 1) * (states_number - 1) + 1][0]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 1][2];
			highmaxmin_z = time_min + d_time_temp;
			d_time_temp = 100000;
			d_time_temp2 = 100000;
		}

		float highmaxmax_x = obstacles[(mone - 1) * (states_number - 1) + 0][1]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 0][2];
		float highmaxmax_y = obstacles[(mone - 1) * (states_number - 1) + 1][1]
				+ d_time * obstacles[(mone - 1) * (states_number - 1) + 1][2];
		float highmaxmax_z = time_max;
		if (highmaxmax_x < plain[0][0])
			d_time_temp = (plain[0][0] - lowmaxmax_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highmaxmax_x > plain[0][1])
			d_time_temp = (plain[0][1] - lowmaxmax_x)
					/ obstacles[(mone - 1) * (states_number - 1) + 0][2];
		if (highmaxmax_y < plain[1][0])
			d_time_temp2 = (plain[1][0] - lowmaxmax_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (highmaxmax_y > plain[1][1])
			d_time_temp2 = (plain[1][1] - lowmaxmax_y)
					/ obstacles[(mone - 1) * (states_number - 1) + 1][2];
		if (d_time_temp2 < d_time_temp)
			d_time_temp = d_time_temp2;
		if (d_time_temp < d_time) {
			highmaxmax_x =
					obstacles[(mone - 1) * (states_number - 1) + 0][1]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 0][2];
			highmaxmax_y =
					obstacles[(mone - 1) * (states_number - 1) + 1][1]
							+ d_time_temp
									* obstacles[(mone - 1) * (states_number - 1)
											+ 1][2];
			highmaxmax_z = time_min + d_time_temp;
			d_time_temp = 100000;
			d_time_temp2 = 100000;
		}

		glBegin(GL_POLYGON); //down
		glColor3f(obst_color[0], obst_color[1], obst_color[2]);
		glVertex3f(lowminmin_x, lowminmin_y, lowminmin_z);
		glVertex3f(lowminmax_x, lowminmax_y, lowminmax_z);
		glVertex3f(lowmaxmax_x, lowmaxmax_y, lowmaxmax_z);
		glVertex3f(lowmaxmin_x, lowmaxmin_y, lowmaxmin_z);
		glEnd();

		glBegin(GL_POLYGON); //up
		glColor3f(obst_color[0], obst_color[1], obst_color[2]);
		glVertex3f(highminmin_x, highminmin_y, highminmin_z);
		glVertex3f(highminmax_x, highminmax_y, highminmax_z);
		glVertex3f(highmaxmax_x, highmaxmax_y, highmaxmax_z);
		glVertex3f(highmaxmin_x, highmaxmin_y, highmaxmin_z);

		glEnd();

		glBegin(GL_POLYGON); //x-z back
		glColor3f(obst_color[0], obst_color[1], obst_color[2]);
		glVertex3f(lowminmin_x, lowminmin_y, lowminmin_z);
		glVertex3f(lowmaxmin_x, lowmaxmin_y, lowmaxmin_z);
		glVertex3f(highmaxmin_x, highmaxmin_y, highmaxmin_z);
		glVertex3f(highminmin_x, highminmin_y, highminmin_z);
		glEnd();

		glBegin(GL_POLYGON); //x-z front
		glColor3f(obst_color[0], obst_color[1], obst_color[2]);
		glVertex3f(lowminmax_x, lowminmax_y, lowminmax_z);
		glVertex3f(lowmaxmax_x, lowmaxmax_y, lowmaxmax_z);
		glVertex3f(highmaxmax_x, highmaxmax_y, highmaxmax_z);
		glVertex3f(highminmax_x, highminmax_y, highminmax_z);
		glEnd();

		glBegin(GL_POLYGON); //y-z back
		glColor3f(obst_color[0], obst_color[1], obst_color[2]);
		glVertex3f(lowminmin_x, lowminmin_y, lowminmin_y);
		glVertex3f(lowminmax_x, lowminmax_y, lowminmax_z);
		glVertex3f(highminmax_x, highminmax_y, highminmax_z);
		glVertex3f(highminmin_x, highminmin_y, highminmin_z);
		glEnd();

		glBegin(GL_POLYGON); //y-z front
		glColor3f(obst_color[0], obst_color[1], obst_color[2]);
		glVertex3f(lowmaxmin_x, lowmaxmin_y, lowmaxmin_z);
		glVertex3f(lowmaxmax_x, lowmaxmax_y, lowmaxmax_z);
		glVertex3f(highmaxmax_x, highmaxmax_y, highmaxmax_z);
		glVertex3f(highmaxmin_x, highmaxmin_y, highmaxmin_z);
		glEnd();
	}
}

void drawing_full(void) {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();

	//	gluLookAt(	x, 1.0f, z,
	//			x+lx, 1.0f,  z+lz,
	//			0.W0f, 1.0f,  0.0f);
	gluLookAt(x_cam, y_cam, z_cam, x_center, y_center, z_center, x_up, y_up,
			z_up);

	//	const float plain[states_number][plain_col] = { { 0, 20 },	//x1 - defines edges of the potential search space
	//	{ 0, 20 },		//x2 - plain row i is made of two numbers that state the min and max
	//	{ 0, 40 } };	//x3 = time - allowed values of the i state

	draw_floor();
	draw_plain();
	draw_obstacles();

//	plot_mode //  1 - only tree  2- only best solution  3 - full tree and solution
	if (solution_mode == 1) {
		if (plot_mode == 1)
			draw_RRT_graph(vertices_final, edges_final, traj_final);
		if (plot_mode == 2) {
			draw_solution_path(vertices_final, edges_final, traj_final);
		}
		if ((plot_mode != 1) && (plot_mode != 2)) {
			draw_RRT_graph(vertices_final, edges_final, traj_final);
			draw_solution_path(vertices_final, edges_final, traj_final);
		}
	}

	//	plot_mode //  1 - all relaxed dubins solutions  2 - only best relaxed dubins solution
	//                3 - all full dubins solutions 	4 - only best full dubins solution
	//				  5 - best euclidean, relaxed and full solutions with optimal solutions
	if (solution_mode == 2) {
		draw_monte_carlo_solutions(solutions_states_relaxed_Dubins_final,
				solutions_marks_relaxed_Dubins_final,
				solutions_states_full_Dubins_final,
				solutions_marks_full_Dubins_final,
				solutions_TTH_Relaxed_Dubins_matrix_final,
				solutions_TTH_full_Dubins_matrix_final,
				monte_carlo_avg_Relaxed_Dubins_final,
				monte_carlo_avg_full_Dubins_final, convergence_values_final);
	}

	glutSwapBuffers();
}

void processNormalKeys(unsigned char key, int x, int y) {

	if (key == 27)
		exit(0);
}

//void computePos(float deltaMove) {
//
//	x += deltaMove * lx * 0.0f;
//	z += deltaMove * lz * 0.1f;
//}

void renderScene() {

	//	if (deltaMove)
	//		computePos(deltaMove);

	// Clear Color and Depth Buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();
	// Set the camera
	//	gluLookAt(	x, 1.0f, z,
	//			x+lx, 1.0f,  z+lz,
	//			0.0f, 1.0f,  0.0f);

	glPushMatrix();
	//	glTranslatef(i*10.0,0,j * 10.0);
	drawing_full();
	glPopMatrix();
	glutSwapBuffers();
}

void void_func(void) {

}

unsigned int get_msec(void) {
	static struct timeval timeval, first_timeval;

	gettimeofday(&timeval, 0);

	if (first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000
			+ (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}

//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Dubins functions//////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

////////
vector<float> vector_subtraction(vector<float> a, vector<float> b) {
	vector<float> c;
	int a_size = a.size();
	int b_size = b.size();

	int c_size;
	if (a_size < b_size)
		c_size = a_size;
	else
		c_size = b_size;

	for (int i = 0; i < c_size; i++) {
		c.push_back(a[i] - b[i]);
	}
	return c;
}
////////
vector<float> vector_adding(vector<float> a, vector<float> b) {
	vector<float> c;
	int a_size = a.size();
	int b_size = b.size();

	int c_size;
	if (a_size < b_size)
		c_size = a_size;
	else
		c_size = b_size;

	for (int i = 0; i < c_size; i++) {
		c.push_back(a[i] + b[i]);
	}
	return c;
}
////////
vector<float> vector_multiplication(vector<float> a, vector<float> b) {
	vector<float> c;
	int a_size = a.size();
	int b_size = b.size();

	int c_size;
	if (a_size < b_size)
		c_size = a_size;
	else
		c_size = b_size;

	for (int i = 0; i < c_size; i++) {
		c.push_back(a[i] * b[i]);
	}
	return c;
}
////////
vector<float> vector_division(vector<float> a, vector<float> b) {
	vector<float> c;
	int a_size = a.size();
	int b_size = b.size();

	int c_size;
	if (a_size < b_size)
		c_size = a_size;
	else
		c_size = b_size;

	for (int i = 0; i < c_size; i++) {
		if (b[i] == 0)
			c.push_back(-999);
		else
			c.push_back(a[i] / b[i]);
	}
	return c;
}

////////
vector<float> vector_scalar_addition(float scalar, vector<float> a) {
	vector<float> c;
	int a_size = a.size();

	for (int i = 0; i < a_size; i++) {
		c.push_back(scalar + a[i]);
	}
	return c;
}
////////
vector<float> vector_scalar_multiplication(float scalar, vector<float> a) {
	vector<float> c;
	int a_size = a.size();

	for (int i = 0; i < a_size; i++) {
		c.push_back(scalar * a[i]);
	}
	return c;
}
////////
vector<float> vector_scalar_power(vector<float> a, float scalar) {
	vector<float> c;
	int a_size = a.size();

	for (int i = 0; i < a_size; i++) {
		c.push_back(pow(a[i], scalar));
	}
	return c;
}

////////
float vector_sum(vector<float> a) {
	float sum = 0;
	int a_size = a.size();

	for (int i = 0; i < a_size; i++) {
		sum = sum + a[i];
	}
	return sum;
}
////////
float vector_norm(vector<float> a) {
	vector<float> b = vector_scalar_power(a, 2);

	float sum = vector_sum(b);

	float sum2 = sqrt(sum);

	return sum2;
}

////////
float atan2_0_2pi(float y, float x) {
	float th;
	if (atan2(y, x) < 0)
		th = 2 * PI + atan2(y, x);
	else
		th = atan2(y, x);
	return th;
}

float second_order_pol_solution(float A, float B, float C, float xmin,
		float xmax)

// this function gives the smallest soluiton for the equation Ax ^ 2 + Bx + C = 0
// in the interval[xmin, xmax] if one exists.
// If not the function returns - 999
		{
	float delta = pow(B, 2) - 4 * A * C;

	float x = -999;
	if (delta < 0) {
		return x;
	}

	float x1 = (-B + pow(delta, 0.5)) / 2 / A;
	float x2 = (-B - pow(delta, 0.5)) / 2 / A;

	if ((x1 >= xmin) && (x1 <= xmax))
		x = x1;

	if ((x2 >= xmin) && (x2 <= xmax)) {
		if ((x == x1) && (x2 < x1))
			x = x2;
		if ((x == -999))
			x = x2;
	}
	return x;
}

vector<float> second_order_equ_calc(float x1, float y1, float x2, float y2,
		float x3, float y3)

// this function gives coefficients for the equation y=Ax ^ 2 + Bx + C
// using 3 points (x1,y1) (x2,y2) (x3,y3)

		{
	vector<float> abc_vec;
	float A = ((y2 - y1) * (x3 - x1) - (y3 - y1) * (x2 - x1))
			/ ((pow(x2, 2) - pow(x1, 2)) * (x3 - x1)
					- (pow(x3, 2) - pow(x1, 2)) * (x2 - x1));

	float B = ((y2 - y1) * (pow(x3, 2) - pow(x1, 2))
			- (y3 - y1) * (pow(x2, 2) - pow(x1, 2)))
			/ ((x2 - x1) * (pow(x3, 2) - pow(x1, 2))
					- (x3 - x1) * (pow(x2, 2) - pow(x1, 2)));

	float C = y1 - A * pow(x1, 2) - B * x1;

	abc_vec.push_back(A);
	abc_vec.push_back(B);
	abc_vec.push_back(C);

	return abc_vec;
}

float dubins(vector<float> pinit, float thinit, vector<float> pf, float thf,
		float R) {

	vector<float> temp1;
	temp1.push_back(-999);
	temp1.push_back(-999);
	vector<float> temp2;
	temp2.push_back(-999);
	temp2.push_back(-999);
	vector<float> temp3;
	temp3.push_back(-999);
	temp3.push_back(-999);

	vector<float> vi;
	vi.push_back(cos(thinit));
	vi.push_back(sin(thinit));
	//cout << "vi = [" << vi[0] << "," << vi[1] << "]" << endl;

	vector<float> vf;
	vf.push_back(cos(thf));
	vf.push_back(sin(thf));
	//cout << "vf = [" << vf[0] << "," << vf[1] << "]" << endl;

	vector<float> initialCircleLeft_centre;
	initialCircleLeft_centre.push_back(pinit[0] - R * vi[1]); // minimum - radius circles
	initialCircleLeft_centre.push_back(pinit[1] + R * vi[0]);
	//cout << "initialCircleLeft_centre = [" << initialCircleLeft_centre[0] << "," << initialCircleLeft_centre[1] << "]" << endl;

	vector<float> initialCircleRight_centre;
	initialCircleRight_centre.push_back(pinit[0] + R * vi[1]);
	initialCircleRight_centre.push_back(pinit[1] - R * vi[0]);
	//cout << "initialCircleRight_centre = [" << initialCircleRight_centre[0] << "," << initialCircleRight_centre[1] << "]" << endl;

	vector<float> finalCircleLeft_centre;
	finalCircleLeft_centre.push_back(pf[0] - R * vf[1]);
	finalCircleLeft_centre.push_back(pf[1] + R * vf[0]);
	//cout << "finalCircleLeft_centre = [" << finalCircleLeft_centre[0] << "," << finalCircleLeft_centre[1] << "]" << endl;

	vector<float> finalCircleRight_centre;
	finalCircleRight_centre.push_back(pf[0] + R * vf[1]);
	finalCircleRight_centre.push_back(pf[1] - R * vf[0]);
	//cout << "finalCircleRight_centre = [" << finalCircleRight_centre[0] << "," << finalCircleRight_centre[1] << "]" << endl;

	float initialCircleLeft_radius = R;
	float initialCircleRight_radius = R;
	float finalCircleLeft_radius = R;
	float finalCircleRight_radius = R;
	//cout << "initialCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "initialCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;

	vector<float> initialCircleLeft_orientation;
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(1);
	//cout << "initialCircleLeft_orientation = [" << initialCircleLeft_orientation[0] << ","
	//											<< initialCircleLeft_orientation[1] << ","
	//											<< initialCircleLeft_orientation[2] << "]" << endl;

	vector<float> initialCircleRight_orientation;
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(-1);
	//cout << "initialCircleRight_orientation = ["	<< initialCircleRight_orientation[0] << ","
	//												<< initialCircleRight_orientation[1] << ","
	//												<< initialCircleRight_orientation[2] << "]" << endl;

	vector<float> finalCircleLeft_orientation;
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(1);
	//cout << "finalCircleLeft_orientation = ["	<< finalCircleLeft_orientation[0] << ","
	//											<< finalCircleLeft_orientation[1] << ","
	//											<< finalCircleLeft_orientation[2] << "]" << endl;

	vector<float> finalCircleRight_orientation;
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(-1);
	//cout << "finalCircleRight_orientation = ["	<< finalCircleRight_orientation[0] << ","
	//											<< finalCircleRight_orientation[1] << ","
	//											<< finalCircleRight_orientation[2] << "]" << endl;

	//AdmissiblePaths = cell(0);                                                  % Preallocate an empty cell array to have a consistent form for the commands that update it with new paths.

	// --------------------------------------------------------------------
	// START: RSR
	// --------------------------------------------------------------------

	float CSC_RR_length = 0;

	vector<float> v;
	v = vector_subtraction(finalCircleRight_centre, initialCircleRight_centre);
	CSC_RR_length = CSC_RR_length + vector_norm(v);

	vector<float> offset;
	offset.push_back(-R * v[1] / CSC_RR_length);
	offset.push_back(R * v[0] / CSC_RR_length);

	vector<float> initialTouchPoint = vector_adding(initialCircleRight_centre,
			offset); // point of contact of the initial circle with the common tangent
	vector<float> finalTouchPoint = vector_adding(finalCircleRight_centre,
			offset); // point of contact of the final circle with the common tangent
	//vector_print_float(2, initialTouchPoint);
	//vector_print_float(2, finalTouchPoint);

	// Find the initial arc

	vector<float> R_pinit = vector_subtraction(pinit,
			initialCircleRight_centre);
	float phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
	vector<float> R_contact = vector_subtraction(initialTouchPoint,
			initialCircleRight_centre);
	float phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

	float initialArcAngle;
	if (phi_pinit >= phi_contact)
		initialArcAngle = phi_pinit - phi_contact;
	else
		initialArcAngle = 2 * PI - (phi_contact - phi_pinit);

	CSC_RR_length = CSC_RR_length + R * initialArcAngle;

	// Find the final arc

	R_contact = vector_subtraction(finalTouchPoint, finalCircleRight_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
	vector<float> R_pf = vector_subtraction(pf, finalCircleRight_centre);
	float phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

	float finalArcAngle;
	if (phi_contact >= phi_pf)
		finalArcAngle = phi_contact - phi_pf;
	else
		finalArcAngle = 2 * PI - (phi_pf - phi_contact);

	CSC_RR_length = CSC_RR_length + R * finalArcAngle;

	float RSR_length = CSC_RR_length;
//		cout << "RSR_length = " << RSR_length << endl;

	// --------------------------------------------------------------------
	// END: RSR
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: RLR
	// --------------------------------------------------------------------

	float CCC_RR_length = 0;

	vector<float> d = vector_subtraction(pf, pinit);
	float phi = atan2_0_2pi(d[1], d[0]);

	float th1;
	if (thinit > phi)
		th1 = thinit - phi;
	else
		th1 = 2 * PI - (phi - thinit);

	float th2;
	if (thf > phi)
		th2 = thf - phi;
	else
		th2 = 2 * PI - (phi - thf);

	//float check1 = vector_norm(d);
	//float check2 = sqrt(4 - pow(fabs(cos(th1)) + fabs(cos(th2)), 2)) + fabs(sin(th1)) + abs(sin(th2));
	//float check3 = vector_norm(v);
	//float check4 = 4 * R;
	//float check5 = 2 * PI - acos((8 * pow(R, 2) - (vector_sum(vector_multiplication(v, v)))) / 8 / pow(R, 2));
	//float check6 = PI;
	//cout << "check1=" << check1 << endl;
	//cout << "check2=" << check2 << endl;
	//cout << "check3=" << check3 << endl;
	//cout << "check4=" << check4 << endl;
	//cout << "check5=" << check5 << endl;
	//cout << "check6=" << check6 << endl;

	if ((vector_norm(d)
			< sqrt(4 - pow(fabs(cos(th1)) + fabs(cos(th2)), 2)) + fabs(sin(th1))
					+ abs(sin(th2))) && (vector_norm(v) <= 4 * R)
			&& ((2 * PI
					- acos(
							(8 * pow(R, 2)
									- (vector_sum(vector_multiplication(v, v))))
									/ 8 / pow(R, 2))) > PI)) {

		// Location of the centre of the lower connecting circle

		vector<float> v_unit = vector_scalar_multiplication(
				(1 / vector_norm(v)), v);
		//vector_print_float(2, v_unit);
		temp1 = vector_scalar_multiplication(0.5, v);
		temp2[0] = v_unit[1];
		temp2[1] = -v_unit[0];
		temp3 = vector_scalar_multiplication(
				sqrt(
						4 * pow(R, 2)
								- vector_sum(vector_multiplication(v, v)) / 4),
				temp2);
		vector<float> connectingCircle_RR_centre = vector_adding(
				initialCircleRight_centre, vector_adding(temp1, temp3));
		//vector_print_float(2, connectingCircle_RR_centre);

		////////////////////////
		// Find the middle arc.
		////////////////////////
		vector<float> R_initial = vector_scalar_multiplication(0.5,
				(vector_subtraction(initialCircleRight_centre,
						connectingCircle_RR_centre)));
		float phi_initial = atan2_0_2pi(R_initial[1], R_initial[0]);
		vector<float> R_final = vector_scalar_multiplication(0.5,
				vector_subtraction(finalCircleRight_centre,
						connectingCircle_RR_centre));
		float phi_final = atan2_0_2pi(R_final[1], R_final[0]);

		float middleArcAngle;
		if (phi_initial >= phi_final)
			middleArcAngle = 2 * PI - (phi_initial - phi_final);
		else
			middleArcAngle = phi_final - phi_initial;

		CCC_RR_length = CCC_RR_length + R * middleArcAngle;
		//cout << "CCC_RR_length=" << CCC_RR_length << endl;

		////////////////////////
		// Find the initial arc.
		////////////////////////

		R_contact = vector_scalar_multiplication(0.5,
				vector_subtraction(connectingCircle_RR_centre,
						initialCircleRight_centre));
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		if (phi_pinit >= phi_contact)
			initialArcAngle = phi_pinit - phi_contact;
		else
			initialArcAngle = 2 * PI - (phi_contact - phi_pinit);

		if (initialArcAngle <= middleArcAngle) // Another necessary condition for optimality of CCC paths.
				{
			CCC_RR_length = CCC_RR_length + R * initialArcAngle;
			//cout << "CCC_RR_length=" << CCC_RR_length << endl;

			////////////////////////
			// Find the final arc.
			////////////////////////

			R_contact = vector_scalar_multiplication(0.5,
					vector_subtraction(connectingCircle_RR_centre,
							finalCircleRight_centre));
			phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
			R_pf = vector_subtraction(pf, finalCircleRight_centre);
			phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

			if (phi_contact >= phi_pf)
				finalArcAngle = phi_contact - phi_pf;
			else
				finalArcAngle = 2 * PI - (phi_pf - phi_contact);

			if ((finalArcAngle <= middleArcAngle)
					&& ((finalArcAngle < (middleArcAngle - PI))
							|| (initialArcAngle < (middleArcAngle - PI)))) // Necessary conditions for optimality of CCC paths.
					{
				CCC_RR_length = CCC_RR_length + R * finalArcAngle;
				//cout << "CCC_RR_length=" << CCC_RR_length << endl;
			} else {
				CCC_RR_length = -999;
			}
		} else {
			CCC_RR_length = -999;
		}
	} else {
		CCC_RR_length = -999;
	}

	float RLR_length = CCC_RR_length;
//		cout << "RLR_length = " << RLR_length << endl;

	// --------------------------------------------------------------------
	// END: RLR
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: RSL
	// --------------------------------------------------------------------

	float CSC_RL_length = 0;

	v = vector_subtraction(finalCircleLeft_centre, initialCircleRight_centre);
	float v_angle = atan2_0_2pi(v[1], v[0]);
	if (vector_norm(v) >= 2 * R) {
		CSC_RL_length = CSC_RL_length
				+ sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2));
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;

		phi = v_angle
				+ atan2_0_2pi(
						sqrt(
								vector_sum(vector_multiplication(v, v))
										- 4 * pow(R, 2)), 2 * R);

		temp1[0] = R * cos(phi);
		temp1[1] = R * sin(phi);
		vector<float> CSC_RL_S_initialpoint = vector_adding(
				initialCircleRight_centre, temp1); // point of contact of the initial circle with the common tangent
		//vector_print_float(2, CSC_RL_S_initialpoint);
		temp2[0] = sin(phi);
		temp2[1] = -cos(phi);
		temp3 = vector_scalar_multiplication(
				sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2)),
				temp2);
		vector<float> CSC_RL_S_finalpoint = vector_adding(CSC_RL_S_initialpoint,
				temp3); // point of contact of the final circle with the common tangent
		//vector_print_float(2, CSC_RL_S_finalpoint);

		// Find the initial arc

		R_pinit = vector_subtraction(pinit, initialCircleRight_centre);
		phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
		R_contact = vector_subtraction(CSC_RL_S_initialpoint,
				initialCircleRight_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		if (phi_pinit >= phi_contact)
			initialArcAngle = phi_pinit - phi_contact;
		else
			initialArcAngle = 2 * PI - (phi_contact - phi_pinit);

		CSC_RL_length = CSC_RL_length + R * initialArcAngle;
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;

		// Find the final arc

		R_contact = vector_subtraction(CSC_RL_S_finalpoint,
				finalCircleLeft_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
		R_pf = vector_subtraction(pf, finalCircleLeft_centre);
		phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

		if (phi_contact >= phi_pf)
			finalArcAngle = 2 * PI - (phi_contact - phi_pf);
		else
			finalArcAngle = phi_pf - phi_contact;

		CSC_RL_length = CSC_RL_length + R * finalArcAngle;
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;
	} else {
		CSC_RL_length = -999;
	}

	float RSL_length = CSC_RL_length;
//		cout << "RSL_length = " << RSL_length << endl;

	// --------------------------------------------------------------------
	// END: RSL
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: LSL
	// --------------------------------------------------------------------

	float CSC_LL_length = 0;

	v = vector_subtraction(finalCircleLeft_centre, initialCircleLeft_centre);

	CSC_LL_length = CSC_LL_length + vector_norm(v);
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	offset[0] = R * v[1] / CSC_LL_length;
	offset[1] = -R * v[0] / CSC_LL_length;
	initialTouchPoint = vector_adding(initialCircleLeft_centre, offset); // point of contact of the initial circle with the common tangent
	finalTouchPoint = vector_adding(finalCircleLeft_centre, offset); // point of contact of the final circle with the common tangent

	// Find the initial arc

	R_pinit = vector_subtraction(pinit, initialCircleLeft_centre);
	phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
	R_contact = vector_subtraction(initialTouchPoint, initialCircleLeft_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

	if (phi_pinit >= phi_contact)
		initialArcAngle = 2 * PI - (phi_pinit - phi_contact);
	else
		initialArcAngle = phi_contact - phi_pinit;

	CSC_LL_length = CSC_LL_length + R * initialArcAngle;
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	// Find the final arc

	R_contact = vector_subtraction(finalTouchPoint, finalCircleLeft_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
	R_pf = vector_subtraction(pf, finalCircleLeft_centre);
	phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

	if (phi_contact >= phi_pf)
		finalArcAngle = 2 * PI - (phi_contact - phi_pf);
	else
		finalArcAngle = phi_pf - phi_contact;

	CSC_LL_length = CSC_LL_length + R * finalArcAngle;
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	float LSL_length = CSC_LL_length;
//		cout << "LSL_length = " << LSL_length << endl;

	// --------------------------------------------------------------------
	// END: LSL
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: LRL
	// --------------------------------------------------------------------

	float CCC_LL_length = 0;

	d = vector_subtraction(pf, pinit);
	phi = atan2_0_2pi(d[1], d[0]);

	if (thinit > phi)
		th1 = thinit - phi;
	else
		th1 = 2 * PI - (phi - thinit);

	if (thf > phi)
		th2 = thf - phi;
	else
		th2 = 2 * PI - (phi - thf);

	v = vector_subtraction(finalCircleLeft_centre, initialCircleLeft_centre);

	//float check1 = vector_norm(d);
	//float check2 = sqrt(4 - pow(fabs(cos(th1)) + fabs(cos(th2)), 2)) + fabs(sin(th1)) + abs(sin(th2));
	//float check3 = vector_norm(v);
	//float check4 = 4 * R;
	//float check5 = 2 * PI - acos((8 * pow(R, 2) - (vector_sum(vector_multiplication(v, v)))) / 8 / pow(R, 2));
	//float check6 = PI;
	//cout << "check1=" << check1 << endl;
	//cout << "check2=" << check2 << endl;
	//cout << "check3=" << check3 << endl;
	//cout << "check4=" << check4 << endl;
	//cout << "check5=" << check5 << endl;
	//cout << "check6=" << check6 << endl;

	if ((vector_norm(d)
			< sqrt(4 - pow(fabs(cos(th1)) + fabs(cos(th2)), 2)) + fabs(sin(th1))
					+ abs(sin(th2))) && (vector_norm(v) <= 4 * R)
			&& ((2 * PI
					- acos(
							(8 * pow(R, 2)
									- (vector_sum(vector_multiplication(v, v))))
									/ 8 / pow(R, 2))) > PI)) {
		vector<float> v_unit = vector_scalar_multiplication(1 / vector_norm(v),
				v);
		//vector_print_float(2, v_unit);
		temp1 = vector_scalar_multiplication(0.5, v);
		temp2[0] = -v_unit[1];
		temp2[1] = v_unit[0];
		temp3 = vector_scalar_multiplication(
				sqrt(
						4 * pow(R, 2)
								- vector_sum(vector_multiplication(v, v)) / 4),
				temp2);
		vector<float> connectingCircle_LL_centre = vector_adding(
				initialCircleLeft_centre, vector_adding(temp1, temp3));
		//vector_print_float(2, connectingCircle_LL_centre);

		////////////////////////
		// Find the middle arc.
		////////////////////////

		vector<float> R_initial = vector_scalar_multiplication(0.5,
				vector_subtraction(initialCircleLeft_centre,
						connectingCircle_LL_centre));
		float phi_initial = atan2_0_2pi(R_initial[1], R_initial[0]);
		vector<float> R_final = vector_scalar_multiplication(0.5,
				vector_subtraction(finalCircleLeft_centre,
						connectingCircle_LL_centre));
		float phi_final = atan2_0_2pi(R_final[1], R_final[0]);

		float middleArcAngle;
		if (phi_initial >= phi_final)
			middleArcAngle = phi_initial - phi_final;
		else
			middleArcAngle = 2 * PI - (phi_final - phi_initial);

		CCC_LL_length = CCC_LL_length + R * middleArcAngle;
		//cout << "CCC_LL_length=" << CCC_LL_length << endl;

		////////////////////////
		// Find the initial arc.
		////////////////////////

		R_contact = vector_scalar_multiplication(0.5,
				vector_subtraction(connectingCircle_LL_centre,
						initialCircleLeft_centre));
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		if (phi_pinit >= phi_contact)
			initialArcAngle = 2 * PI - (phi_pinit - phi_contact);
		else
			initialArcAngle = phi_contact - phi_pinit;

		if (initialArcAngle <= middleArcAngle) // Another necessary condition for optimality of CCC paths.
				{
			CCC_LL_length = CCC_LL_length + R * initialArcAngle;
			//cout << "CCC_LL_length=" << CCC_LL_length << endl;

			////////////////////////
			// Find the final arc.
			////////////////////////

			R_contact = vector_scalar_multiplication(0.5,
					vector_subtraction(connectingCircle_LL_centre,
							finalCircleLeft_centre));
			phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
			R_pf = vector_subtraction(pf, finalCircleLeft_centre);
			phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

			if (phi_contact >= phi_pf)
				finalArcAngle = 2 * PI - (phi_contact - phi_pf);
			else
				finalArcAngle = phi_pf - phi_contact;

			if ((finalArcAngle <= middleArcAngle)
					&& ((finalArcAngle < (middleArcAngle - PI))
							|| (initialArcAngle < (middleArcAngle - PI)))) // Necessary conditions for optimality of CCC paths.
					{
				CCC_LL_length = CCC_LL_length + R * finalArcAngle;
				//cout << "CCC_LL_length=" << CCC_LL_length << endl;
			} else {
				CCC_LL_length = -999;
			}
		} else {
			CCC_LL_length = -999;
		}
	} else {
		CCC_LL_length = -999;
	}

	float LRL_length = CCC_LL_length;
//		cout << "LRL_length = " << LRL_length << endl;

	// --------------------------------------------------------------------
	// END: LRL
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: LSR
	// --------------------------------------------------------------------

	float CSC_LR_length = 0;

	v = vector_subtraction(finalCircleRight_centre, initialCircleLeft_centre);
	v_angle = atan2_0_2pi(v[1], v[0]);
	if (vector_norm(v) >= 2 * R) {
		CSC_LR_length = CSC_LR_length
				+ sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2));
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;

		phi = v_angle
				- atan2_0_2pi(
						sqrt(
								vector_sum(vector_multiplication(v, v))
										- 4 * pow(R, 2)), 2 * R);

		temp1[0] = R * cos(phi);
		temp1[1] = R * sin(phi);
		vector<float> CSC_LR_S_initialpoint = vector_adding(
				initialCircleLeft_centre, temp1); // point of contact of the initial circle with the common tangent
		temp2[0] = -sin(phi);
		temp2[1] = cos(phi);
		temp3 = vector_scalar_multiplication(
				sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2)),
				temp2);
		vector<float> CSC_LR_S_finalpoint = vector_adding(CSC_LR_S_initialpoint,
				temp3); // point of contact of the final circle with the common tangent
		//vector_print_float(2, CSC_LR_S_initialpoint);
		//vector_print_float(2, CSC_LR_S_finalpoint);

		// Find the initial arc
		R_pinit = vector_subtraction(pinit, initialCircleLeft_centre);
		phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
		R_contact = vector_subtraction(CSC_LR_S_initialpoint,
				initialCircleLeft_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		if (phi_pinit >= phi_contact)
			initialArcAngle = 2 * PI - (phi_pinit - phi_contact);
		else
			initialArcAngle = phi_contact - phi_pinit;

		CSC_LR_length = CSC_LR_length + R * initialArcAngle;
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;

		// Find the final arc

		R_contact = vector_subtraction(CSC_LR_S_finalpoint,
				finalCircleRight_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
		R_pf = vector_subtraction(pf, finalCircleRight_centre);
		phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

		if (phi_contact >= phi_pf)
			finalArcAngle = phi_contact - phi_pf;
		else
			finalArcAngle = 2 * PI - (phi_pf - phi_contact);

		CSC_LR_length = CSC_LR_length + R * finalArcAngle;
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;
	} else {
		CSC_LR_length = -999;
	}

	float LSR_length = CSC_LR_length;
//		cout << "LSR_length = " << LSR_length << endl;

	// --------------------------------------------------------------------
	// END: LSR
	// --------------------------------------------------------------------

	//RSR_length - 1
	//RLR_length - 2
	//RSL_length - 3
	//LSL_length - 4
	//LRL_length - 5
	//LSR_length - 6

	float best_solution = RSR_length;
	float solution_mark = 1;
	if (((best_solution > RLR_length) && (RLR_length != -999))
			|| (best_solution == -999)) {
		best_solution = RLR_length;
		solution_mark = 2;
	}
	if (((best_solution > RSL_length) && (RSL_length != -999))
			|| (best_solution == -999)) {
		best_solution = RSL_length;
		solution_mark = 3;
	}
	if (((best_solution > LSL_length) && (LSL_length != -999))
			|| (best_solution == -999)) {
		best_solution = LSL_length;
		solution_mark = 4;
	}
	if (((best_solution > LRL_length) && (LRL_length != -999))
			|| (best_solution == -999)) {
		best_solution = LRL_length;
		solution_mark = 5;
	}
	if (((best_solution > LSR_length) && (LSR_length != -999))
			|| (best_solution == -999)) {
		best_solution = LSR_length;
		solution_mark = 6;
	}

	return solution_mark;

}

float dubins_four(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R) {

	vector<float> temp1;
	temp1.push_back(-999);
	temp1.push_back(-999);
	vector<float> temp2;
	temp2.push_back(-999);
	temp2.push_back(-999);
	vector<float> temp3;
	temp3.push_back(-999);
	temp3.push_back(-999);

	vector<float> vi;
	vi.push_back(cos(thinit));
	vi.push_back(sin(thinit));
	//cout << "vi = [" << vi[0] << "," << vi[1] << "]" << endl;

	vector<float> vf;
	vf.push_back(cos(thf));
	vf.push_back(sin(thf));
	//cout << "vf = [" << vf[0] << "," << vf[1] << "]" << endl;

	vector<float> initialCircleLeft_centre;
	initialCircleLeft_centre.push_back(pinit[0] - R * vi[1]); // minimum - radius circles
	initialCircleLeft_centre.push_back(pinit[1] + R * vi[0]);
	//cout << "initialCircleLeft_centre = [" << initialCircleLeft_centre[0] << "," << initialCircleLeft_centre[1] << "]" << endl;

	vector<float> initialCircleRight_centre;
	initialCircleRight_centre.push_back(pinit[0] + R * vi[1]);
	initialCircleRight_centre.push_back(pinit[1] - R * vi[0]);
	//cout << "initialCircleRight_centre = [" << initialCircleRight_centre[0] << "," << initialCircleRight_centre[1] << "]" << endl;

	vector<float> finalCircleLeft_centre;
	finalCircleLeft_centre.push_back(pf[0] - R * vf[1]);
	finalCircleLeft_centre.push_back(pf[1] + R * vf[0]);
	//cout << "finalCircleLeft_centre = [" << finalCircleLeft_centre[0] << "," << finalCircleLeft_centre[1] << "]" << endl;

	vector<float> finalCircleRight_centre;
	finalCircleRight_centre.push_back(pf[0] + R * vf[1]);
	finalCircleRight_centre.push_back(pf[1] - R * vf[0]);
	//cout << "finalCircleRight_centre = [" << finalCircleRight_centre[0] << "," << finalCircleRight_centre[1] << "]" << endl;

	float initialCircleLeft_radius = R;
	float initialCircleRight_radius = R;
	float finalCircleLeft_radius = R;
	float finalCircleRight_radius = R;
	//cout << "initialCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "initialCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;

	vector<float> initialCircleLeft_orientation;
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(1);
	//cout << "initialCircleLeft_orientation = [" << initialCircleLeft_orientation[0] << ","
	//											<< initialCircleLeft_orientation[1] << ","
	//											<< initialCircleLeft_orientation[2] << "]" << endl;

	vector<float> initialCircleRight_orientation;
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(-1);
	//cout << "initialCircleRight_orientation = ["	<< initialCircleRight_orientation[0] << ","
	//												<< initialCircleRight_orientation[1] << ","
	//												<< initialCircleRight_orientation[2] << "]" << endl;

	vector<float> finalCircleLeft_orientation;
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(1);
	//cout << "finalCircleLeft_orientation = ["	<< finalCircleLeft_orientation[0] << ","
	//											<< finalCircleLeft_orientation[1] << ","
	//											<< finalCircleLeft_orientation[2] << "]" << endl;

	vector<float> finalCircleRight_orientation;
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(-1);
	//cout << "finalCircleRight_orientation = ["	<< finalCircleRight_orientation[0] << ","
	//											<< finalCircleRight_orientation[1] << ","
	//											<< finalCircleRight_orientation[2] << "]" << endl;

	//AdmissiblePaths = cell(0);                                                  % Preallocate an empty cell array to have a consistent form for the commands that update it with new paths.

	// --------------------------------------------------------------------
	// START: RSR
	// --------------------------------------------------------------------

	float CSC_RR_length = 0;

	vector<float> v;
	v = vector_subtraction(finalCircleRight_centre, initialCircleRight_centre);
	CSC_RR_length = CSC_RR_length + vector_norm(v);

	vector<float> offset;
	offset.push_back(-R * v[1] / CSC_RR_length);
	offset.push_back(R * v[0] / CSC_RR_length);

	vector<float> initialTouchPoint = vector_adding(initialCircleRight_centre,
			offset); // point of contact of the initial circle with the common tangent
	vector<float> finalTouchPoint = vector_adding(finalCircleRight_centre,
			offset); // point of contact of the final circle with the common tangent
	//vector_print_float(2, initialTouchPoint);
	//vector_print_float(2, finalTouchPoint);

	// Find the initial arc

	vector<float> R_pinit = vector_subtraction(pinit,
			initialCircleRight_centre);
	float phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
	vector<float> R_contact = vector_subtraction(initialTouchPoint,
			initialCircleRight_centre);
	float phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

	float initialArcAngle;
	if (phi_pinit >= phi_contact)
		initialArcAngle = phi_pinit - phi_contact;
	else
		initialArcAngle = 2 * PI - (phi_contact - phi_pinit);

	CSC_RR_length = CSC_RR_length + R * initialArcAngle;

	// Find the final arc

	R_contact = vector_subtraction(finalTouchPoint, finalCircleRight_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
	vector<float> R_pf = vector_subtraction(pf, finalCircleRight_centre);
	float phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

	float finalArcAngle;
	if (phi_contact >= phi_pf)
		finalArcAngle = phi_contact - phi_pf;
	else
		finalArcAngle = 2 * PI - (phi_pf - phi_contact);

	CSC_RR_length = CSC_RR_length + R * finalArcAngle;

	float RSR_length = CSC_RR_length;
	cout << "RSR_length = " << RSR_length << endl;

	// --------------------------------------------------------------------
	// END: RSR
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: RSL
	// --------------------------------------------------------------------

	float CSC_RL_length = 0;

	v = vector_subtraction(finalCircleLeft_centre, initialCircleRight_centre);
	float v_angle = atan2_0_2pi(v[1], v[0]);
	if (vector_norm(v) >= 2 * R) {
		CSC_RL_length = CSC_RL_length
				+ sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2));
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;

		float phi = v_angle
				+ atan2_0_2pi(
						sqrt(
								vector_sum(vector_multiplication(v, v))
										- 4 * pow(R, 2)), 2 * R);

		temp1[0] = R * cos(phi);
		temp1[1] = R * sin(phi);
		vector<float> CSC_RL_S_initialpoint = vector_adding(
				initialCircleRight_centre, temp1); // point of contact of the initial circle with the common tangent
		//vector_print_float(2, CSC_RL_S_initialpoint);
		temp2[0] = sin(phi);
		temp2[1] = -cos(phi);
		temp3 = vector_scalar_multiplication(
				sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2)),
				temp2);
		vector<float> CSC_RL_S_finalpoint = vector_adding(CSC_RL_S_initialpoint,
				temp3); // point of contact of the final circle with the common tangent
		//vector_print_float(2, CSC_RL_S_finalpoint);

		// Find the initial arc

		R_pinit = vector_subtraction(pinit, initialCircleRight_centre);
		phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
		R_contact = vector_subtraction(CSC_RL_S_initialpoint,
				initialCircleRight_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		if (phi_pinit >= phi_contact)
			initialArcAngle = phi_pinit - phi_contact;
		else
			initialArcAngle = 2 * PI - (phi_contact - phi_pinit);

		CSC_RL_length = CSC_RL_length + R * initialArcAngle;
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;

		// Find the final arc

		R_contact = vector_subtraction(CSC_RL_S_finalpoint,
				finalCircleLeft_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
		R_pf = vector_subtraction(pf, finalCircleLeft_centre);
		phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

		if (phi_contact >= phi_pf)
			finalArcAngle = 2 * PI - (phi_contact - phi_pf);
		else
			finalArcAngle = phi_pf - phi_contact;

		CSC_RL_length = CSC_RL_length + R * finalArcAngle;
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;
	} else {
		CSC_RL_length = -999;
	}

	float RSL_length = CSC_RL_length;
	cout << "RSL_length = " << RSL_length << endl;

	// --------------------------------------------------------------------
	// END: RSL
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: LSL
	// --------------------------------------------------------------------

	float CSC_LL_length = 0;

	v = vector_subtraction(finalCircleLeft_centre, initialCircleLeft_centre);

	CSC_LL_length = CSC_LL_length + vector_norm(v);
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	offset[0] = R * v[1] / CSC_LL_length;
	offset[1] = -R * v[0] / CSC_LL_length;
	initialTouchPoint = vector_adding(initialCircleLeft_centre, offset); // point of contact of the initial circle with the common tangent
	finalTouchPoint = vector_adding(finalCircleLeft_centre, offset); // point of contact of the final circle with the common tangent

	// Find the initial arc

	R_pinit = vector_subtraction(pinit, initialCircleLeft_centre);
	phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
	R_contact = vector_subtraction(initialTouchPoint, initialCircleLeft_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

	if (phi_pinit >= phi_contact)
		initialArcAngle = 2 * PI - (phi_pinit - phi_contact);
	else
		initialArcAngle = phi_contact - phi_pinit;

	CSC_LL_length = CSC_LL_length + R * initialArcAngle;
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	// Find the final arc

	R_contact = vector_subtraction(finalTouchPoint, finalCircleLeft_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
	R_pf = vector_subtraction(pf, finalCircleLeft_centre);
	phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

	if (phi_contact >= phi_pf)
		finalArcAngle = 2 * PI - (phi_contact - phi_pf);
	else
		finalArcAngle = phi_pf - phi_contact;

	CSC_LL_length = CSC_LL_length + R * finalArcAngle;
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	float LSL_length = CSC_LL_length;
	cout << "LSL_length = " << LSL_length << endl;

	// --------------------------------------------------------------------
	// END: LSL
	// --------------------------------------------------------------------

	// --------------------------------------------------------------------
	// START: LSR
	// --------------------------------------------------------------------

	float CSC_LR_length = 0;

	v = vector_subtraction(finalCircleRight_centre, initialCircleLeft_centre);
	v_angle = atan2_0_2pi(v[1], v[0]);
	if (vector_norm(v) >= 2 * R) {
		CSC_LR_length = CSC_LR_length
				+ sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2));
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;

		float phi = v_angle
				- atan2_0_2pi(
						sqrt(
								vector_sum(vector_multiplication(v, v))
										- 4 * pow(R, 2)), 2 * R);

		temp1[0] = R * cos(phi);
		temp1[1] = R * sin(phi);
		vector<float> CSC_LR_S_initialpoint = vector_adding(
				initialCircleLeft_centre, temp1); // point of contact of the initial circle with the common tangent
		temp2[0] = -sin(phi);
		temp2[1] = cos(phi);
		temp3 = vector_scalar_multiplication(
				sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2)),
				temp2);
		vector<float> CSC_LR_S_finalpoint = vector_adding(CSC_LR_S_initialpoint,
				temp3); // point of contact of the final circle with the common tangent
		//vector_print_float(2, CSC_LR_S_initialpoint);
		//vector_print_float(2, CSC_LR_S_finalpoint);

		// Find the initial arc
		R_pinit = vector_subtraction(pinit, initialCircleLeft_centre);
		phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
		R_contact = vector_subtraction(CSC_LR_S_initialpoint,
				initialCircleLeft_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		if (phi_pinit >= phi_contact)
			initialArcAngle = 2 * PI - (phi_pinit - phi_contact);
		else
			initialArcAngle = phi_contact - phi_pinit;

		CSC_LR_length = CSC_LR_length + R * initialArcAngle;
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;

		// Find the final arc

		R_contact = vector_subtraction(CSC_LR_S_finalpoint,
				finalCircleRight_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
		R_pf = vector_subtraction(pf, finalCircleRight_centre);
		phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

		if (phi_contact >= phi_pf)
			finalArcAngle = phi_contact - phi_pf;
		else
			finalArcAngle = 2 * PI - (phi_pf - phi_contact);

		CSC_LR_length = CSC_LR_length + R * finalArcAngle;
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;
	} else {
		CSC_LR_length = -999;
	}

	float LSR_length = CSC_LR_length;
	cout << "LSR_length = " << LSR_length << endl;

	// --------------------------------------------------------------------
	// END: LSR
	// --------------------------------------------------------------------

	//RSR_length - 1
	//RSL_length - 3
	//LSL_length - 4
	//LSR_length - 6

	float best_solution = RSR_length;
	if (((best_solution > RSL_length) && (RSL_length != -999))
			|| (best_solution == -999)) {
		best_solution = RSL_length;
	}
	if (((best_solution > LSL_length) && (LSL_length != -999))
			|| (best_solution == -999)) {
		best_solution = LSL_length;
	}
	if (((best_solution > LSR_length) && (LSR_length != -999))
			|| (best_solution == -999)) {
		best_solution = LSR_length;
	}

	return best_solution;

}

vector<float> dubins_LSL(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R) {
	vector<float> result_vector;

	vector<float> temp1;
	temp1.push_back(-999);
	temp1.push_back(-999);
	vector<float> temp2;
	temp2.push_back(-999);
	temp2.push_back(-999);
	vector<float> temp3;
	temp3.push_back(-999);
	temp3.push_back(-999);

	vector<float> vi;
	vi.push_back(cos(thinit));
	vi.push_back(sin(thinit));
	//cout << "vi = [" << vi[0] << "," << vi[1] << "]" << endl;

	vector<float> vf;
	vf.push_back(cos(thf));
	vf.push_back(sin(thf));
	//cout << "vf = [" << vf[0] << "," << vf[1] << "]" << endl;

	vector<float> initialCircleLeft_centre;
	initialCircleLeft_centre.push_back(pinit[0] - R * vi[1]); // minimum - radius circles
	initialCircleLeft_centre.push_back(pinit[1] + R * vi[0]);
	//cout << "initialCircleLeft_centre = [" << initialCircleLeft_centre[0] << "," << initialCircleLeft_centre[1] << "]" << endl;

	vector<float> finalCircleLeft_centre;
	finalCircleLeft_centre.push_back(pf[0] - R * vf[1]);
	finalCircleLeft_centre.push_back(pf[1] + R * vf[0]);
	//cout << "finalCircleLeft_centre = [" << finalCircleLeft_centre[0] << "," << finalCircleLeft_centre[1] << "]" << endl;

	float initialCircleLeft_radius = R;
	float finalCircleLeft_radius = R;
	//cout << "initialCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "initialCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;

	vector<float> initialCircleLeft_orientation;
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(1);
	//cout << "initialCircleLeft_orientation = [" << initialCircleLeft_orientation[0] << ","
	//											<< initialCircleLeft_orientation[1] << ","
	//											<< initialCircleLeft_orientation[2] << "]" << endl;

	vector<float> finalCircleLeft_orientation;
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(1);
	//cout << "finalCircleLeft_orientation = ["	<< finalCircleLeft_orientation[0] << ","
	//											<< finalCircleLeft_orientation[1] << ","
	//											<< finalCircleLeft_orientation[2] << "]" << endl;

	// --------------------------------------------------------------------
	// START: LSL
	// --------------------------------------------------------------------

	float CSC_LL_length = 0;

	vector<float> v = vector_subtraction(finalCircleLeft_centre,
			initialCircleLeft_centre);

	float v_norm = vector_norm(v);
	CSC_LL_length = CSC_LL_length + v_norm;
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	vector<float> offset;
	offset.push_back(R * v[1] / CSC_LL_length);
	offset.push_back(-R * v[0] / CSC_LL_length);
	vector<float> initialTouchPoint = vector_adding(initialCircleLeft_centre,
			offset); // point of contact of the initial circle with the common tangent
	vector<float> finalTouchPoint = vector_adding(finalCircleLeft_centre,
			offset); // point of contact of the final circle with the common tangent

	// Find the initial arc

	vector<float> R_pinit = vector_subtraction(pinit, initialCircleLeft_centre);
	float phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
	vector<float> R_contact = vector_subtraction(initialTouchPoint,
			initialCircleLeft_centre);
	float phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

	float initialArcAngle;
	if (phi_pinit >= phi_contact)
		initialArcAngle = 2 * PI - (phi_pinit - phi_contact);
	else
		initialArcAngle = phi_contact - phi_pinit;

	CSC_LL_length = CSC_LL_length + R * initialArcAngle;
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	// Find the final arc

	R_contact = vector_subtraction(finalTouchPoint, finalCircleLeft_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
	vector<float> R_pf = vector_subtraction(pf, finalCircleLeft_centre);
	float phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

	float finalArcAngle;
	if (phi_contact >= phi_pf)
		finalArcAngle = 2 * PI - (phi_contact - phi_pf);
	else
		finalArcAngle = phi_pf - phi_contact;

	CSC_LL_length = CSC_LL_length + R * finalArcAngle;
	//cout << "CSC_LL_length=" << CSC_LL_length << endl;

	float LSL_length = CSC_LL_length;
	//cout << "LSL_length = " << LSL_length << endl;

	result_vector.push_back(LSL_length);
	result_vector.push_back(initialCircleLeft_centre[0]);
	result_vector.push_back(initialCircleLeft_centre[1]);
	result_vector.push_back(initialArcAngle);
	result_vector.push_back(phi_pinit);
	result_vector.push_back(v_norm);
	result_vector.push_back(finalCircleLeft_centre[0]);
	result_vector.push_back(finalCircleLeft_centre[1]);
	result_vector.push_back(finalArcAngle);
	result_vector.push_back(phi_contact);

	return result_vector;
	// --------------------------------------------------------------------
	// END: LSL
	// --------------------------------------------------------------------

}

vector<float> dubins_RSR(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R) {
	vector<float> result_vector;

	vector<float> temp1;
	temp1.push_back(-999);
	temp1.push_back(-999);
	vector<float> temp2;
	temp2.push_back(-999);
	temp2.push_back(-999);
	vector<float> temp3;
	temp3.push_back(-999);
	temp3.push_back(-999);

	vector<float> vi;
	vi.push_back(cos(thinit));
	vi.push_back(sin(thinit));
	//cout << "vi = [" << vi[0] << "," << vi[1] << "]" << endl;

	vector<float> vf;
	vf.push_back(cos(thf));
	vf.push_back(sin(thf));
	//cout << "vf = [" << vf[0] << "," << vf[1] << "]" << endl;

	vector<float> initialCircleRight_centre;
	initialCircleRight_centre.push_back(pinit[0] + R * vi[1]);
	initialCircleRight_centre.push_back(pinit[1] - R * vi[0]);
	//cout << "initialCircleRight_centre = [" << initialCircleRight_centre[0] << "," << initialCircleRight_centre[1] << "]" << endl;

	vector<float> finalCircleRight_centre;
	finalCircleRight_centre.push_back(pf[0] + R * vf[1]);
	finalCircleRight_centre.push_back(pf[1] - R * vf[0]);
	//cout << "finalCircleRight_centre = [" << finalCircleRight_centre[0] << "," << finalCircleRight_centre[1] << "]" << endl;

	float initialCircleRight_radius = R;
	float finalCircleRight_radius = R;
	//cout << "initialCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "initialCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;

	vector<float> initialCircleRight_orientation;
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(-1);
	//cout << "initialCircleRight_orientation = ["	<< initialCircleRight_orientation[0] << ","
	//												<< initialCircleRight_orientation[1] << ","
	//												<< initialCircleRight_orientation[2] << "]" << endl;

	vector<float> finalCircleRight_orientation;
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(-1);
	//cout << "finalCircleRight_orientation = ["	<< finalCircleRight_orientation[0] << ","
	//											<< finalCircleRight_orientation[1] << ","
	//											<< finalCircleRight_orientation[2] << "]" << endl;

	// --------------------------------------------------------------------
	// START: RSR
	// --------------------------------------------------------------------

	float CSC_RR_length = 0;

	vector<float> v;
	v = vector_subtraction(finalCircleRight_centre, initialCircleRight_centre);
	float v_norm = vector_norm(v);
	CSC_RR_length = CSC_RR_length + v_norm;

	vector<float> offset;
	offset.push_back(-R * v[1] / CSC_RR_length);
	offset.push_back(R * v[0] / CSC_RR_length);

	vector<float> initialTouchPoint = vector_adding(initialCircleRight_centre,
			offset); // point of contact of the initial circle with the common tangent
	vector<float> finalTouchPoint = vector_adding(finalCircleRight_centre,
			offset); // point of contact of the final circle with the common tangent
	//vector_print_float(2, initialTouchPoint);
	//vector_print_float(2, finalTouchPoint);

	// Find the initial arc

	vector<float> R_pinit = vector_subtraction(pinit,
			initialCircleRight_centre);
	float phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
	vector<float> R_contact = vector_subtraction(initialTouchPoint,
			initialCircleRight_centre);
	float phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

	float initialArcAngle;
	if (phi_pinit >= phi_contact)
		initialArcAngle = phi_pinit - phi_contact;
	else
		initialArcAngle = 2 * PI - (phi_contact - phi_pinit);

	CSC_RR_length = CSC_RR_length + R * initialArcAngle;

	// Find the final arc

	R_contact = vector_subtraction(finalTouchPoint, finalCircleRight_centre);
	phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
	vector<float> R_pf = vector_subtraction(pf, finalCircleRight_centre);
	float phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

	float finalArcAngle;
	if (phi_contact >= phi_pf)
		finalArcAngle = phi_contact - phi_pf;
	else
		finalArcAngle = 2 * PI - (phi_pf - phi_contact);

	CSC_RR_length = CSC_RR_length + R * finalArcAngle;

	float RSR_length = CSC_RR_length;
	//cout << "RSR_length = " << RSR_length << endl;

	// --------------------------------------------------------------------
	// END: RSR
	// --------------------------------------------------------------------

	result_vector.push_back(RSR_length);
	result_vector.push_back(initialCircleRight_centre[0]);
	result_vector.push_back(initialCircleRight_centre[1]);
	result_vector.push_back(initialArcAngle);
	result_vector.push_back(phi_pinit);
	result_vector.push_back(v_norm);
	result_vector.push_back(finalCircleRight_centre[0]);
	result_vector.push_back(finalCircleRight_centre[1]);
	result_vector.push_back(finalArcAngle);
	result_vector.push_back(phi_contact);

	return result_vector;
}

vector<float> dubins_LSR(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R) {
	vector<float> result_vector;

	vector<float> temp1;
	temp1.push_back(-999);
	temp1.push_back(-999);
	vector<float> temp2;
	temp2.push_back(-999);
	temp2.push_back(-999);
	vector<float> temp3;
	temp3.push_back(-999);
	temp3.push_back(-999);

	vector<float> vi;
	vi.push_back(cos(thinit));
	vi.push_back(sin(thinit));
	//cout << "vi = [" << vi[0] << "," << vi[1] << "]" << endl;

	vector<float> vf;
	vf.push_back(cos(thf));
	vf.push_back(sin(thf));
	//cout << "vf = [" << vf[0] << "," << vf[1] << "]" << endl;

	vector<float> initialCircleLeft_centre;
	initialCircleLeft_centre.push_back(pinit[0] - R * vi[1]); // minimum - radius circles
	initialCircleLeft_centre.push_back(pinit[1] + R * vi[0]);
	//cout << "initialCircleLeft_centre = [" << initialCircleLeft_centre[0] << "," << initialCircleLeft_centre[1] << "]" << endl;

	vector<float> finalCircleRight_centre;
	finalCircleRight_centre.push_back(pf[0] + R * vf[1]);
	finalCircleRight_centre.push_back(pf[1] - R * vf[0]);
	//cout << "finalCircleRight_centre = [" << finalCircleRight_centre[0] << "," << finalCircleRight_centre[1] << "]" << endl;

	float initialCircleLeft_radius = R;
	float finalCircleRight_radius = R;
	//cout << "initialCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "initialCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;

	vector<float> initialCircleLeft_orientation;
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(0);
	initialCircleLeft_orientation.push_back(1);
	//cout << "initialCircleLeft_orientation = [" << initialCircleLeft_orientation[0] << ","
	//											<< initialCircleLeft_orientation[1] << ","
	//											<< initialCircleLeft_orientation[2] << "]" << endl;

	vector<float> finalCircleRight_orientation;
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(0);
	finalCircleRight_orientation.push_back(-1);
	//cout << "finalCircleRight_orientation = ["	<< finalCircleRight_orientation[0] << ","
	//											<< finalCircleRight_orientation[1] << ","
	//											<< finalCircleRight_orientation[2] << "]" << endl;

	// --------------------------------------------------------------------
	// START: LSR
	// --------------------------------------------------------------------

	float CSC_LR_length = 0;

	vector<float> v = vector_subtraction(finalCircleRight_centre,
			initialCircleLeft_centre);
	float v_angle = atan2_0_2pi(v[1], v[0]);
	if (vector_norm(v) >= 2 * R) {
		float v_norm = sqrt(
				vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2));
		CSC_LR_length = CSC_LR_length + v_norm;
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;

		float phi = v_angle
				- atan2_0_2pi(
						sqrt(
								vector_sum(vector_multiplication(v, v))
										- 4 * pow(R, 2)), 2 * R);

		temp1[0] = R * cos(phi);
		temp1[1] = R * sin(phi);
		vector<float> CSC_LR_S_initialpoint = vector_adding(
				initialCircleLeft_centre, temp1); // point of contact of the initial circle with the common tangent
		temp2[0] = -sin(phi);
		temp2[1] = cos(phi);
		temp3 = vector_scalar_multiplication(
				sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2)),
				temp2);
		vector<float> CSC_LR_S_finalpoint = vector_adding(CSC_LR_S_initialpoint,
				temp3); // point of contact of the final circle with the common tangent
		//vector_print_float(2, CSC_LR_S_initialpoint);
		//vector_print_float(2, CSC_LR_S_finalpoint);

		// Find the initial arc
		vector<float> R_pinit = vector_subtraction(pinit,
				initialCircleLeft_centre);
		float phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
		vector<float> R_contact = vector_subtraction(CSC_LR_S_initialpoint,
				initialCircleLeft_centre);
		float phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		float initialArcAngle;
		if (phi_pinit >= phi_contact)
			initialArcAngle = 2 * PI - (phi_pinit - phi_contact);
		else
			initialArcAngle = phi_contact - phi_pinit;

		CSC_LR_length = CSC_LR_length + R * initialArcAngle;
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;

		// Find the final arc

		R_contact = vector_subtraction(CSC_LR_S_finalpoint,
				finalCircleRight_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
		vector<float> R_pf = vector_subtraction(pf, finalCircleRight_centre);
		float phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

		float finalArcAngle;
		if (phi_contact >= phi_pf)
			finalArcAngle = phi_contact - phi_pf;
		else
			finalArcAngle = 2 * PI - (phi_pf - phi_contact);

		CSC_LR_length = CSC_LR_length + R * finalArcAngle;
		//cout << "CSC_LR_length=" << CSC_LR_length << endl;

		float LSR_length = CSC_LR_length;

		result_vector.push_back(LSR_length);
		result_vector.push_back(initialCircleLeft_centre[0]);
		result_vector.push_back(initialCircleLeft_centre[1]);
		result_vector.push_back(initialArcAngle);
		result_vector.push_back(phi_pinit);
		result_vector.push_back(v_norm);
		result_vector.push_back(finalCircleRight_centre[0]);
		result_vector.push_back(finalCircleRight_centre[1]);
		result_vector.push_back(finalArcAngle);
		result_vector.push_back(phi_contact);

		return result_vector;
	} else {
		CSC_LR_length = -999;
		float LSR_length = CSC_LR_length;
		result_vector.push_back(LSR_length);
		return result_vector;
	}

	//cout << "LSR_length = " << LSR_length << endl;

	// --------------------------------------------------------------------
	// END: LSR
	// --------------------------------------------------------------------

}

vector<float> dubins_RSL(vector<float> pinit, float thinit, vector<float> pf,
		float thf, float R) {
	vector<float> result_vector;

	vector<float> temp1;
	temp1.push_back(-999);
	temp1.push_back(-999);
	vector<float> temp2;
	temp2.push_back(-999);
	temp2.push_back(-999);
	vector<float> temp3;
	temp3.push_back(-999);
	temp3.push_back(-999);

	vector<float> vi;
	vi.push_back(cos(thinit));
	vi.push_back(sin(thinit));
	//cout << "vi = [" << vi[0] << "," << vi[1] << "]" << endl;

	vector<float> vf;
	vf.push_back(cos(thf));
	vf.push_back(sin(thf));
	//cout << "vf = [" << vf[0] << "," << vf[1] << "]" << endl;

	vector<float> initialCircleRight_centre;
	initialCircleRight_centre.push_back(pinit[0] + R * vi[1]);
	initialCircleRight_centre.push_back(pinit[1] - R * vi[0]);
	//cout << "initialCircleRight_centre = [" << initialCircleRight_centre[0] << "," << initialCircleRight_centre[1] << "]" << endl;

	vector<float> finalCircleLeft_centre;
	finalCircleLeft_centre.push_back(pf[0] - R * vf[1]);
	finalCircleLeft_centre.push_back(pf[1] + R * vf[0]);
	//cout << "finalCircleLeft_centre = [" << finalCircleLeft_centre[0] << "," << finalCircleLeft_centre[1] << "]" << endl;

	float initialCircleRight_radius = R;
	float finalCircleLeft_radius = R;
	//cout << "initialCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "initialCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleLeft_radius = [" << initialCircleLeft_radius << "]" << endl;
	//cout << "finalCircleRight_radius = [" << initialCircleLeft_radius << "]" << endl;

	vector<float> initialCircleRight_orientation;
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(0);
	initialCircleRight_orientation.push_back(-1);
	//cout << "initialCircleRight_orientation = ["	<< initialCircleRight_orientation[0] << ","
	//												<< initialCircleRight_orientation[1] << ","
	//												<< initialCircleRight_orientation[2] << "]" << endl;

	vector<float> finalCircleLeft_orientation;
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(0);
	finalCircleLeft_orientation.push_back(1);
	//cout << "finalCircleLeft_orientation = ["	<< finalCircleLeft_orientation[0] << ","
	//											<< finalCircleLeft_orientation[1] << ","
	//											<< finalCircleLeft_orientation[2] << "]" << endl;

	// --------------------------------------------------------------------
	// START: RSL
	// --------------------------------------------------------------------

	float CSC_RL_length = 0;

	vector<float> v = vector_subtraction(finalCircleLeft_centre,
			initialCircleRight_centre);
	float v_angle = atan2_0_2pi(v[1], v[0]);
	if (vector_norm(v) >= 2 * R) {
		float v_norm = sqrt(
				vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2));
		CSC_RL_length = CSC_RL_length + v_norm;
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;

		float phi = v_angle
				+ atan2_0_2pi(
						sqrt(
								vector_sum(vector_multiplication(v, v))
										- 4 * pow(R, 2)), 2 * R);

		temp1[0] = R * cos(phi);
		temp1[1] = R * sin(phi);
		vector<float> CSC_RL_S_initialpoint = vector_adding(
				initialCircleRight_centre, temp1); // point of contact of the initial circle with the common tangent
		//vector_print_float(2, CSC_RL_S_initialpoint);
		temp2[0] = sin(phi);
		temp2[1] = -cos(phi);
		temp3 = vector_scalar_multiplication(
				sqrt(vector_sum(vector_multiplication(v, v)) - 4 * pow(R, 2)),
				temp2);
		vector<float> CSC_RL_S_finalpoint = vector_adding(CSC_RL_S_initialpoint,
				temp3); // point of contact of the final circle with the common tangent
		//vector_print_float(2, CSC_RL_S_finalpoint);

		// Find the initial arc

		vector<float> R_pinit = vector_subtraction(pinit,
				initialCircleRight_centre);
		float phi_pinit = atan2_0_2pi(R_pinit[1], R_pinit[0]);
		vector<float> R_contact = vector_subtraction(CSC_RL_S_initialpoint,
				initialCircleRight_centre);
		float phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);

		float initialArcAngle;
		if (phi_pinit >= phi_contact)
			initialArcAngle = phi_pinit - phi_contact;
		else
			initialArcAngle = 2 * PI - (phi_contact - phi_pinit);

		CSC_RL_length = CSC_RL_length + R * initialArcAngle;
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;

		// Find the final arc

		R_contact = vector_subtraction(CSC_RL_S_finalpoint,
				finalCircleLeft_centre);
		phi_contact = atan2_0_2pi(R_contact[1], R_contact[0]);
		vector<float> R_pf = vector_subtraction(pf, finalCircleLeft_centre);
		float phi_pf = atan2_0_2pi(R_pf[1], R_pf[0]);

		float finalArcAngle;
		if (phi_contact >= phi_pf)
			finalArcAngle = 2 * PI - (phi_contact - phi_pf);
		else
			finalArcAngle = phi_pf - phi_contact;

		CSC_RL_length = CSC_RL_length + R * finalArcAngle;
		//cout << "CSC_RL_length=" << CSC_RL_length << endl;
		float RSL_length = CSC_RL_length;

		result_vector.push_back(RSL_length);
		result_vector.push_back(initialCircleRight_centre[0]);
		result_vector.push_back(initialCircleRight_centre[1]);
		result_vector.push_back(initialArcAngle);
		result_vector.push_back(phi_pinit);
		result_vector.push_back(v_norm);
		result_vector.push_back(finalCircleLeft_centre[0]);
		result_vector.push_back(finalCircleLeft_centre[1]);
		result_vector.push_back(finalArcAngle);
		result_vector.push_back(phi_contact);

		return result_vector;

	} else {
		CSC_RL_length = -999;
		float RSL_length = CSC_RL_length;
		result_vector.push_back(RSL_length);
		return result_vector;
	}

	//cout << "RSL_length = " << RSL_length << endl;

	// --------------------------------------------------------------------
	// END: RSL
	// --------------------------------------------------------------------
}

float eq_LSL(float alpha, float beta, float d, float r, float k) {
	float p = r
			* pow(
					2 + pow(d / r, 2) - 2 * cos(alpha - beta)
							+ 2 * (d / r) * (sin(alpha) - sin(beta)), 0.5);

	float b_m_a = beta - alpha;

	while (b_m_a >= 2 * PI) {
		b_m_a = b_m_a - 2 * PI;
	}

	while (b_m_a < 0) {
		b_m_a = b_m_a + 2 * PI;
	}

	float L = (b_m_a) * r + p + k * r;

	return L;
}

float eq_RSR(float alpha, float beta, float d, float r, float k) {
	float p = r
			* pow(
					2 + pow(d / r, 2) - 2 * cos(alpha - beta)
							+ 2 * (d / r) * (sin(beta) - sin(alpha)), 0.5);

	float a_m_b = alpha - beta;

	while (a_m_b >= 2 * PI) {
		a_m_b = a_m_b - 2 * PI;
	}
	while (a_m_b < 0) {
		a_m_b = a_m_b + 2 * PI;
	}

	float L = (alpha - beta) * r + p + k * r;

	return L;
}

float eq_LSR(float alpha, float beta, float d, float r, float k) {
	float L;

	float terminal_point_eq = pow(d / r, 2) - 2 + 2 * cos(alpha - beta)
			+ 2 * (d / r) * (sin(alpha) + sin(beta));

	if (terminal_point_eq < 0) {
		L = -999;
		return L;
	}

	float p = r * pow(terminal_point_eq, 0.5);

	float t = r
			* (-alpha
					+ atan2(-cos(alpha) - cos(beta),
							d / r + sin(alpha) + sin(beta)) - atan2(-2, p / r));

	L = (alpha - beta) * r + 2 * t + p + k * r;

	return L;
}

float eq_RSL(float alpha, float beta, float d, float r, float k) {
	float L;
	float terminal_point_eq = pow(d / r, 2) - 2 + 2 * cos(alpha - beta)
			- 2 * (d / r) * (sin(alpha) + sin(beta));

	if (terminal_point_eq < 0) {
		L = -999;
		return L;
	}

	float p = r * pow(terminal_point_eq, 0.5);

	float t = r
			* (alpha
					- atan2(cos(alpha) + cos(beta),
							d / r - sin(alpha) - sin(beta)) + atan2(2, p / r));

	L = (-alpha + beta) * r + 2 * t + p + k * r;

	return L;
}

float find_terminal_point(float alpha, float beta, float R_min, float R_max,
		float d, char method_str) {
	// This function finds the terminal LSR and RSL radiuses of opposite direction circles of Dubins curves
	// in order to find maximum relevant R value for the requierd L.
	// The function calculates the second order polynom for the solution and
	// solves it using the 'second_order_pol_solution' function.

	float A = -2 + 2 * cos(alpha - beta);
	float B;
	switch (method_str) {
	case 'L': // 'LSR'
		B = 2 * d * (sin(alpha) + sin(beta));
		break;
	case 'R': // 'RSL'
		B = -2 * d * (sin(alpha) + sin(beta));
		break;
	}
	float C = pow(d, 2);

	float R_target = second_order_pol_solution(A, B, C, R_min, R_max);

	return R_target;
}

float interp_similar_circles(float alpha, float beta, float R_min, float R_max,
		float d, float L_target, char method_str, float k) {

	// This function solves the function of similar direction circles of Dubins curves
	// in order to find R value for the requierd L.
	// The function calculates the second order polynom for the solution and
	// solves it using the 'second_order_pol_solution' function.

	float L_min, L_max;

	switch (method_str) {
	case 'L': // 'LSL'
		L_min = eq_LSL(alpha, beta, d, R_min, k);
		L_max = eq_LSL(alpha, beta, d, R_max, k);
		break;
	case 'R': // 'RSR'
		L_min = eq_RSR(alpha, beta, d, R_min, k);
		L_max = eq_RSR(alpha, beta, d, R_max, k);
		break;
	}

	if ((L_target < L_min) || (L_target > L_max) || (L_max < L_min))
		return -999;

	float A, B, C;
	switch (method_str) {
	case 'L': // 'LSL'
		A = pow(beta - alpha + k, 2) - 2 + 2 * cos(alpha - beta);
		B = -2 * L_target * (beta - alpha + k)
				- 2 * d * (sin(alpha) - sin(beta));
		break;
	case 'R': // 'RSR'
		A = pow(alpha - beta + k, 2) - 2 + 2 * cos(alpha - beta);
		B = -2 * L_target * (alpha - beta + k)
				- 2 * d * (sin(beta) - sin(alpha));
		break;
	}
	C = pow(L_target, 2) - pow(d, 2);

	float R_target = second_order_pol_solution(A, B, C, R_min, R_max);

	return R_target;
}

float interp_opp_circles(float alpha, float beta, float R_min, float R_max,
		float d, float L_target, float dL, char method_str, float k) {

	// This function interpulates the function of opposite circles Dubins curves
	// in order to find R value for the requierd L.
	// The function begins by calculates 3 point and inserts them to R_vec and
	// L_vec, those are : R = 0, min_R and max_R and their corespondig L values: d, L_min and L_max.
	// the interpulation is calculationg a parabolic equation astimation and redefines the R values
	// until reaching L valus which is far from L_target by no more than dL.

	int max_inter = 10;

	float L_min, L_max;

//	cout << dL << endl;

	switch (method_str) {
	case 'L': // 'LSR'
		L_min = eq_LSR(alpha, beta, d, R_min, k);
		L_max = eq_LSR(alpha, beta, d, R_max, k);
		if (L_max == -999) {
			R_max = find_terminal_point(alpha, beta, R_min, R_max, d, 'L')
					- 0.0001;
			L_max = eq_LSR(alpha, beta, d, R_max, k);
		}
		break;
	case 'R': // 'RSL'
		L_min = eq_RSL(alpha, beta, d, R_min, k);
		L_max = eq_RSL(alpha, beta, d, R_max, k);
		if (L_max == -999) {
			R_max = find_terminal_point(alpha, beta, R_min, R_max, d, 'R')
					- 0.0001;
			L_max = eq_RSL(alpha, beta, d, R_max, k);
		}
		break;
	}

	int n_iter = 1;

	if ((L_target < L_min) || (L_target > L_max) || (L_max < L_min))
		return -999;

	vector<float> abc_vec = second_order_equ_calc(0, d, R_min, L_min, R_max,
			L_max);
	float R_new = second_order_pol_solution(abc_vec[0], abc_vec[1],
			abc_vec[2] - L_target, R_min, R_max);

	if ((R_new == -999) || (R_new < R_min))
		return -999;

	float L_new;
	switch (method_str) {
	case 'L': // 'LSR'
		L_new = eq_LSR(alpha, beta, d, R_new, k);
		break;
	case 'R': // 'RSL'
		L_new = eq_RSL(alpha, beta, d, R_new, k);
		break;
	}

	float R_left = 0;
	float R_right = R_max;
	float R_mid = R_min;

	float L_left = d;
	float L_right = L_max;
	float L_mid = L_min;

	while ((fabs(L_new - L_target) > dL) && (n_iter < max_inter)) {
		n_iter = n_iter + 1;

		if (R_new > R_mid) {
			R_left = R_mid;
			R_mid = R_new;
			L_left = L_mid;
			L_mid = L_new;
		} else {
			if (R_new < R_mid) {
				R_right = R_mid;
				R_mid = R_new;
				L_right = L_mid;
				L_mid = L_new;
			}
		}
		vector<float> abc_vec = second_order_equ_calc(R_left, L_left, R_mid,
				L_mid, R_right, L_right);
		R_new = second_order_pol_solution(abc_vec[0], abc_vec[1],
				abc_vec[2] - L_target, R_min, R_max);

		switch (method_str) {
		case 'L': // 'LSR'
			L_new = eq_LSR(alpha, beta, d, R_new, k);
			break;
		case 'R': // 'RSL'
			L_new = eq_RSL(alpha, beta, d, R_new, k);
			break;
		}

	}
	float R_target = R_new;
	if (n_iter >= max_inter)
		cout << "n_iter=" << n_iter << endl;

	if ((R_target < R_min) || (R_target > R_max) || (n_iter >= max_inter))
		return -999;

	return R_target;
}

vector<float> Dubins_state_time_steer(vector<float> vector_start,
		vector<float> vector_end) {

	float dL_obs_S = sqrt(pow(ds_obs, 2) * (1 + 1 / pow(V_max, 2)));
//					cout << dL_obs_S << endl;
//					vector_print_float(states_number, vector_start);
//					vector_print_float(states_number, vector_end);

	vector<float> result;
	vector<float> path_result_new;
	vector<float> path_result_temp;

	vector<float> pinit, pf;
	float thinit, thf;
	float R;

	pinit.push_back(vector_start[0]);
	pinit.push_back(vector_start[1]);
	thinit = vector_start[2];
	pf.push_back(vector_end[0]);
	pf.push_back(vector_end[1]);
	thf = vector_end[2];

	float d = vector_norm(vector_subtraction(pf, pinit));
	float LOS_ang = atan2(pf[1] - pinit[1], pf[0] - pinit[0]);
	float alpha = thinit - LOS_ang;
	float beta = thf - LOS_ang;

	alpha = angle_to_pi_interval(alpha);
	beta = angle_to_pi_interval(beta);
	//	cout << "d=" << d << endl;
	//	cout << "alpha=" << alpha << endl;
	//	cout << "beta=" << beta << endl;
	float L_target = (vector_end[3] - vector_start[3]) * V_max;

	float R_new, R_temp;
	float c_new = 1; // C_newL: 1-LSL,0*PI, 2-RSR,0*PI, 3-LSR,0*PI, 4-RSL,0*PI, 11-LSL,2*PI, 12-RSR,2*PI, 13-LSR,2*PI, 14-RSL,2*PI
	path_result_new.push_back(-999);

	// attemps to solve using similar circles path
	R_new = interp_similar_circles(alpha, beta, Sim_R_min, Sim_R_max, d,
			L_target, 'L', 0); // R='RSR' ,  L='LSL'
	if (R_new != -999) {
		path_result_temp = dubins_LSL(pinit, thinit, pf, thf, R_new);
//						vector_print_float(10, path_result_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_new = -999;
		else
			path_result_new = vector_copy_float(path_result_temp.size(),
					path_result_temp);
	}

	R_temp = interp_similar_circles(alpha, beta, Sim_R_min, Sim_R_max, d,
			L_target, 'R', 0); // R='RSR' ,  L='LSL'
	if (R_temp != -999) {
		path_result_temp = dubins_RSR(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}

	result.push_back(R_new);
	result.push_back(c_new);
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 2;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}

	// attemps to solve using opposite circles path
	R_temp = interp_opp_circles(alpha, beta, Sim_R_min, Sim_R_max, d, L_target,
			dL_obs_S / 10, 'L', 0); // R='RSL' ,  L='LSR'
	if (R_temp != -999) {
		path_result_temp = dubins_LSR(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 3;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	R_temp = interp_opp_circles(alpha, beta, Sim_R_min, Sim_R_max, d, L_target,
			dL_obs_S / 10, 'R', 0); // R='RSL' ,  L='LSR'
	if (R_temp != -999) {
		path_result_temp = dubins_RSL(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 4;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	if (R_new != -999) {
		result[0] = R_new;
		result[1] = c_new;
		result.push_back(path_result_new[0]);
		result.push_back(path_result_new[1]);
		result.push_back(path_result_new[2]);
		result.push_back(path_result_new[3]);
		result.push_back(path_result_new[4]);
		result.push_back(path_result_new[5]);
		result.push_back(path_result_new[6]);
		result.push_back(path_result_new[7]);
		result.push_back(path_result_new[8]);
		result.push_back(path_result_new[9]);
		return result;
	}

	// attemps to solve using similar circles path and a 2_pi path
	R_temp = interp_similar_circles(alpha, beta, Sim_R_min, Sim_R_max, d,
			L_target, 'L', 2 * PI); // R='RSR' ,  L='LSL'

	if (R_temp != -999) {
		path_result_temp = dubins_LSL(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 11;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	R_temp = interp_similar_circles(alpha, beta, Sim_R_min, Sim_R_max, d,
			L_target, 'R', 2 * PI); // R='RSR' ,  L='LSL'
	if (R_temp != -999) {
		path_result_temp = dubins_RSR(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 12;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	if (R_new != -999) {
		result[0] = R_new;
		result[1] = c_new;
		result.push_back(path_result_new[0]);
		result.push_back(path_result_new[1]);
		result.push_back(path_result_new[2]);
		result.push_back(path_result_new[3]);
		result.push_back(path_result_new[4]);
		result.push_back(path_result_new[5]);
		result.push_back(path_result_new[6]);
		result.push_back(path_result_new[7]);
		result.push_back(path_result_new[8]);
		result.push_back(path_result_new[9]);
		return result;
	}

	// attemps to solve using opposite circles path and a 2_pi path
	R_temp = interp_opp_circles(alpha, beta, Sim_R_min, Sim_R_max, d, L_target,
			dL_obs_S / 2, 'L', 2 * PI); // R='RSL' ,  L='LSR'
	if (R_temp != -999) {
		path_result_temp = dubins_LSR(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 13;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	R_temp = interp_opp_circles(alpha, beta, Sim_R_min, Sim_R_max, d, L_target,
			dL_obs_S / 2, 'R', 2 * PI); // R='RSL' ,  L='LSR'
	if (R_temp != -999) {
		path_result_temp = dubins_RSL(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 14;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	if (R_new != -999) {
		result[0] = R_new;
		result[1] = c_new;
		result.push_back(path_result_new[0]);
		result.push_back(path_result_new[1]);
		result.push_back(path_result_new[2]);
		result.push_back(path_result_new[3]);
		result.push_back(path_result_new[4]);
		result.push_back(path_result_new[5]);
		result.push_back(path_result_new[6]);
		result.push_back(path_result_new[7]);
		result.push_back(path_result_new[8]);
		result.push_back(path_result_new[9]);
		return result;
	}

	// attemps to solve using similar circles path and a 4_pi path
	R_temp = interp_similar_circles(alpha, beta, Sim_R_min, Sim_R_max, d,
			L_target, 'L', 4 * PI); // R='RSR' ,  L='LSL'

	if (R_temp != -999) {
		path_result_temp = dubins_LSL(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 21;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	R_temp = interp_similar_circles(alpha, beta, Sim_R_min, Sim_R_max, d,
			L_target, 'R', 4 * PI); // R='RSR' ,  L='LSL'
	if (R_temp != -999) {
		path_result_temp = dubins_RSR(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 22;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	if (R_new != -999) {
		result[0] = R_new;
		result[1] = c_new;
		result.push_back(path_result_new[0]);
		result.push_back(path_result_new[1]);
		result.push_back(path_result_new[2]);
		result.push_back(path_result_new[3]);
		result.push_back(path_result_new[4]);
		result.push_back(path_result_new[5]);
		result.push_back(path_result_new[6]);
		result.push_back(path_result_new[7]);
		result.push_back(path_result_new[8]);
		result.push_back(path_result_new[9]);
		return result;
	}

	// attemps to solve using opposite circles path and a 2_pi path
	R_temp = interp_opp_circles(alpha, beta, Sim_R_min, Sim_R_max, d, L_target,
			dL_obs_S / 2, 'L', 4 * PI); // R='RSL' ,  L='LSR'
	if (R_temp != -999) {
		path_result_temp = dubins_LSR(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 23;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	R_temp = interp_opp_circles(alpha, beta, Sim_R_min, Sim_R_max, d, L_target,
			dL_obs_S / 2, 'R', 4 * PI); // R='RSL' ,  L='LSR'
	if (R_temp != -999) {
		path_result_temp = dubins_RSL(pinit, thinit, pf, thf, R_temp);
		//cout << path_result_temp[0] << endl;
		if (fabs(path_result_temp[0] - L_target) > dL_obs_S)
			R_temp = -999;
	}
	if (((R_new == -999) && (R_temp != -999))
			|| ((R_temp < R_new) && (R_temp != -999))) {
		R_new = R_temp;
		c_new = 24;
		path_result_new = vector_copy_float(path_result_temp.size(),
				path_result_temp);
	}
	if (R_new != -999) {
		result[0] = R_new;
		result[1] = c_new;
		result.push_back(path_result_new[0]);
		result.push_back(path_result_new[1]);
		result.push_back(path_result_new[2]);
		result.push_back(path_result_new[3]);
		result.push_back(path_result_new[4]);
		result.push_back(path_result_new[5]);
		result.push_back(path_result_new[6]);
		result.push_back(path_result_new[7]);
		result.push_back(path_result_new[8]);
		result.push_back(path_result_new[9]);
		return result;
	}

	return result;

	//	/////////////   dubins algorithm check	///////////////
	//	float best_solution = dubins(pinit, thinit, pf, thf, R);
	//	cout << "full dubins best_solution=" << best_solution << endl;
	//
	//	best_solution = dubins_four(pinit, thinit, pf, thf, R);
	//	cout << "four dubins best_solution=" << best_solution << endl;
	//
	//	best_solution = dubins_LSL(pinit, thinit, pf, thf, R);
	//	cout << "LSL dubins best_solution=" << best_solution << endl;
	//
	//	best_solution = dubins_RSR(pinit, thinit, pf, thf, R);
	//	cout << "RSR dubins best_solution=" << best_solution << endl;
	//
	//	best_solution = dubins_LSR(pinit, thinit, pf, thf, R);
	//	cout << "LSR dubins best_solution=" << best_solution << endl;
	//
	//	best_solution = dubins_RSL(pinit, thinit, pf, thf, R);
	//	cout << "RSL dubins best_solution=" << best_solution << endl;

	///////////////   second_order_equ_calc check	///////////////
	//float x1 = 2;
	//float y1 = 5;
	//float x2 = 1;
	//float y2 = 7;
	//float x3 = -3;
	//float y3 = 6;
	//vector<float> abc_vec = second_order_equ_calc(x1, y1, x2, y2, x3, y3);
	//cout << "A=" << abc_vec[0] << " B=" << abc_vec[1] << " C=" << abc_vec[2] << endl;
	//cout << "equation is: " << abc_vec[0] << "x^2+" << abc_vec[1] << "x+" << abc_vec[2] << endl;

	/////////////////   interp_opp_circles check	///////////////
	//float L_target = 1.75;
	//float R_target = interp_opp_circles(alpha, beta, 0.2, 0.55, d, L_target, 0.0001, 'R', 0); // R='RSL' ,  L='LSR'
	//cout << "R_target=" << R_target << endl;
	//float L_target_new = eq_RSL(alpha, beta, d, R_target, 0);
	//cout << "L_target_new=" << L_target_new << endl;

	/////////////////   interp_similar_circles check	///////////////
	//float L_target = 10;
	//float R_target = interp_similar_circles(alpha, beta, 0.5, 1.6, d, L_target, 'L', 2*PI); // R='RSR' ,  L='LSL'
	//cout << "R_target=" << R_target << endl;
	//float L_target_new = eq_LSL(alpha, beta, d, R_target, 2 * PI);
	//cout << "L_target_new=" << L_target_new << endl;

}

//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// main /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	int age;
//
//	int row_num=10000;
//	int column_num=10000;
//
//	for (int i=0;i<10000;i++)
//	{
//		cout << i << endl;
//		float** consumer_matrix;
//		consumer_matrix = matrix_creat_float(row_num, column_num);
//		matrix_delete_float(consumer_matrix,row_num);
//		consumer_matrix = matrix_creat_float(row_num, column_num);
//
//	}

	int vertices_col = steps + 1;
	edges_col = steps;
	traj_col = steps;

	vertices_final = matrix_creat_float(states_number, vertices_col);
	edges_final = matrix_creat_int(edges_row, edges_col);
	traj_final = matrix_creat_float(traj_row, traj_col);

	if (solution_mode == 1) {
		main_RRT_sim(vertices_final, edges_final, traj_final);

		cout
				<< " insert plot_mode key for tree graph (1) - only tree  (2)- only best solution  (3) - full tree and solution:";
		cin >> age;

		plot_mode = age;

		if (plot_mode == 2) {
			if ((hit_index_RD == -999) && (hit_index == -999)) {
				cout << " ERROR!!! - No solution was found!" << endl;
			}
			if ((hit_index_RD != -999) && (hit_index == -999)) {
				Dubins_sol_mode = 1;
				cout << " Solution presented for relaxed Dubins case" << endl;
			}
			if ((hit_index != -999) && (hit_index_RD == -999)) {
				Dubins_sol_mode = 2;
				cout << " Solution presented for full Dubins case" << endl;
			}
			if ((hit_index_RD != -999) && (hit_index != -999)) {
				cout
						<< " Select solution type (1) - Relaxed Dubins  (2)- Full Dubins: ";
				int age;
				cin >> age;
				Dubins_sol_mode = age;
			}

		}

		comm_mov.push_back(0);
		comm_mov.push_back(0);
		comm_mov.push_back(0);

		// init GLUT and create window
		glutInit(&argc, argv);
		glutInitWindowPosition(window_pos_left, window_pos_hight);
		glutInitWindowSize(window_width, window_height);
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
		glutCreateWindow("main1");

		// register callbacks
		glutDisplayFunc(renderScene);
		glutReshapeFunc(changeSize);
		glutIdleFunc(renderScene);

		glClearColor(1.0, 1.0, 1.0, 1.0);

		glutIgnoreKeyRepeat(1);
		glutMouseFunc(mouseButton);
		glutMotionFunc(mouseMove);

		glutKeyboardFunc(processNormalKeys);

		// enter GLUT event processing cycle
		glutMainLoop();
	}

	if (solution_mode == 2) {
		cout
				<< " Choose a solution method (1) - new RRT simulation  (2) - read results from 'results_file_1000000_50_1.txt' file:";
		int age;
		cin >> age;
		if (age == 1)
			main_RRT_sim(vertices_final, edges_final, traj_final);

		sizes_values_final = read_solution_file_sizes();

//								cout << sizes_values_final[0] << " " << sizes_values_final[1] << " " << sizes_values_final[2] << " " << sizes_values_final[3] << endl;

		solutions_states_relaxed_Dubins_final = matrix_creat_float(
				monte_carlo_number * states_number, sizes_values_final[0]);
		solutions_marks_relaxed_Dubins_final = matrix_creat_float(
				monte_carlo_number, sizes_values_final[1]);
		solutions_states_full_Dubins_final = matrix_creat_float(
				monte_carlo_number * states_number, sizes_values_final[2]);
		solutions_marks_full_Dubins_final = matrix_creat_float(
				monte_carlo_number, sizes_values_final[3]);

		solutions_TTH_Relaxed_Dubins_matrix_final = matrix_creat_float(
				monte_carlo_number, convergence_size);
		solutions_TTH_full_Dubins_matrix_final = matrix_creat_float(
				monte_carlo_number, convergence_size);
		float** vectors_matrix_temp = matrix_creat_float(3, convergence_size);

		read_solution_file(solutions_states_relaxed_Dubins_final,
				solutions_marks_relaxed_Dubins_final,
				solutions_states_full_Dubins_final,
				solutions_marks_full_Dubins_final,
				solutions_TTH_Relaxed_Dubins_matrix_final,
				solutions_TTH_full_Dubins_matrix_final, vectors_matrix_temp,
				sizes_values_final);

		//			      	  	  	  	  matrix_print_float(12, 46,solutions_states_relaxed_Dubins);
		//			      	  	  	  	  matrix_print_float(3, 540,solutions_marks_relaxed_Dubins);
		//			      	  	  	  	  matrix_print_float(12, 107,solutions_states_full_Dubins);
		//			      	  	  	  	  matrix_print_float(3, 1272,solutions_marks_full_Dubins);
		//			      	  	  	  	  matrix_print_float(3, 8,solutions_TTH_Relaxed_Dubins_matrix);
		//			      	  	  	  	  matrix_print_float(3, 8,solutions_TTH_full_Dubins_matrix);
		//			      	  	  	  	  matrix_print_float(3, 8,vectors_matrix_temp);

		for (int i = 0; i < convergence_size; i++) {
			monte_carlo_avg_Relaxed_Dubins_final.push_back(
					vectors_matrix_temp[0][i]);
			monte_carlo_avg_full_Dubins_final.push_back(
					vectors_matrix_temp[1][i]);
			convergence_values_final.push_back(vectors_matrix_temp[2][i]);
		}

		//	cout << "hit index in main:" << endl;
		//	cout << hit_index << endl;

		if (solution_mode == 2) {
			cout
					<< " insert plot_mode key for tree graph (1) - all relaxed dubins solutions (2)- only best relaxed dubins solution,"
					<< endl;
			cout
					<< " 									  (3) - all full dubins solutions (4)- only best full dubins solution,"
					<< endl;
			cout
					<< " 									  (5) - best euclidean, relaxed and full solutions with optimal solutions:"
					<< endl;
			cin >> age;

			plot_mode = age;

			color_convergence_matrix = matrix_creat_float(monte_carlo_number,
					3);
			for (int i = 0; i < monte_carlo_number; i++) {
				for (int j = 0; j < 3; j++) {
					color_convergence_matrix[i][j] = uniform(0.0, 1.0);
				}
			}

			comm_mov.push_back(0);
			comm_mov.push_back(0);
			comm_mov.push_back(0);

			// init GLUT and create window
			glutInit(&argc, argv);
			glutInitWindowPosition(window_pos_left, window_pos_hight);
			glutInitWindowSize(window_width, window_height);
			glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
			glutCreateWindow("main1");

			// register callbacks
			glutDisplayFunc(renderScene);
			glutReshapeFunc(changeSize);
			glutIdleFunc(renderScene);

			glClearColor(1.0, 1.0, 1.0, 1.0);

			glutIgnoreKeyRepeat(1);
			glutMouseFunc(mouseButton);
			glutMotionFunc(mouseMove);

			glutKeyboardFunc(processNormalKeys);

			// enter GLUT event processing cycle
			glutMainLoop();
		}
	}

	cout << "sim finished" << endl;
	cin >> age;
	return 1;

}


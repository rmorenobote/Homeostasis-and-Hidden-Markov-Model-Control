

#define _CRT_SECURE_NO_WARNINGS //used only for f_open, but check if other warnings arise
#include <iostream>
#include <random>

/* function declaration */
double* beta(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2,
	int door_x1, double p_y_food, double p_y_no_food, double prob_switch);

double trajectory_beta(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1, 
	double p_y_food, double p_y_no_food, double* beta_vec, int x1_0, int x2_0, double prob_switch);

double* value_function(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2,
	int door_x1, double p_y_food, double p_y_no_food, double alpha, double prob_switch);

double trajectory_value(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1,
	double p_y_food, double p_y_no_food, double* value_vec, int x1_0, int x2_0, double alpha, double prob_switch);

int* trajectory_beta_samples(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1,
	double p_y_food, double p_y_no_food, double* beta_vec, int x1_0, int x2_0, double prob_switch);

int* trajectory_value_samples(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1,
	double p_y_food, double p_y_no_food, double* value_vec, int x1_0, int x2_0, double alpha, double prob_switch);



int main(void)
{
	//Parameters

	const int     repetitions = 10000; //repetitions to compute averages

	const int	  T = 20; //horizon

	const int     grid_size_x1 = 10; //grid size
	const int     grid_size_x2 = 10;

	//make sure that the food source is is a location where all 5 movements are possible; codes below are only valid in this case
	const int     x1_food = 7; //food location, with range (0, grid_size_x1 - 1)
	const int     x2_food = 5;

	const double  prob_switch = 0.2; //prob of swithing from 0 to the other 4 cells

	const double  prob_switch_ini = 0;
	const double  prob_switch_end = 0.5;
	const int     num_prob_switch = 20;

	//we add walls at middle_wal x2, with a door at (door_x1, middle_wall_x2); actions are fobidden towards the walls, 
	// including from the door space
	const int     door_x1 = 1;
	const int     middle_wall_x2 = 3;
	
	const int     x1_0 = 3; //initial condition trajectory
	const int     x2_0 = 0;

	const double  alpha = 200;

	const double  alpha_ini = 1;
	const double  alpha_end = 300;
	const int     num_alpha = 20;

	const double  p_y_food = 1.00; //y is binary 1, and 0, p_y_food is the prob of y=1 at the food, independent of actions in our case
	const double  p_y_no_food = 0.9; //prob of y=1 when not at the food

	const int     seed_r = 369;
	
	//Definitions
	int	     rep, t, i1, i2, i;

	double   beta_aux, beta_ini, value_aux;
	double   beta_acum_HMM, beta_acum_variational, beta_aver_HMM, beta_aver_variational;

	double   beta2_aver_HMM, beta2_aver_variational, error_beta_aver_HMM, error_beta_aver_variational;

	double   alpha0, d_alpha, prob_switch0, d_prob_switch;

	double   *beta_vec, *value_vec;
	int* trajectory_beta_vec, * trajectory_value_vec;

	beta_vec = (double*)malloc(grid_size_x1 * grid_size_x2 * T * sizeof(double));
	value_vec = (double*)malloc(grid_size_x1 * grid_size_x2 * T * sizeof(double));
	trajectory_beta_vec = (int*)malloc(2 * T * sizeof(int));
	trajectory_value_vec = (int*)malloc(2 * T * sizeof(int));

	FILE* trajectories_HMM, * trajectories_var_inf;
	FILE *beta_values_alpha, * beta_values_prob_switch, *bellman_values;

	beta_values_alpha = fopen("beta_values_alpha.m", "w");
	beta_values_prob_switch = fopen("beta_values_prob_switch.m", "w");
	bellman_values = fopen("bellman_values.m", "w");

	trajectories_HMM = fopen("trajectories_HMM.m", "w");
	trajectories_var_inf = fopen("trajectories_var_inf.m", "w");

	//random generator
	std::mt19937 gen(seed_r); // Standard mersenne_twister_engine seeded with constant seed_r
	//std::normal_distribution<> d{ 0., sigma }; //normal distribution with mean mu and std sigma
	//std::uniform_int_distribution<> dis_unif(1, 6); //discrete uniform, range of values used in experiment
	std::uniform_real_distribution<> uniform(0.0, 1.0);
	
	//ploting here just one trajectory for each model, for just 100 trajectories
	//betas, computing prob p(y_{t:T}=y_d);
	beta_vec = beta(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1, p_y_food, p_y_no_food, prob_switch);

	value_vec = value_function(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1, p_y_food, p_y_no_food, 
		alpha, prob_switch);

	for (rep = 0; rep < 100; rep++) {
		trajectory_beta_vec = trajectory_beta_samples(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1,
			p_y_food, p_y_no_food, beta_vec, x1_0, x2_0, prob_switch);

		trajectory_value_vec = trajectory_value_samples(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1,
			p_y_food, p_y_no_food, value_vec, x1_0, x2_0, alpha, prob_switch);

		for (t = 0; t < T; t++) {
			fprintf(trajectories_HMM, "%d %d %d \n", t, trajectory_beta_vec[2 * t], trajectory_beta_vec[2 * t + 1]);
			fprintf(trajectories_var_inf, "%d %d %d \n", t, trajectory_value_vec[2 * t], trajectory_value_vec[2 * t + 1]);
		}
	}

	for (t = T - 1; t >= 0; t--) { //backwards from t=T to t=1
		for (i1 = 0; i1 < grid_size_x1; i1++) {
			for (i2 = 0; i2 < grid_size_x2; i2++) {
				beta_aux = beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t];
				value_aux = value_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t];

				printf("%d %d %d %f %f \n", t, i1, i2, beta_aux, value_aux);
			}
		}
	}

	//theoretical, prob of all y_t=1, theoretical; beta_acum is the empirical one
	beta_ini = beta_vec[x1_0 + grid_size_x1 * x2_0 + grid_size_x1 * grid_size_x2 * 0];

	beta_aver_HMM = 0;
	beta_aver_variational = 0;
	beta2_aver_HMM = 0;
	beta2_aver_variational = 0;
	for (rep = 0; rep < repetitions; rep++) {
		//empirical beta in one trajectory, or prob(y=y_d), computed from the HMM agent (note averages over y are already implement per trajectory)
		beta_acum_HMM = trajectory_beta(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1, 
			p_y_food, p_y_no_food, beta_vec, x1_0, x2_0, prob_switch);

		//empirical beta in one trajectory, or prob(y=y_d), computed from the variational agent
		beta_acum_variational = trajectory_value(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1,
			p_y_food, p_y_no_food, value_vec, x1_0, x2_0, alpha, prob_switch);

		beta_aver_HMM = beta_aver_HMM + beta_acum_HMM;
		beta_aver_variational = beta_aver_variational + beta_acum_variational;

		beta2_aver_HMM = beta2_aver_HMM + pow(beta_acum_HMM ,2);
		beta2_aver_variational = beta2_aver_variational + pow(beta_acum_variational, 2);
	}

	beta_aver_HMM = beta_aver_HMM / repetitions;
	beta_aver_variational = beta_aver_variational / repetitions;

	beta2_aver_HMM = beta2_aver_HMM / repetitions;
	beta2_aver_variational = beta2_aver_variational / repetitions;

	error_beta_aver_HMM = sqrt(abs(beta2_aver_HMM - pow(beta_aver_HMM, 2)) / repetitions);
	error_beta_aver_variational = sqrt(abs(beta2_aver_variational - pow(beta_aver_variational, 2)) / repetitions);

	printf("%f %f %f %f %f %f %f \n", beta_ini, beta_aver_HMM, error_beta_aver_HMM, beta_ini / beta_aver_HMM,
		-log10(beta_ini), -log10(beta_aver_HMM), -log10(error_beta_aver_HMM));

	printf("%f %f %f %f %f %f %f \n", beta_ini, beta_aver_variational, error_beta_aver_variational, beta_ini / beta_aver_variational,
		-log10(beta_ini), -log10(beta_aver_variational), -log10(error_beta_aver_variational));


	//moving alpha
	d_alpha = (alpha_end - alpha_ini) / num_alpha;
	for (i = 0; i <= num_alpha; i++) {

		alpha0 = alpha_ini + i * d_alpha;  

		beta_vec = beta(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1, p_y_food, p_y_no_food, prob_switch);

		value_vec = value_function(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1, p_y_food, p_y_no_food,
			alpha0, prob_switch);

		//theoretical, prob of all y_t=1, theoretical; beta_acum is the empirical one
		beta_ini = beta_vec[x1_0 + grid_size_x1 * x2_0 + grid_size_x1 * grid_size_x2 * 0];

		beta_aver_HMM = 0;
		beta_aver_variational = 0;
		beta2_aver_HMM = 0;
		beta2_aver_variational = 0;
		for (rep = 0; rep < repetitions; rep++) {
			//empirical beta in one trajectory, or prob(y=y_d), computed from the HMM agent (note averages over y are already implement per trajectory)
			beta_acum_HMM = trajectory_beta(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1,
				p_y_food, p_y_no_food, beta_vec, x1_0, x2_0, prob_switch);

			//empirical beta in one trajectory, or prob(y=y_d), computed from the variational agent
			beta_acum_variational = trajectory_value(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1,
				p_y_food, p_y_no_food, value_vec, x1_0, x2_0, alpha0, prob_switch);

			beta_aver_HMM = beta_aver_HMM + beta_acum_HMM;
			beta_aver_variational = beta_aver_variational + beta_acum_variational;

			beta2_aver_HMM = beta2_aver_HMM + pow(beta_acum_HMM, 2);
			beta2_aver_variational = beta2_aver_variational + pow(beta_acum_variational, 2);
		}

		beta_aver_HMM = beta_aver_HMM / repetitions;
		beta_aver_variational = beta_aver_variational / repetitions;

		beta2_aver_HMM = beta2_aver_HMM / repetitions;
		beta2_aver_variational = beta2_aver_variational / repetitions;

		error_beta_aver_HMM = sqrt(abs(beta2_aver_HMM - pow(beta_aver_HMM, 2)) / repetitions);
		error_beta_aver_variational = sqrt(abs(beta2_aver_variational - pow(beta_aver_variational, 2)) / repetitions);

		printf("%f %f %f %f %f %f %f %f %f \n", alpha0, beta_ini,
			beta_aver_HMM, -log(beta_ini),
			-log(beta_aver_HMM), error_beta_aver_HMM,
			beta_aver_variational, -log(beta_aver_variational), error_beta_aver_variational);

		fprintf(beta_values_alpha, "%f %f %f %f %f %f \n", alpha0, beta_ini,
			beta_aver_HMM, error_beta_aver_HMM,
			beta_aver_variational, error_beta_aver_variational);
	}

	//moving prob_switch
	d_prob_switch = (prob_switch_end - prob_switch_ini) / num_prob_switch;
	for (i = 0; i <= num_prob_switch; i++) {

		prob_switch0 = prob_switch_ini + i * d_prob_switch;

		beta_vec = beta(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1, p_y_food, p_y_no_food, prob_switch0);

		value_vec = value_function(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1, p_y_food, p_y_no_food,
			alpha, prob_switch0);

		//theoretical, prob of all y_t=1, theoretical; beta_acum is the empirical one
		beta_ini = beta_vec[x1_0 + grid_size_x1 * x2_0 + grid_size_x1 * grid_size_x2 * 0];

		beta_aver_HMM = 0;
		beta_aver_variational = 0;
		beta2_aver_HMM = 0;
		beta2_aver_variational = 0;
		for (rep = 0; rep < repetitions; rep++) {
			//empirical beta in one trajectory, or prob(y=y_d), computed from the HMM agent (note averages over y are already implement per trajectory)
			beta_acum_HMM = trajectory_beta(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1,
				p_y_food, p_y_no_food, beta_vec, x1_0, x2_0, prob_switch0);

			//empirical beta in one trajectory, or prob(y=y_d), computed from the variational agent
			beta_acum_variational = trajectory_value(grid_size_x1, grid_size_x2, T, x1_food, x2_food, middle_wall_x2, door_x1,
				p_y_food, p_y_no_food, value_vec, x1_0, x2_0, alpha, prob_switch0);

			beta_aver_HMM = beta_aver_HMM + beta_acum_HMM;
			beta_aver_variational = beta_aver_variational + beta_acum_variational;

			beta2_aver_HMM = beta2_aver_HMM + pow(beta_acum_HMM, 2);
			beta2_aver_variational = beta2_aver_variational + pow(beta_acum_variational, 2);
		}

		beta_aver_HMM = beta_aver_HMM / repetitions;
		beta_aver_variational = beta_aver_variational / repetitions;

		beta2_aver_HMM = beta2_aver_HMM / repetitions;
		beta2_aver_variational = beta2_aver_variational / repetitions;

		error_beta_aver_HMM = sqrt(abs(beta2_aver_HMM - pow(beta_aver_HMM, 2)) / repetitions);
		error_beta_aver_variational = sqrt(abs(beta2_aver_variational - pow(beta_aver_variational, 2)) / repetitions);

		printf("%f %f %f %f %f %f \n", prob_switch0, beta_ini,
			beta_aver_HMM, error_beta_aver_HMM,
			beta_aver_variational, error_beta_aver_variational);

		fprintf(beta_values_prob_switch, "%f %f %f %f %f %f \n", prob_switch0, beta_ini,
			beta_aver_HMM, error_beta_aver_HMM,
			beta_aver_variational, error_beta_aver_variational);
	}

	fclose(beta_values_alpha);
	fclose(beta_values_prob_switch);
	fclose(bellman_values);
	
	return 0;
}



//Functions

//betas, computing prob p(y_{t:T}=y_d);
double* beta(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, 
	int door_x1, double p_y_food, double p_y_no_food, double prob_switch) {

	int    t, i1, i2;

	double beta_max, beta_aux, beta_aver;

	double* beta_vec;

	beta_vec = (double*)malloc(grid_size_x1 * grid_size_x2 * T * sizeof(double)); //time-dependent beta_values (multiplicative)

	//backpropagated beta values using multiplicative beta max operator

	//for t=T-1, boundary conditions betas
	for (i1 = 0; i1 < grid_size_x1; i1++) {
		for (i2 = 0; i2 < grid_size_x2; i2++) {
			if (i1 == food_x1 && i2 == food_x2) {
				beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (T - 1)] = p_y_food;
			}
			else if (i1 == food_x1 - 1 && i2 == food_x2) {//dangerous cell right on the left of the food source
				beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (T - 1)] = 1e-200 * p_y_no_food;
			}
			else {
				beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (T - 1)] = p_y_no_food;
			}	
		}
	}

	//for t<T-1
	for (t = T-2; t >= 0; t--) { //backwards from t=T-1 to t=0
		for (i1 = 0; i1 < grid_size_x1; i1++) {
			for (i2 = 0; i2 < grid_size_x2; i2++) {

				//computing the best transition; transitions are deterministic in this example
				// 5 possible actions:
				
				

				//make sure that the food source is at a location where all 5 movements are possible
			
				//action 0, Not moving
				beta_max = beta_vec[i1 + 0 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)];
				//action 1, moving North
				if (i2 + 1 == grid_size_x2 || (i2 + 1 == middle_wall_x2 && i1 != door_x1) ) { //if next to the north wall, or middle wall but not on the door
					//not moving in this case because next to wall
					beta_aux = beta_vec[i1 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)]; //not moving against the right wall
				}
				else {
					//move in this case because away from wall
					beta_aux = beta_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t + 1)]; //moving towards the right wall
				}
				if (beta_aux > beta_max) {
					beta_max = beta_aux;
				}
				//action 2, moving East
				if (i1 + 1 == grid_size_x1 || (i1 == door_x1 && i2 == middle_wall_x2) ) { //if next to the right wall, or right at the door
					beta_aux = beta_vec[(i1 + 0) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)]; //not moving against the right wall
				}
				else {
					beta_aux = beta_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)]; //moving towards the right wall
				}
				if (beta_aux > beta_max) {
					beta_max = beta_aux;
				}
				//action 3, moving South
				if (i2 - 1 < 0 || (i2 - 1 == middle_wall_x2 && i1 != door_x1)) {
					beta_aux = beta_vec[i1 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				else {
					beta_aux = beta_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				if (beta_aux > beta_max) {
					beta_max = beta_aux;
				}
				//action 4, moving West
				if (i1 - 1 < 0 || (i1 == door_x1 && i2 == middle_wall_x2)) {
					beta_aux = beta_vec[(i1 + 0) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				else {
					beta_aux = beta_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				if (beta_aux > beta_max) {
					beta_max = beta_aux;
				}
				
				if (i1 == food_x1 && i2 == food_x2) {
					//no choice of action if in the food source, diffuse randomly, undoing the above calculation
					beta_max = (1. - prob_switch) * beta_vec[(i1 + 0) + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 + 0) + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 + 1) + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 + 0) + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 - 1) + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				
				//final update
				//we assume that if in the food source, there is a small probability of not staying if this action is chosen
				if (i1 == food_x1 && i2 == food_x2) {
					beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t] = p_y_food * beta_max;
				}
				else if (i1 == food_x1 - 1 && i2 == food_x2) {//dangerous cell right on the left of the food source
					beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t] = 1e-200 * p_y_no_food * beta_max;
				}
				else {
					beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t] = p_y_no_food * beta_max;
				}
				
			}
		}	
	}

	return beta_vec;

}


//generating one trajectory from beta_vec
double trajectory_beta(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1, 
	double p_y_food, double p_y_no_food, double* beta_vec, int x1_0, int x2_0, double prob_switch) {

	int    t, i1, i2, x1_new, x2_new;
	double beta_max, beta_acum, beta_ini;
	double r, prob1, prob2, prob3, prob4, prob5;

	double *trajectory_vec;

	trajectory_vec = (double*)malloc(2 * T * sizeof(double));

	//random generator
	std::random_device rd;
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with constant seed_r
	std::uniform_real_distribution<> uniform(0.0, 1.0);

	//for t=0
	i1 = x1_0;
	i2 = x2_0;
	if (i1 == food_x1 && i2 == food_x2) {
		beta_acum = p_y_food;
	}
	else if (i1 == food_x1 - 1 && i2 == food_x2) {
		beta_acum = 1e-200 * p_y_no_food;
	}
	else {
		beta_acum = p_y_no_food;
	}

	for (t = 1; t < T; t++) { //trajectory from t=1 to t=T
		//computing optimal actions, 5 options
		//no motion
		x1_new = i1;
		x2_new = i2;
		
		//not moving
		beta_max = beta_vec[i1 + 0 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t - 0)];
		//moving North
		if (i2 + 1 < grid_size_x2 && (i2 + 1 != middle_wall_x2 || i1 == door_x1) &&
			beta_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1;
			x2_new = i2 + 1;
			beta_max = beta_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
		}
		//moving East
		if (i1 + 1 < grid_size_x1 && (i2 != middle_wall_x2) &&
			beta_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1 + 1;
			x2_new = i2;
			beta_max = beta_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
		}
		//moving South
		if (i2 - 1 >= 0 && (i2 - 1 != middle_wall_x2 || i1 == door_x1) &&
			beta_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1;
			x2_new = i2 - 1;
			beta_max = beta_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
		}
		//moving West
		if (i1 - 1 >= 0 && (i2 != middle_wall_x2) &&
			beta_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1 - 1;
			x2_new = i2;
			beta_max = beta_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
		}

		//if in the food then move randomly in any direction, no possibility of choice
		//no possibility of choosing action in this case
		if (i1 == food_x1 && i2 == food_x2) {

			//cumulative probs
			prob1 = 1. - prob_switch;
			prob2 = prob1 + prob_switch / 4.;
			prob3 = prob2 + prob_switch / 4.;
			prob4 = prob3 + prob_switch / 4.;
			prob5 = prob4 + prob_switch / 4.;

			//sampling action
			r = uniform(gen);

			//default
			x1_new = i1;
			x2_new = i2;
			if (r < prob1) {
				x1_new = i1;
				x2_new = i2;
			}
			if (r < prob2 && r >= prob1) {
				x1_new = i1;
				x2_new = i2 + 1;
			}
			if (r < prob3 && r >= prob2) {
				x1_new = i1 + 1;
				x2_new = i2;
			}
			if (r < prob4 && r >= prob3) {
				x1_new = i1;
				x2_new = i2 - 1;
			}
			if (r <= prob5 && r >= prob4) {
				x1_new = i1 - 1;
				x2_new = i2;
			}
		}

		//final update
		i1 = x1_new;
		i2 = x2_new;

		//probability p(y=y_d)
		if (i1 == food_x1 && i2 == food_x2) {
			beta_acum = beta_acum * p_y_food;
		}
		else if (i1 == food_x1 - 1 && i2 == food_x2) {//dangerous cell right on the left of the food source
			beta_acum = beta_acum * 1e-200 * p_y_no_food;
		}
		else {
			beta_acum = beta_acum * p_y_no_food;
		}
		

		trajectory_vec[2 * t] = i1;
		trajectory_vec[2 * t + 1] = i2;
	}

	//prob of all y_t=1, theoretical; beta_acum is the empirical one
	beta_ini = beta_vec[x1_0 + grid_size_x1 * x2_0 + grid_size_x1 * grid_size_x2 * 0]; 

	//printf("%f %f %f %f %f \n", beta_ini, beta_acum, beta_ini/beta_acum, -log10(beta_ini), -log10(beta_acum));

	return beta_acum;

}


//value function, Bellman recursion, variational inference;
double* value_function(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2,
	int door_x1, double p_y_food, double p_y_no_food, double alpha, double prob_switch) {

	int    t, i1, i2;

	double beta1, beta2, beta3, beta4, beta5, beta_aver, log_sum_exp; //q values, each for each action

	double* beta_vec; //value function

	beta_vec = (double*)malloc(grid_size_x1 * grid_size_x2 * T * sizeof(double)); //time-dependent beta_values (multiplicative)

	//backpropagated beta values using optimal value Bellman recursion

	//for t=T-1, boundary conditions betas
	for (i1 = 0; i1 < grid_size_x1; i1++) {
		for (i2 = 0; i2 < grid_size_x2; i2++) {
			if (i1 == food_x1 && i2 == food_x2) {
				beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (T - 1)] = alpha * log(p_y_food) + log(5.);
			}
			else if (i1 == food_x1 - 1 && i2 == food_x2) {//dangerous cell right on the left of the food source
				beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (T - 1)] = alpha * log(1e-200 * p_y_no_food) + log(5.);
			}
			else {
				beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (T - 1)] = alpha * log(p_y_no_food) + log(5.);
			}
		}
	}

	//for t<T-1
	for (t = T - 2; t >= 0; t--) { //backwards from t=T-1 to t=0
		for (i1 = 0; i1 < grid_size_x1; i1++) {
			for (i2 = 0; i2 < grid_size_x2; i2++) {

				//computing the best transition; transitions are deterministic in this example
				// 5 possible actions:
				//action 0, Not moving
				beta1 = beta_vec[i1 + 0 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)];	
				//action 1, moving North
				if (i2 + 1 == grid_size_x2 || (i2 + 1 == middle_wall_x2 && i1 != door_x1)) { //if next to the north wall, or middle wall but not on the door
					//not moving in this case because next to wall
					beta2 = beta_vec[i1 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)]; //not moving against the right wall
				}
				else {
					//move in this case because away from wall
					beta2 = beta_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t + 1)]; //moving towards the right wall
				}
				//action 2, moving East
				if (i1 + 1 == grid_size_x1 || (i1 == door_x1 && i2 == middle_wall_x2)) { //if next to the right wall, or right at the door
					beta3 = beta_vec[(i1 + 0) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)]; //not moving against the right wall
				}
				else {
					beta3 = beta_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)]; //moving towards the right wall
				}
				//action 3, moving South
				if (i2 - 1 < 0 || (i2 - 1 == middle_wall_x2 && i1 != door_x1)) {
					beta4 = beta_vec[i1 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				else {
					beta4 = beta_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				//action 4, moving West
				if (i1 - 1 < 0 || (i1 == door_x1 && i2 == middle_wall_x2)) {
					beta5 = beta_vec[(i1 + 0) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)];
				}
				else {
					beta5 = beta_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t + 1)];
				}

				log_sum_exp = log(exp(beta1) + exp(beta2) + exp(beta3) + exp(beta4) + exp(beta5));

				//if in the food source, no choice of action, so random diffusion
				//make sure that the food source is at a location where all 5 movements are possible
				if (i1 == food_x1 && i2 == food_x2) {
					beta1 = (1. - prob_switch) * beta_vec[(i1 + 0) + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 + 0) + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 + 1) + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 + 0) + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t + 1)]
						+ prob_switch / 4. * beta_vec[(i1 - 1) + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t + 1)];

					log_sum_exp = log(exp(beta1)); //no action allowed
				}

				//final update
				//we assume that if in the food source, there is a small probability of not staying if this action is chosen
				if (i1 == food_x1 && i2 == food_x2) {
					beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t] = alpha * log(p_y_food) + log_sum_exp;
				}
				else if (i1 == food_x1 - 1 && i2 == food_x2) {//dangerous cell right on the left of the food source
					beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t] = alpha * log(1e-200 * p_y_no_food) + log_sum_exp;
				}
				else {
					beta_vec[i1 + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * t] = alpha * log(p_y_no_food) + log_sum_exp;
				}
			}
		}
	}

	return beta_vec;
}

//generating one trajectory from value_vec
double trajectory_value(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1,
	double p_y_food, double p_y_no_food, double* value_vec, int x1_0, int x2_0, double alpha, double prob_switch) {

	int    t, i1, i2, x1_new, x2_new;
	double beta_acum, r;

	double beta1, beta2, beta3, beta4, beta5, sum_exp; //q values, each for each action
	double prob1, prob2, prob3, prob4, prob5; //prob each for each action

	double* trajectory_vec;

	trajectory_vec = (double*)malloc(2 * T * sizeof(double));

	//random generator
	std::random_device rd;
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with constant seed_r
	std::uniform_real_distribution<> uniform(0.0, 1.0);

	//for t=0
	i1 = x1_0;
	i2 = x2_0;
	if (i1 == food_x1 && i2 == food_x2) {
		beta_acum = p_y_food;
	}
	else if (i1 == food_x1 - 1 && i2 == food_x2) {
		beta_acum = 1e-200 * p_y_no_food;
	}
	else {
		beta_acum = p_y_no_food;
	}

	for (t = 1; t < T; t++) { //trajectory from t=1 to t=T
		
		
		//if in the food, then move randomly in any direction
		if (i1 == food_x1 && i2 == food_x2) {

			//cumulative probs
			prob1 = 1. - prob_switch;
			prob2 = prob1 + prob_switch / 4.;
			prob3 = prob2 + prob_switch / 4.;
			prob4 = prob3 + prob_switch / 4.;
			prob5 = prob4 + prob_switch / 4.;
		}
		else {
			//computing probability of actions, 5 options

			//no motion
			//make sure that the food source is at a location where all 5 movements are possible
			beta1 = value_vec[i1 + 0 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t - 0)];
			//moving North
			if (i2 + 1 < grid_size_x2 && (i2 + 1 != middle_wall_x2 || i1 == door_x1))
			{
				beta2 = value_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta2 = beta1;
			}
			//moving East
			if (i1 + 1 < grid_size_x1 && (i2 != middle_wall_x2))
			{
				beta3 = value_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta3 = beta1;
			}
			//moving South
			if (i2 - 1 >= 0 && (i2 - 1 != middle_wall_x2 || i1 == door_x1))
			{
				beta4 = value_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta4 = beta1;
			}
			//moving West
			if (i1 - 1 >= 0 && (i2 != middle_wall_x2))
			{
				beta5 = value_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta5 = beta1;
			}

			sum_exp = exp(beta1) + exp(beta2) + exp(beta3) + exp(beta4) + exp(beta5);

			//cumulative probs
			prob1 = exp(beta1) / sum_exp;
			prob2 = prob1 + exp(beta2) / sum_exp;
			prob3 = prob2 + exp(beta3) / sum_exp;
			prob4 = prob3 + exp(beta4) / sum_exp;
			prob5 = prob4 + exp(beta5) / sum_exp;

		}

		//sampling action
		r = uniform(gen);

		//default
		x1_new = i1;
		x2_new = i2;
		if (r < prob1) {
			x1_new = i1;
			x2_new = i2;
		}
		if (r < prob2 && r >= prob1) {
			if (i2 + 1 < grid_size_x2 && (i2 + 1 != middle_wall_x2 || i1 == door_x1)) {
				x1_new = i1;
				x2_new = i2 + 1;
			}
		}
		if (r < prob3 && r >= prob2) {
			if (i1 + 1 < grid_size_x1 && (i2 != middle_wall_x2)) {
				x1_new = i1 + 1;
				x2_new = i2;
			}
		}
		if (r < prob4 && r >= prob3) {
			if (i2 - 1 >= 0 && (i2 - 1 != middle_wall_x2 || i1 == door_x1)) {
				x1_new = i1;
				x2_new = i2 - 1;
			}
		}
		if (r <= prob5 && r >= prob4) {
			if (i1 - 1 >= 0 && (i2 != middle_wall_x2)) {
				x1_new = i1 - 1;
				x2_new = i2;
			}
		}

		//final update
		i1 = x1_new;
		i2 = x2_new;

		//probability
		if (i1 == food_x1 && i2 == food_x2) {
			beta_acum = beta_acum * p_y_food;
		}
		else if (i1 == food_x1 - 1 && i2 == food_x2) {
			beta_acum = beta_acum * 1e-200 * p_y_no_food;
		}
		else {
			beta_acum = beta_acum * p_y_no_food;
		}

		trajectory_vec[2 * t] = i1;
		trajectory_vec[2 * t + 1] = i2;

		//printf("%d %d %d \n", t, i1, i2);
	}

	//printf("%f %f \n", beta_acum, -log10(beta_acum));

	return beta_acum;

}


//generating one trajectory from beta_vec
int* trajectory_beta_samples(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1,
	double p_y_food, double p_y_no_food, double* beta_vec, int x1_0, int x2_0, double prob_switch) {

	int    t, i1, i2, x1_new, x2_new;
	double beta_max, beta_acum, beta_ini;
	double r, prob1, prob2, prob3, prob4, prob5;

	int* trajectory_vec;

	trajectory_vec = (int*)malloc(2 * T * sizeof(int));

	//random generator
	std::random_device rd;
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with constant seed_r
	std::uniform_real_distribution<> uniform(0.0, 1.0);

	//for t=0
	i1 = x1_0;
	i2 = x2_0;
	trajectory_vec[0] = i1;
	trajectory_vec[1] = i2;

	if (i1 == food_x1 && i2 == food_x2) {
		beta_acum = p_y_food;
	}
	else if (i1 == food_x1 - 1 && i2 == food_x2) {//dangerous cell right on the left of the food source
		beta_acum = 1e-200 * p_y_no_food;
	}
	else {
		beta_acum = p_y_no_food;
	}

	for (t = 1; t < T; t++) { //trajectory from t=1 to t=T
		//computing optimal actions, 5 options
		//no motion
		x1_new = i1;
		x2_new = i2;

		//not moving
		beta_max = beta_vec[i1 + 0 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t - 0)];
		//moving North
		if (i2 + 1 < grid_size_x2 && (i2 + 1 != middle_wall_x2 || i1 == door_x1) &&
			beta_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1;
			x2_new = i2 + 1;
			beta_max = beta_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
		}
		//moving East
		if (i1 + 1 < grid_size_x1 && (i2 != middle_wall_x2) &&
			beta_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1 + 1;
			x2_new = i2;
			beta_max = beta_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
		}
		//moving South
		if (i2 - 1 >= 0 && (i2 - 1 != middle_wall_x2 || i1 == door_x1) &&
			beta_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1;
			x2_new = i2 - 1;
			beta_max = beta_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
		}
		//moving West
		if (i1 - 1 >= 0 && (i2 != middle_wall_x2) &&
			beta_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)] > beta_max)
		{
			x1_new = i1 - 1;
			x2_new = i2;
			beta_max = beta_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
		}

		//if in the food then move randomly in any direction, no possibility of choice
		//no possibility of choosing action in this case
		if (i1 == food_x1 && i2 == food_x2) {

			//cumulative probs
			prob1 = 1. - prob_switch;
			prob2 = prob1 + prob_switch / 4.;
			prob3 = prob2 + prob_switch / 4.;
			prob4 = prob3 + prob_switch / 4.;
			prob5 = prob4 + prob_switch / 4.;

			//sampling action
			r = uniform(gen);

			//default
			x1_new = i1;
			x2_new = i2;
			if (r < prob1) {
				x1_new = i1;
				x2_new = i2;
			}
			if (r < prob2 && r >= prob1) {
				x1_new = i1;
				x2_new = i2 + 1;
			}
			if (r < prob3 && r >= prob2) {
				x1_new = i1 + 1;
				x2_new = i2;
			}
			if (r < prob4 && r >= prob3) {
				x1_new = i1;
				x2_new = i2 - 1;
			}
			if (r <= prob5 && r >= prob4) {
				x1_new = i1 - 1;
				x2_new = i2;
			}
		}

		//final update
		i1 = x1_new;
		i2 = x2_new;

		//probability p(y=y_d)
		if (i1 == food_x1 && i2 == food_x2) {
			beta_acum = beta_acum * p_y_food;
		}
		else if (i1 == food_x1 - 1 && i2 == food_x2) {//dangerous cell right on the left of the food source
			beta_acum = beta_acum * 1e-200 * p_y_no_food;
		}
		else {
			beta_acum = beta_acum * p_y_no_food;
		}


		trajectory_vec[2 * t] = i1;
		trajectory_vec[2 * t + 1] = i2;
	}

	//prob of all y_t=1, theoretical; beta_acum is the empirical one
	beta_ini = beta_vec[x1_0 + grid_size_x1 * x2_0 + grid_size_x1 * grid_size_x2 * 0];

	printf("%f %f %f %f %f %f \n", beta_ini, beta_acum, beta_ini / beta_acum, -log10(beta_ini), -log10(beta_acum), 1. * r);

	return trajectory_vec;

}


//generating one trajectory from value_vec
int* trajectory_value_samples(int grid_size_x1, int grid_size_x2, int T, int food_x1, int food_x2, int middle_wall_x2, int door_x1,
	double p_y_food, double p_y_no_food, double* value_vec, int x1_0, int x2_0, double alpha, double prob_switch) {

	int    t, i1, i2, x1_new, x2_new;
	double beta_acum, r;

	double beta1, beta2, beta3, beta4, beta5, sum_exp; //q values, each for each action
	double prob1, prob2, prob3, prob4, prob5; //prob each for each action

	int* trajectory_vec;

	trajectory_vec = (int*)malloc((2 * T + 1) * sizeof(int));

	//random generator
	std::random_device rd;
	std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with constant seed_r
	std::uniform_real_distribution<> uniform(0.0, 1.0);

	//for t=0
	i1 = x1_0;
	i2 = x2_0;
	trajectory_vec[0] = i1;
	trajectory_vec[1] = i2;

	if (i1 == food_x1 && i2 == food_x2) {
		beta_acum = p_y_food;
	}
	else if (i1 == food_x1 - 1 && i2 == food_x2) {
		beta_acum = 1e-200 * p_y_no_food;
	}
	else {
		beta_acum = p_y_no_food;
	}

	for (t = 1; t < T; t++) { //trajectory from t=1 to t=T
		//computing probability of actions, 5 options


		//if in the food, then move randomly in any direction
		if (i1 == food_x1 && i2 == food_x2) {

			//cumulative probs
			prob1 = 1. - prob_switch;
			prob2 = prob1 + prob_switch / 4.;
			prob3 = prob2 + prob_switch / 4.;
			prob4 = prob3 + prob_switch / 4.;
			prob5 = prob4 + prob_switch / 4.;
		}
		else {

			//no motion
			beta1 = value_vec[i1 + 0 + grid_size_x1 * (i2 + 0) + grid_size_x1 * grid_size_x2 * (t - 0)];
			//moving North
			if (i2 + 1 < grid_size_x2 && (i2 + 1 != middle_wall_x2 || i1 == door_x1))
			{
				beta2 = value_vec[i1 + grid_size_x1 * (i2 + 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta2 = beta1;
			}
			//moving East
			if (i1 + 1 < grid_size_x1 && (i2 != middle_wall_x2))
			{
				beta3 = value_vec[(i1 + 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta3 = beta1;
			}
			//moving South
			if (i2 - 1 >= 0 && (i2 - 1 != middle_wall_x2 || i1 == door_x1))
			{
				beta4 = value_vec[i1 + grid_size_x1 * (i2 - 1) + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta4 = beta1;
			}
			//moving West
			if (i1 - 1 >= 0 && (i2 != middle_wall_x2))
			{
				beta5 = value_vec[(i1 - 1) + grid_size_x1 * i2 + grid_size_x1 * grid_size_x2 * (t - 0)];
			}
			else {
				beta5 = beta1;
			}

			sum_exp = exp(beta1) + exp(beta2) + exp(beta3) + exp(beta4) + exp(beta5);

			//cumulative probs
			prob1 = exp(beta1) / sum_exp;
			prob2 = prob1 + exp(beta2) / sum_exp;
			prob3 = prob2 + exp(beta3) / sum_exp;
			prob4 = prob3 + exp(beta4) / sum_exp;
			prob5 = prob4 + exp(beta5) / sum_exp;

		}

		//sampling action
		r = uniform(gen);

		//default
		x1_new = i1;
		x2_new = i2;
		if (r < prob1) {
			x1_new = i1;
			x2_new = i2;
		}
		if (r < prob2 && r >= prob1) {
			if (i2 + 1 < grid_size_x2 && (i2 + 1 != middle_wall_x2 || i1 == door_x1)) {
				x1_new = i1;
				x2_new = i2 + 1;
			}
		}
		if (r < prob3 && r >= prob2) {
			if (i1 + 1 < grid_size_x1 && (i2 != middle_wall_x2)) {
				x1_new = i1 + 1;
				x2_new = i2;
			}
		}
		if (r < prob4 && r >= prob3) {
			if (i2 - 1 >= 0 && (i2 - 1 != middle_wall_x2 || i1 == door_x1)) {
				x1_new = i1;
				x2_new = i2 - 1;
			}
		}
		if (r <= prob5 && r >= prob4) {
			if (i1 - 1 >= 0 && (i2 != middle_wall_x2)) {
				x1_new = i1 - 1;
				x2_new = i2;
			}
		}

		//final update
		i1 = x1_new;
		i2 = x2_new;

		//probability
		if (i1 == food_x1 && i2 == food_x2) {
			beta_acum = beta_acum * p_y_food;
		}
		else if (i1 == food_x1 - 1 && i2 == food_x2) {
			beta_acum = beta_acum * 1e-200 * p_y_no_food;
		}
		else {
			beta_acum = beta_acum * p_y_no_food;
		}

		trajectory_vec[2 * t] = i1;
		trajectory_vec[2 * t + 1] = i2;

		//printf("%d %d %d \n", t, i1, i2);
	}

	return trajectory_vec;

}
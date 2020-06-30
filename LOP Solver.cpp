// LOP Solver.cpp : Defines the entry point for the console application.
//

#include <ilcplex/ilocplex.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include <windows.h>
#include <vector>
#include <random>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
using namespace std;

ILOSTLBEGIN
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<IloRangeArray> ConstraintMatrix;

typedef vector<int > IntVector;
typedef vector< vector<int > > TwoIntMatrix;
typedef vector< vector< vector<int > > > ThreeIntMatrix;
typedef vector<double > DoubleVector;
typedef vector< vector<double > > TwoDMatrix;
typedef vector< vector< vector<double > > > ThreeDMatrix;

void lp_out(IloEnv &env, IloCplex cplex, ConstraintMatrix c);
// Functions for reading data
TwoDMatrix readdata(char datafile[], int*);
TwoIntMatrix readmissingdata(char missingfile[], int* error);
configuration readconfig(char configname[], int * error);
// Functions for creating a starting solution
ThreeIntMatrix start(int, TwoDMatrix);
// Functions for setting up the LP Model
IloNumVarArray LP_vars_cr(IloEnv &env, ThreeIntMatrix Preferences);
IloRange LP_Sum(IloEnv &env, IloNumVarArray vars, ThreeIntMatrix Preferences);
ConstraintMatrix LP_Cons(IloEnv &env, IloNumVarArray vars, ThreeIntMatrix Preferences, TwoDMatrix P, int size);
IloModel LP_model(IloEnv &env, IloNumVarArray vars, IloObjective obj, IloRange sum, ConstraintMatrix c, int size);
// Functions for setting up the IP Model
NumVarMatrix IP_vars_cr(IloEnv &env, int size);
IloModel IP_model(IloEnv &env, NumVarMatrix IP_vars, IloNumVar target, IloObjective obj, int size);
IloModel IP_model(IloEnv &env, NumVarMatrix ip_vars, IloObjective obj, TwoDMatrix P, int size);
// Functions for the pricing problem
void IP_Update(IloEnv &env, IloCplex cplex_lp, ConstraintMatrix c_lp, IloRange sum, NumVarMatrix ip_vars,
	IloNumVar target_var, IloObjective ip_obj, int size);
ThreeIntMatrix exact_price(IloCplex cplex_ip, NumVarMatrix ip_vars, int size, double* ip_sol);
ThreeIntMatrix best_insertion(IloEnv &env, IloCplex cplex_LP, ConstraintMatrix c, IloRange sum, int size, double* ip_sol);
IntVector dominant_alternative(IntVector ranking, TwoDMatrix obj_coef, int size, IloEnv &env, IloCplex cplex_LP, ConstraintMatrix c);
ThreeIntMatrix time_limit_price(IloEnv &env, IloModel ip_model, IloCplex cplex_ip, NumVarMatrix ip_vars, int size);
// Outputfunctions
void outputpatterns(ThreeIntMatrix Preferences, int size);
// Function for finding FDIs
ThreeDMatrix Find_Valid_Inequality(ThreeDMatrix VIneq, IloEnv &env, IloCplex cplex_lp, ConstraintMatrix c_lp, IloRange sum,
	int size, int large_loop);
bool check_VIneq(ThreeDMatrix VIneq, TwoDMatrix point, int size);
bool check_Triangle(TwoDMatrix point, int size);
void output_VIneq(ThreeDMatrix VIneq, int size);

int main(int argc, char  * argv[])
{
	int i, j, k;
	int error = 0;
	int large_loop;
	int size = 0;
	int patterns = 0;
	int FDI = 0;
	int inpoly = 0; int outpoly = 0;
	int count_CPLEX = 0; int count_nontriangle = 0;
	double help_double;
	double randomfromunif;
	int random_number;
	int credit;
	int var_in_use; // Number of variables actually being used in the restricted master (taking into account overwrites).
	int var_replace; // Variable to be overwritten because it is not used.
	bool replace_flag;
	vector<int> var_counters;
	vector<IloNum> vals;
	//ofstream outputdetail;
	//ofstream outputLP;
	//ofstream outputstat;
	ofstream outputgeneral;
	//ofstream outputIP;
	//outputLP.open("outputLP.txt");
	//outputIP.open("OutputIP.txt");
	//outputdetail.open("Outputdetail.txt");
	//outputstat.open("Point Results.txt");

	boost::random::mt19937 generator;
	boost::random::uniform_01<double> distribution;

	TwoDMatrix P, Point;
	TwoIntMatrix Missingdata;
	ThreeIntMatrix Preferences; // A matrix which saves all preferences / columns generated pre and during column generation.
	ThreeDMatrix FDIs;
	ThreeDMatrix VIneq;

	// Time variables
	vector<double> price_time;
	vector<double> lp_value;


	int nr_tests = std::atoi(argv[5]);
	if (nr_tests == 0)
		cout << "NOTE; NUMBER OF REPETITIONS SET TO 0" << endl;

	cout << "Reading data" << endl;

	P = readdata(argv[1], &size);
	Missingdata = readmissingdata(argv[2], &error);
	credit = 15; // If LP variables have not been used in the restricted master this number of times, they can be discarded.
	Point = P;
	cout << "Starting Heuristic subfunction" << endl;
	Preferences = start(size, P);
	patterns = Preferences.size();
	var_in_use = patterns;
	// Write header for the outputfile.
	//outputdetail << "Pattern" << "\t\t" << "Value \t\t Time";
	// Initialize the output vectors for time needed and solution obtained.
	price_time.resize(patterns, 0);
	lp_value.resize(patterns, 0);
	var_counters.resize(patterns, credit);
	//cout << "Setting up LP problem" << endl;

	// Time Measurement
	time_t start_time;
	time(&start_time);

	// Setup of the LP Model, creation of all variables, constraints and objective function. Adding them to the model and setting 
	// solution parameters.
	IloEnv env;
	IloNumVarArray lp_vars = LP_vars_cr(env, Preferences);
	IloObjective lp_obj = (IloMinimize(env, lp_vars[0]));
	IloRange sum = LP_Sum(env, lp_vars, Preferences);
	ConstraintMatrix c = LP_Cons(env, lp_vars, Preferences, P, size);
	IloModel LP = LP_model(env, lp_vars, lp_obj, sum, c, size);

	IloCplex cplex_LP(LP);
	cplex_LP.setParam(IloCplex::NumericalEmphasis, 1);
	cplex_LP.setOut(env.getNullStream());

	// Setup of the IP Model, creation of all variables, constraints and objective function. Adding them to the model and setting 
	// solution parameters. Also creation of time-limited copy of the model.
	IloNumVar target_var(env, 1, 1, ILOINT);
	NumVarMatrix ip_vars = IP_vars_cr(env, size);
	IloObjective ip_obj = IloMaximize(env);
	IloModel IP = IP_model(env, ip_vars, ip_obj, P, size);

	IloCplex cplex_IP(IP);
	cplex_IP.setParam(IloCplex::NumericalEmphasis, 1);

	// Setup of the Column Generation
	// First, we declare some of the variables which will be used to store the value of current solutions
	bool improv;
	double lp_sol, ip_sol;
	double target;

	std::string str = string(argv[4]);
	boost::random::seed_seq seed(str.begin(), str.end());

	generator.seed(seed);
	std::cout << "Your seed produced: " << generator() << std::endl;

	// This loop is for sampling and rerunning using starting sets and hybrid method.
	for (large_loop = 0; large_loop < nr_tests; large_loop++)
	{
		std::cout << large_loop << endl;
		if (large_loop % 10000000 == 0)
			std::cout << large_loop << "\t" << inpoly << endl;
		for (i = 0; i < size; i++) /*Sampling of a datapoint */
		{
			for (j = i + 1; j < size; j++)
			{
				if (Missingdata[i][j] == 0)
				{
					randomfromunif = distribution(generator);
					help_double = randomfromunif;
					boost::math::beta_distribution<> dist(P[i][j] + 1, P[j][i] + 1);
					help_double = boost::math::quantile(dist, randomfromunif);
					c[i][j].setLB(help_double); Point[i][j] = help_double;
					c[j][i].setLB(1 - help_double); Point[j][i] = 1 - help_double;
				}
				else
				{
					c[i][j].setLB(0); Point[i][j] = 0;
					c[j][i].setLB(0); Point[j][i] = 0;
				}
			}
		}
		// We include a function which already separates the triangle (3-Dicycle) facets
		if (check_Triangle(Point, size) == 1)
		{
			lp_sol = 100;
			goto skip;
		}
		// First, we look through the already found FDIs, if any are violated, we can skip the CG and mark the data as outside the polytope.
		if (check_VIneq(VIneq, Point, size) == 1)
		{
			lp_sol = 10;
			goto skip;
		}
		// Now we start the main loop, it will run as long as the lp does not have a solution equal to 1 (with some tolerance) or until no more
		// improvements may be found.
		for (lp_sol = 1, ip_sol = 1; lp_sol >= 0.000001 && ip_sol >= 0.0000001;)
		{
			cplex_LP.solve();
			
			vals.resize(var_in_use + 1);
			for (i = 0; i <= var_in_use; i++)
			{
				vals[i] = cplex_LP.getValue(lp_vars[i]);
			}
			
			lp_sol = cplex_LP.getObjValue();
			// If no solution equal to one (with tolerance) is found, we go to the improvement step.
			if (lp_sol >= 0.000001)
			{				
				// Check which vars in LP are in use and change counters.
				for (i = 1; i <= var_in_use; i++) // vals[0] corresponds to z-variable.
				{
					if (vals[i] <= 0 && var_counters[i - 1] > 0)
						var_counters[i - 1] = var_counters[i - 1] - 1;
					if (vals[i] > 0)
						var_counters[i - 1] = credit;
				}
				
				// Find the first var to be overwritten
				for (i = 0, replace_flag = 0; i < var_in_use; i++)
				{
					if (var_counters[i] == 0)
					{
						var_replace = i + 1; // +1 due to offset of z-variable.
						replace_flag = 1;
						break;
					}
				}

				//Update the ip_problem
				IP_Update(env, cplex_LP, c, sum, ip_vars, target_var, ip_obj, size);

				// Find new columns
				ThreeIntMatrix add_col;

				// Option 0. Heuristic Approach (best insertion). If no improvement, move to option 1.
				ip_sol = 0;
				add_col = best_insertion(env, cplex_LP, c, sum, size, &ip_sol);

				// Option 1. Solve cplex_IP exactly and save the solution for a new column
				if (ip_sol < 0.0000001)
				{
					//cout << "Start CPLEX" << endl;
					add_col = exact_price(cplex_IP, ip_vars, size, &ip_sol);
					count_CPLEX++;
				}

				// Add these new columns to the linear program
				if (ip_sol >= 0.0000001)
				{
					//cout << "Add column" << endl;
					if (replace_flag == 0) // If no var to replace, add an extra one
					{
						for (i = 0; i < add_col.size(); i++)
						{
							var_counters.push_back(credit);
							lp_vars.add(IloNumVar(env, 0, 1, ILOFLOAT));
							patterns++;
							var_in_use++;
							sum.setLinearCoef(lp_vars[var_in_use], -1);
							for (j = 0; j < size; j++)
							{
								for (k = 0; k < size; k++)
								{
									if (add_col[i][j][k] == 1)
									{
										c[j][k].setLinearCoef(lp_vars[var_in_use], 1);
									}
								}
							}
						}
					}
					if (replace_flag == 1) // Else replace the var
					{
						patterns++;
						var_counters[var_replace - 1] = credit;
						for (j = 0; j < size; j++)
						{
							for (k = 0; k < size; k++)
							{
								if (add_col[0][j][k] == 1)
								{
									c[j][k].setLinearCoef(lp_vars[var_replace], 1);
								}
								else
								{
									c[j][k].setLinearCoef(lp_vars[var_replace], 0);
								}
							}
						}
					}
				}
			}
			// If CG will end in the next iteration without a solution, fire routine to find FDIs.
			if (lp_sol > 0.000001 && ip_sol < 0.0000001)
			{
				count_nontriangle++; VIneq = Find_Valid_Inequality(VIneq, env, cplex_LP, c, sum, size, large_loop);
			}
		}
	skip:
		if (lp_sol < 0.000001)
		{
			inpoly++;
		}
		else
		{
			outpoly++;
		}
	}

	time_t current_time;
	time(&current_time);
	outputgeneral.open(argv[3]);
	outputgeneral << "Time elapsed = " << difftime(current_time, start_time) << " s" << endl;
	outputgeneral << "Patterns created = " << patterns << endl;
	outputgeneral << "Active variables = " << var_in_use << endl;
	outputgeneral << "Points in polytope = " << inpoly << endl;
	outputgeneral << "Points outside of polytope = " << outpoly << endl;
	outputgeneral << "Times CPLEX called = " << count_CPLEX << endl;
	outputgeneral << "No Triangle violated, but outside polytope = " << count_nontriangle << endl;

	//outputpatterns(Preferences, size);
	//output_VIneq(VIneq, size);
	return 0;
}

void output_VIneq(ThreeDMatrix VIneq, int size)
{
	int i, j, k;
	ofstream Output_VIneq;

	Output_VIneq.open("OutputVIneq.txt");

	for (i = 0; i < VIneq.size(); i++)
	{
		Output_VIneq << "VIneq " << i << "\t from point " << VIneq[i][0][0] << endl;
		for (j = 0; j < size; j++)
		{
			Output_VIneq << "[ ";
			for (k = 0; k < size; k++)
			{
				if (j == 0 && k == 0)
				{
					Output_VIneq << "0, ";
				}
				else if (k != size)
				{
					Output_VIneq << VIneq[i][j][k] << ", ";
				}
				else
				{
					Output_VIneq << VIneq[i][j][k];
				}
			}
			Output_VIneq << "]" << endl;
		}
		Output_VIneq << endl;
	}

}

bool check_VIneq(ThreeDMatrix VIneq, TwoDMatrix point, int size)
{
	int i, j, k;
	double LHS;
	bool violated_VIneq = 0;

	for (i = 0; i < VIneq.size(); i++)
	{
		LHS = 0;
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				LHS = LHS + point[j][k] * VIneq[i][j][k];
			}
		}
		if (LHS > VIneq[i][size - 1][size - 1])
		{
			violated_VIneq = 1;
			break;
		}
	}
	return violated_VIneq;
}

bool check_Triangle(TwoDMatrix point, int size)
{
	int i, j, k;
	bool violated_VIneq = 0;

	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				if (point[i][j] + point[j][k] + point[k][i] > 2)
				{
					violated_VIneq = 1;
					break;
				}
			}
			if (violated_VIneq == 1)
				break;
		}
		if (violated_VIneq == 1)
			break;
	}
	return violated_VIneq;
}

ThreeDMatrix Find_Valid_Inequality(ThreeDMatrix VIneq, IloEnv &env, IloCplex cplex_lp, ConstraintMatrix c_lp, IloRange sum,
	int size, int large_loop)
{
	int i, j;
	TwoDMatrix new_VIneq;

	new_VIneq.resize(size);
	for (i = 0; i < size; i++)
	{
		new_VIneq[i].resize(size);
	}

	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			new_VIneq[i][j] = cplex_lp.getDual(c_lp[i][j]);
			if (new_VIneq[i][j] < 0.000001 && new_VIneq[i][j] > -0.000001)
			{
				new_VIneq[i][j] = 0;
			}
		}
	}
	new_VIneq[size - 1][size - 1] = cplex_lp.getDual(sum); // The target is saved to the final cell in the matrix
	new_VIneq[0][0] = large_loop; // The data point where the new valid inequality is found is saved

	VIneq.push_back(new_VIneq);

	return VIneq;
}

ThreeIntMatrix time_limit_price(IloEnv &env, IloModel ip_model, IloCplex cplex_ip, NumVarMatrix ip_vars, int size)
{
	/* This function solves the integer program with a time limit. If no solution is found which can add an extra column, then the
	IP is solved without a time limit. */
	int i, j;
	ThreeIntMatrix solution;
	IloCplex cplex_time(ip_model);
	cplex_time.setParam(IloCplex::TiLim, 1);

	solution.resize(1);
	solution[0].resize(size);
	for (i = 0; i < size; i++)
	{
		solution[0][i].resize(size);
	}

	cplex_time.solve();
	if (cplex_time.getObjValue() > 0)
	{
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				if (cplex_time.getValue(ip_vars[i][j]) > 0.99)
				{
					solution[0][i][j] = 1;
				}
				else
				{
					solution[0][i][j] = 0;
				}

			}
		}
	}
	else
	{
		cplex_ip.solve();
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				if (cplex_ip.getValue(ip_vars[i][j]) > 0.99)
				{
					solution[0][i][j] = 1;
				}
				else
				{
					solution[0][i][j] = 0;
				}

			}
		}
	}

	cplex_time.end();
	return solution;
}

ThreeIntMatrix exact_price(IloCplex cplex_ip, NumVarMatrix ip_vars, int size, double *ip_sol)
{
	int i, j;
	ThreeIntMatrix solution;

	solution.resize(1);
	solution[0].resize(size);
	for (i = 0; i < size; i++)
	{
		solution[0][i].resize(size);
	}
	cplex_ip.solve();
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			if (cplex_ip.getValue(ip_vars[i][j]) > 0.99)
			{
				solution[0][i][j] = 1;
			}
			else
			{
				solution[0][i][j] = 0;
			}

		}
	}
	*ip_sol = cplex_ip.getObjValue();
	//cout << *ip_sol << endl;
	//system("PAUSE");

	return solution;
}

ThreeIntMatrix best_insertion(IloEnv &env, IloCplex cplex_LP, ConstraintMatrix c, IloRange sum, int size, double* ip_sol)
{
	int i, j, k, l, m, n;
	double insertion_value;
	double solution_value;
	double best_insertion_value;
	double target_var = cplex_LP.getDual(sum);
	bool improve;
	bool reached_moving_alternative;
	int best_insertion_location;
	int best_insertion_alternative;
	TwoDMatrix obj_coef;
	obj_coef.resize(size);
	for (i = 0; i < size; i++)
	{
		obj_coef[i].resize(size);
		for (j = 0; j < size; j++)
		{
			obj_coef[i][j] = cplex_LP.getDual(c[i][j]);
			/*cout << obj_coef[i][j] << "\t";*/
		}
		/*cout << endl;*/
	}
	IntVector ranking;
	ranking.resize(size);
	ThreeIntMatrix output;

	for (i = 1; i < size; i++) // Start at i = 1, we always start with object 0, which is automatically placed first.
	{
		best_insertion_value = 0;
		best_insertion_location = 0;
		for (j = 0; j <= i; j++)
		{
			insertion_value = 0;
			for (k = 0; k < j; k++)
			{
				insertion_value = insertion_value + obj_coef[ranking[k]][i];
			}
			for (k = j; k < i; k++)
			{
				insertion_value = insertion_value + obj_coef[i][ranking[k]];
			}
			if (insertion_value >= best_insertion_value)
			{
				best_insertion_value = insertion_value;
				best_insertion_location = j;
			}
		}
		// Move all alternatives one place back, if the new alternative is placed before them.
		for (k = size - 2; k >= best_insertion_location; k--)
		{
			ranking[k + 1] = ranking[k];
		}
		ranking[best_insertion_location] = i;
	}


	for (improve = 1; improve == 1;)
	{
		improve = 0;
		/*for(k = 0; k < size; k++)
		{
		cout << ranking[k] << "\t";
		}
		cout << endl;*/
		for (i = 0; i < size && improve == 0; i++)
		{
			// First compute the current value associated with the object in a specific location.
			solution_value = 0;
			for (j = 0; j < i; j++)
			{
				solution_value = solution_value + obj_coef[ranking[j]][ranking[i]];
			}
			for (j = i + 1; j < size; j++)
			{
				solution_value = solution_value + obj_coef[ranking[i]][ranking[j]];
			}
			// We have now established the value of the current position of alternative ranking[i].
			// Now calculate the current best position.
			// We test size+1 positions, starting at j=0, i.e. the alternative is placed first.
			// j = 1 means placed between position 0 and position 1
			// j = size + 1 means being placed in last position
			for (j = 0; j < size + 1 && improve == 0; j++)
			{
				insertion_value = 0;
				for (k = 0; k < j; k++)
				{
					insertion_value = insertion_value + obj_coef[ranking[k]][ranking[i]];
				}
				for (k = j; k < size; k++)
				{
					insertion_value = insertion_value + obj_coef[ranking[i]][ranking[k]];
				}
				// Rearange ranking based on new location if better position is found.
				if (insertion_value > solution_value + 0.0000001)
				{
					improve = 1;
					/*cout << "i = " << i << "\t j = " << j << endl;*/
					best_insertion_alternative = ranking[i];
					if (i > j) // The alternative is moved to the front.
					{
						for (k = i; k > j; k--)
						{
							ranking[k] = ranking[k - 1];
						}
						ranking[j] = best_insertion_alternative;
					}
					if (j > i)
					{
						for (k = i; k < j; k++)
						{
							ranking[k] = ranking[k + 1];
						}
						ranking[j - 1] = best_insertion_alternative;
					}
					/*system("PAUSE");*/
				}
			}
		}
		/*for(k = 0; k < size; k++)
		{
		cout << ranking[k] << "\t";
		}*/
		/*cout << endl;*/
		/*system("PAUSE");*/
	}
	// Use Local Search to improve solution through insertion.

	// Adjust final ranking to take into account the P[i][j] values. If high, we would like i before j, even if there is no direct value
	// attached to this in the pricing problem. Obviously, the pricing objective value can not be lowered by making these adjustments.
	ranking = dominant_alternative(ranking, obj_coef, size, env, cplex_LP, c);

	output.resize(1);
	output[0].resize(size);
	for (j = 0; j < size; j++)
	{
		output[0][j].resize(size);
	}

	solution_value = 0;
	// Rewrite output to precedence matrix;
	for (i = 0; i < size; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			output[0][ranking[i]][ranking[j]] = 1;
			solution_value = solution_value + obj_coef[ranking[i]][ranking[j]];
		}
	}
	//cout << solution_value << endl;
	solution_value = solution_value - target_var;
	//cout << solution_value << endl;
	*ip_sol = solution_value;
	return output;
}

void IP_Update(IloEnv &env, IloCplex cplex_lp, ConstraintMatrix c_lp, IloRange sum, NumVarMatrix ip_vars, IloNumVar target_var, IloObjective ip_obj, int size)
{
	int i, j;

	ip_obj.setLinearCoef(target_var, -cplex_lp.getDual(sum));
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			ip_obj.setLinearCoef(ip_vars[i][j], cplex_lp.getDual(c_lp[i][j]));
		}
	}
}

IloModel IP_model(IloEnv &env, NumVarMatrix x, IloObjective obj, TwoDMatrix P, int size)
{
	int i, j, k;
	IloModel model(env);

	// Triangle inequality (enforce transitivity)
	for (i = 0; i<size; i++)
	{
		for (j = i + 1; j<size; j++)
		{
			for (k = j + 1; k<size; k++)
			{
				model.add(x[i][j] + x[j][k] + x[k][i] <= 2);
				model.add(x[j][i] + x[i][k] + x[k][j] <= 2);
			}
		}
	}

	// Enforce anti-symmetry
	for (i = 0; i<size; i++)
	{
		for (j = i; j<size; j++)
		{
			if (i == j)
				model.add(x[i][i] == 0);
			else
				model.add(x[i][j] + x[j][i] == 1);

		}
	}

	// No patterns that can't be used
	/*for(i=0; i<size; i++)
	{
	for(j=0 ; j<size; j++)
	{
	if(P[i][j] == 1 && i != j)
	{model.add(x[i][j] == 1);}
	}
	}*/

	// Add the objective function
	model.add(obj);

	return model;
}

NumVarMatrix IP_vars_cr(IloEnv &env, int size)
{
	int i;
	NumVarMatrix vars(env, size);
	for (i = 0; i < size; i++)
	{
		vars[i] = IloNumVarArray(env, size, 0, 1, ILOINT);
	}
	return vars;
}

IloModel LP_model(IloEnv &env, IloNumVarArray vars, IloObjective obj, IloRange sum, ConstraintMatrix c, int size)
{
	int i;
	IloModel model(env);

	model.add(obj);
	model.add(sum);
	for (i = 0; i < size; i++)
	{
		model.add(c[i]);
	}

	return model;
}

ConstraintMatrix LP_Cons(IloEnv &env, IloNumVarArray vars, ThreeIntMatrix start_col, TwoDMatrix P, int size)
{
	int i, j;
	double value;
	ConstraintMatrix c(env, size);

	for (i = 0; i < size; i++)
	{
		c[i] = IloRangeArray(env, size);
	}

	// The following has only been built for the basic starting pattersn ab...yz - zy...ba.
	for (i = 0; i < size; i++)
	{
		for (j = 0; j< size; j++)
		{
			if (i < j)
			{
				value = P[i][j] / (P[i][j] + P[j][i]);
				c[i][j] = IloRange(vars[0] + vars[1] >= value);
			}
			if (i > j)
			{
				value = P[i][j] / (P[i][j] + P[j][i]);
				c[i][j] = IloRange(vars[0] + vars[2] >= value);
			}
			if (i == j)
			{
				c[i][j] = IloRange(vars[0] >= 0);
			}
		}
	}

	return c;
}

IloRange LP_Sum(IloEnv &env, IloNumVarArray vars, ThreeIntMatrix start_col)
{
	int i;
	IloRange sum;

	sum = IloRange(-vars[1] >= -1);

	// vars 1 to ... are the variables corresponding to patterns. All of these must have -1 as coefficient.
	for (i = 1; i <= start_col.size(); i++)
	{
		sum.setLinearCoef(vars[i], -1);
	}

	return sum;
}

IloNumVarArray LP_vars_cr(IloEnv &env, ThreeIntMatrix start_col)
{
	int i;
	IloNumVarArray vars(env);

	// One variable for "z", then one per pattern created at start
	for (i = 0; i < start_col.size() + 1; i++)
	{
		vars.add(IloNumVar(env, 0, 1, ILOFLOAT));
	}

	return vars;
}

TwoDMatrix readdata(char inputfile[], int *size)
{
	int i, j;
	TwoDMatrix P_sub;
	ifstream inFile(inputfile);

	if (!inFile) // Check for existence of file, if not goto end.
	{
		//cout << endl << "Failed to open file ";
		goto end;
	}

	inFile >> *size;

	// Further definition of P
	P_sub.resize(*size);
	for (i = 0; i < *size; i++)
	{
		P_sub[i].resize(*size);
	}

	for (i = 0; i < *size; i++)
	{
		P_sub[i][i] = 0;
		for (j = 0; j < *size; j++)
		{
			inFile >> P_sub[i][j];
		}
	}
end:;
	return P_sub;
}

TwoIntMatrix readmissingdata(char inputfile[], int *size)
{
	int i, j;
	TwoIntMatrix P_sub;
	ifstream inFile(inputfile);

	if (!inFile) // Check for existence of file, if not goto end.
	{
	    cout << endl << "Failed to open missing data file ";
		goto end;
	}

	inFile >> *size;

	// Further definition of P
	P_sub.resize(*size);
	for (i = 0; i < *size; i++)
	{
		P_sub[i].resize(*size);
	}

	for (i = 0; i < *size; i++)
	{
		for (j = 0; j < *size; j++)
		{
			inFile >> P_sub[i][j];
		}
		P_sub[i][i] = 0;
	}
end:;
	return P_sub;
}

configuration readconfig(char configname[], int * error)
{
	configuration config;
	namespace po = boost::program_options;
	ifstream configfile(configname);
	if (!configfile)
	{
		//cout << "No configuration file found." << endl;
		*error = 1;
	}
	po::options_description file_options(configname);
	file_options.add_options()
		("Repetitions", po::value<int>(&config.repetitions)->required(), "")
		;
	po::variables_map vm;
	ifstream ifs(configname);
	if (!ifs)
	{
		//cout << "Failed to open Configuration File" << endl;
		*error = 1;
		return config;
	}
	else
	{
		try {
			store(po::parse_config_file(ifs, file_options), vm);
			notify(vm);
		}
		catch (exception e)
		{
			printf("Required values not included in config file"); 
			std::cout << "Press enter to continue" << endl; cin.get(); *error = 1;
		}
	}
	return config;
}

ThreeIntMatrix start(int size, TwoDMatrix P)
{
	int i, j;
	TwoDMatrix P_Help;
	ThreeIntMatrix start_col;

	//cout << "Choose Starting Heuristic Here; so far only basic is implemented, it will be chosen by default." << endl;
	P_Help = P;
	start_col.resize(2);
	for (i = 0; i < 2; i++)
	{
		start_col[i].resize(size);
		for (j = 0; j < size; j++)
		{
			start_col[i][j].resize(size);
		}
	}
	for (i = 0; i < size; i++)
	{
		start_col[0][i][i] = 0;
		start_col[1][i][i] = 0;
		for (j = i + 1; j < size; j++)
		{
			start_col[0][i][j] = 1;
			start_col[1][j][i] = 1;
		}
	}
	return start_col;
}

void outputpatterns(ThreeIntMatrix Patterns, int size)
{
	ofstream outputpatterns;
	outputpatterns.open("Outputpatterns.txt");
	int i, j, k;

	for (i = 0; i < Patterns.size(); i++)
	{
		outputpatterns << "Pattern " << i << endl;
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				if (k == 0)
					outputpatterns << "[" << Patterns[i][j][k] << ", ";
				else if (k == size - 1)
					outputpatterns << Patterns[i][j][k] << "]" << endl;
				else
					outputpatterns << Patterns[i][j][k] << ", ";
			}
		}
	}
}

IntVector dominant_alternative(IntVector ranking, TwoDMatrix obj_coef, int size, IloEnv &env, IloCplex cplex_LP, ConstraintMatrix c)
{
	int i, j, k;
	int best_insertion_alternative;
	bool improve;
	double solution_value;
	double insertion_value;

	// First, fix the preferences that contribute to the objective function by setting their coefficient to an arbitrary high number
	for (i = 0; i < size; i++)
	{
		for (j = i + 1; j < size; j++)
		{
			if (obj_coef[ranking[i]][ranking[j]] > 0)
			{
				obj_coef[ranking[i]][ranking[j]] = 100;
			}
		}
	}
	// Then change the coëfficients to give preference to dominant alternatives
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			if (obj_coef[i][j] != 100)
			{
				obj_coef[i][j] = c[i][j].getLB();
			}
		}
	}

	// Now run the improvement algorithm again.
	for (improve = 1; improve == 1;)
	{
		improve = 0;
		for (i = 0; i < size && improve == 0; i++)
		{
			// First compute the current value associated with the object in a specific location.
			solution_value = 0;
			for (j = 0; j < i; j++)
			{
				solution_value = solution_value + obj_coef[ranking[j]][ranking[i]];
			}
			for (j = i + 1; j < size; j++)
			{
				solution_value = solution_value + obj_coef[ranking[i]][ranking[j]];
			}
			// We have now established the value of the current position of alternative ranking[i].
			// Now calculate the current best position.
			// We test size+1 positions, starting at j=0, i.e. the alternative is placed first.
			// j = 1 means placed between position 0 and position 1
			// j = size + 1 means being placed in last position
			for (j = 0; j < size + 1 && improve == 0; j++)
			{
				insertion_value = 0;
				for (k = 0; k < j; k++)
				{
					insertion_value = insertion_value + obj_coef[ranking[k]][ranking[i]];
				}
				for (k = j; k < size; k++)
				{
					insertion_value = insertion_value + obj_coef[ranking[i]][ranking[k]];
				}
				// Rearange ranking based on new location if better position is found.
				if (insertion_value > solution_value + 0.0000001)
				{
					improve = 1;
					/*cout << "i = " << i << "\t j = " << j << endl;*/
					best_insertion_alternative = ranking[i];
					if (i > j) // The alternative is moved to the front.
					{
						for (k = i; k > j; k--)
						{
							ranking[k] = ranking[k - 1];
						}
						ranking[j] = best_insertion_alternative;
					}
					if (j > i)
					{
						for (k = i; k < j; k++)
						{
							ranking[k] = ranking[k + 1];
						}
						ranking[j - 1] = best_insertion_alternative;
					}
					/*system("PAUSE");*/
				}
			}
		}
	}
	return ranking;
}


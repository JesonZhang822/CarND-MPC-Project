#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

using namespace std;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// reference v ,epsi and cte
double ref_v = 100;
double ref_epsi = 0;
double ref_cte = 0;

// set the weight of the cost values
double cte_factor = 1800;//2000;//1500;
double epsi_factor = 1800;//2000;//1500;
double diff_delta_factor = 160;//200;//180;
double delta_factor = 10; //10
double a_factor = 10; //10
double diff_a_factor = 50;//20;

// set the starting index of x,y,psi,v,cte,epsi,delta and a in  vector vars
int x_start = 0;
int y_start = x_start + N;
int psi_start = y_start + N;
int v_start = psi_start + N;
int cte_start = v_start + N;
int epsi_start = cte_start + N;
int delta_start = epsi_start + N;
int a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    //Setp 0: Cost function
    //the cost based on the reference state : cross-track error,heading error and velocity error
    fg[0] = 0;
    for (int i = 0;i < N ;i++){

      fg[0] += cte_factor * CppAD::pow(vars[cte_start + i] - ref_cte,2);
      fg[0] += epsi_factor * CppAD::pow(vars[epsi_start + i] - ref_epsi,2);
      fg[0] += CppAD::pow(vars[v_start + i] - ref_v,2);

    }

    // Minimize the change rate
    for (int i =0;i < N -1 ;i++){

      fg[0] += delta_factor * CppAD::pow(vars[delta_start + i],2);
      fg[0] += a_factor * CppAD::pow(vars[a_start + i],2);
    }

    // Minimize the value gap between sequential actuations
    for (int i= 0;i < N -2 ;i++){
      fg[0] += diff_delta_factor * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i],2);
      fg[0] += diff_a_factor * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i],2);

    }
    
    //Step 1:Initialize the model to the initial state
    //fg[0] is reserved for the cost value
    //fg[x_start + 1] is reserved for the initial value of x
    //fg[x_start + 1 + i] is reserved for the ith of N values of x that the slover computes
    fg[x_start + 1] = vars[x_start];
    fg[y_start + 1] = vars[y_start];
    fg[psi_start + 1] = vars[psi_start];
    fg[v_start + 1] = vars[v_start];
    fg[cte_start + 1] = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];

    for (int i = 0; i < N-1;i++){

      //the state at time i
      AD<double> x0 = vars[x_start + i];
      AD<double> y0 = vars[y_start + i];
      AD<double> psi0 = vars[psi_start + i];
      AD<double> v0 = vars[v_start + i];
      AD<double> cte0 = vars[cte_start + i];
      AD<double> epsi0 = vars[epsi_start + i];

      // Only consider the actuation at time i
      AD<double> delta0 = vars[delta_start + i];
      AD<double> a0 = vars[a_start + i];

      //evaluate at polynomail at time i
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      // calculate the psides ,which is can be calculated as arctan(f'(x))
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0);

      // the state at time i+1
      AD<double> x1 = vars[x_start + i + 1];
      AD<double> y1 = vars[y_start + i + 1];
      AD<double> psi1 = vars[psi_start + i + 1];
      AD<double> v1 = vars[v_start + i + 1];
      AD<double> cte1 = vars[cte_start + i + 1];
      AD<double> epsi1 = vars[epsi_start + i + 1];  

      // constraints
      fg[x_start + i +2] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[y_start + i +2] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[psi_start + i +2] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
      fg[v_start + i + 2] = v1 - (v0 + a0 * dt);
      fg[cte_start + i +2] = cte1 - ((f0 - y0) + v0 * CppAD::sin(epsi0) * dt);
      fg[epsi_start + i +2] = epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);

    }

  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N - 1);
  // TODO: Set the number of constraints
  size_t n_constraints = 6 * N ;

  // get the value of x,y,psi,v,cte,epsi from state vector
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];


  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  //set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // set the range of delta values to [-25,25] in radiants
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  // for acceleration/deacceleration
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = - 1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // set the lower limits for the constraints
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  // set the upper limits for the constraints
  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  //std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;
  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  //return the predict x,y
  for (int i = 0; i < N-1; i++)
  {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }
  return result;
}

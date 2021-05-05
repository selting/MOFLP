#include "FLPModel.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>


FLPModel::FLPModel(std::string filename) {
	readFLP(filename);
	std::cout << *this << std::endl;


}

void FLPModel::readFLP(std::string filename)
{
	std::string line;
	std::ifstream in_file(filename.c_str());
	int cpt = 0;
	filename_ = filename;
	if (!in_file.is_open())
	{
		std::cerr << "can't open file " << filename << std::endl;
		exit(2);
	}
	while (!in_file.eof())
	{
		std::getline(in_file, line);
		std::stringstream tok(line);
		switch (cpt)
		{
		case 0:
			// first line: num customers and num facilities
			tok >> N_;
			tok >> H_;
			tok >> delta_;
			// allocate memory
			q_ = std::vector<int>(N_);
			c_ = std::vector<int>(H_);
			Q_ = std::vector<int>(H_);
			break;

		case 1:
			// customer demand
			for (int i = 0; i < N_; ++i)
			{
				tok >> q_[i];
			}
			break;

		case 2:
			// facility opening cost
			for (int i = 0; i < H_; ++i)
			{
				tok >> c_[i];
			}
			break;

		case 3:
			// facility capacity
			for (int i = 0; i < H_; ++i)
			{
				tok >> Q_[i];
			}
			break;

		default:
			// customer-to-facility distances
			std::vector<int> dvec;
			for (int i = 0; i < H_; ++i)
			{
				int d;
				tok >> d;
				dvec.push_back(d);
			}
			dist_.push_back(dvec);
		}
		cpt++;

	}
}



FLPModel::FLPModel(const FLPModel& original)
{
}

FLPModel::~FLPModel()
{
}

void FLPModel::makeCPLEXmodel()
{
	model_ = IloModel(env_);
	cplex_ = IloCplex(model_);

	// decision variables
	y_ = IloBoolVarArray(env_, H_);
	x_ = IloArray<IloBoolVarArray>(env_, N_);
	for (int i = 0; i < N_; ++i)
	{
		x_[i] = IloBoolVarArray(env_, H_);

	}

	// objective function f1
	IloExpr expr1(env_);
	for (int j = 0; j < H_; ++j)
	{
		expr1 += c_[j] * y_[j];
	}
	obj1_ = IloMinimize(env_, expr1, "min cost");

	// objective function f2
	IloExpr expr2(env_);
	for (int i = 0; i < N_; ++i) {
		for (int j = 0; j < H_; ++j)
		{
			expr2 += dist_[i][j] * x_[i][j];
		}
	}
	obj2_ = IloMinimize(env_, expr2, "min dist");


	// constraint: each customer is assigned to exactly one facility
	for (int i = 0; i < N_; ++i)
	{
		model_.add(IloSum(x_[i]) == 1);
	}

	// facility capacity constraints
	for (int j = 0; j < H_; ++j)
	{
		IloExpr expr(env_);
		for (int i = 0; i < N_; ++i)
		{
			expr += q_[i] * x_[i][j];
		}
		model_.add(expr <= Q_[j]);
	}

	// constraint: customer can only be assigned to facility that is open
	for (int i = 0; i < N_; ++i)
	{
		for (int j = 0; j < H_; ++j)
		{
			model_.add(x_[i][j] <= y_[j]);
		}

	}

	// turn off output
	cplex_.setOut(env_.getNullStream());

}

int FLPModel::solve()
{

	int opt = -1;
	try {
		cplex_.solve();
		if (cplex_.getStatus() == IloAlgorithm::Optimal) {
			opt = cplex_.getObjValue();
			std::cout << std::endl << "Objective=" << cplex_.getObjective().getName() << ", cost=" << get_f1() << ", dist=" << get_f2() << std::endl;
			for (unsigned int j = 0; j < H_; j++) {
				if (cplex_.getValue(y_[j]) > 0.5) {
					std::cout << "\tFacility " << j << ":";
					for (unsigned int i = 1; i < N_; i++) {
						if (cplex_.getValue(x_[i][j]) > 0.5) {
							std::cout << " " << i;
						}
					}
					std::cout << std::endl;
				}
			}
		}
		else {
			std::cout << "No feasible solution found" << std::endl;
		}
	}
	catch (IloException& e) {
		std::cerr << "Error solving model: " << e.getMessage() << std::endl;
		e.end();
	}
	return opt;
}


int FLPModel::pareto_boundary_v1()
{
	std::vector<int> opt_cost;
	std::vector<int> opt_dist;

	int omega = 1;
	int epsilon = 9999;

	// set epsilon constraint on f2
	IloRange ub_f2(env_, 0, obj2_.getExpr(), epsilon);
	model_.add(ub_f2);


	while (true)
	{
		// set f1 as objective
		model_.add(obj1_);
		
		//optimize the first objective
		cplex_.exportModel("FLP.lp");
		int z1 = solve();
		if (cplex_.getStatus() != IloAlgorithm::Optimal)
		{
			break;
		}

		//collect the objective values
		//int f1 = get_f1();
		//int f2 = get_f2();
		//opt_cost.push_back(f1);
		//opt_dist.push_back(f2);

		//set bound on f1
		IloRange f1_bound(env_, 0, obj1_.getExpr(), z1);
		model_.add(f1_bound);

		//change the objective
		model_.remove(obj1_);
		model_.add(obj2_);

		//optimize for the second objective
		int z2 = solve();
		if (cplex_.getStatus() != IloAlgorithm::Optimal)
		{
			break;
		}

		//collect the objective values
		int f1 = get_f1();
		int f2 = get_f2();
		opt_cost.push_back(f1);
		opt_dist.push_back(f2);

		model_.remove(f1_bound);		
		ub_f2.setUB(f2 - omega);

		model_.remove(obj2_);
	}

	// print the pareto frontier
	std::cout << "f1: cost" << "\tf2: dist" << std::endl;
	for (int i = 0; i < opt_cost.size(); ++i)
	{
		std::cout << opt_cost[i] << "\t" << opt_dist[i] << std::endl;
	}
}


int FLPModel::pareto_boundary_v2()
{
	std::vector<int> opt_cost;
	std::vector<int> opt_dist;

	int omega = 1;
	int epsilon = 9999;

	// set epsilon constraint on f2
	IloRange ub_f2(env_, 0, obj2_.getExpr(), epsilon);
	model_.add(ub_f2);
	
	// set f1 as objective
	model_.add(obj1_);

	
	while (true)
	{
		//optimize the first objective
		int z1 = solve();
		if (cplex_.getStatus() != IloAlgorithm::Optimal)
		{
			break;
		}

		//collect the objective values
		int f1 = get_f1();
		int f2 = get_f2();
		opt_cost.push_back(f1);
		opt_dist.push_back(f2);

		//update upper bound constraint on f2
		ub_f2.setUB(f2 - omega);
	}

	// print the pareto frontier
	std::cout << "f1: cost" << "\tf2: dist" << std::endl;
	for (int i = 0; i < opt_cost.size(); ++i)
	{
		std::cout << opt_cost[i] << "\t" << opt_dist[i] << std::endl;
	}
}

int FLPModel::get_f1()
{
	int f1 = 0;
	for (int i = 0; i < H_; ++i)
	{
		f1 += cplex_.getValue(y_[i]) * c_[i];
	}
	return f1;
}

int FLPModel::get_f2()
{
	int f2 = 0;
	for (int i = 0; i < N_; ++i)
	{
		IloNumArray x(env_);
		cplex_.getValues(x, x_[i]);
		for (int j = 0; j < H_; ++j)
		{
			f2 += x[j] * dist_[i][j];
		}
	}
	return f2;
}

std::ostream& operator<<(std::ostream& os, const FLPModel& instance)
{
	os << "Facility Location instance loaded from " << instance.filename()
		<< " with " << instance.H() << " facilities,  " << instance.N() << " Customers, "
		<< "delta = " << instance.delta() << std::endl;
	return os;
}

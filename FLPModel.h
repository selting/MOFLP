#pragma once

#include <iostream>
#include <vector>
#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>



class FLPModel
{
protected:
	// instance name/file location
	std::string filename_;
	// number of customers
	int N_;
	// number of facilities
	int H_;
	// customer demand
	std::vector<int> q_;
	// cost for opening a facility
	std::vector<int> c_;
	// capacity of facilities
	std::vector<int> Q_;
	// distance from customer i to facility j
	std::vector<std::vector<int>> dist_;
	// delta parameter for epsilon-constraints
	int delta_;


	// CPLEX model members
	IloEnv env_;
	IloCplex cplex_;
	IloModel model_;
	// decision var y: is facility j opened?
	IloBoolVarArray y_;
	// decision var y: is customer i assigned to facility j=
	IloArray<IloBoolVarArray> x_;
	// objective 1
	IloObjective obj1_;
	// objective 2
	IloObjective obj2_;

	
	// read an instance from the OPL format
	void readFLP(std::string filename);

	int get_f1();
	int get_f2();



public:
	// constructor that accepts filename
	FLPModel(std::string filename);

	// copy constructor
	FLPModel(const FLPModel& original);
	//destructor
	~FLPModel();

	// getters
	int N() const { return N_; };
	int H() const { return H_; };
	int q(int i) const { return q_[i]; };
	int c(int i) const { return c_[i]; };
	int Q(int i) const { return Q_[i]; };
	int dist(int i, int j) const { return dist_[i][j]; };
	int delta() const { return delta_; };
	std::string filename() const { return filename_; };

	// create a cplex model from the instance
	void makeCPLEXmodel();
	
	// solve model using cplex (uni-objective?)
	int solve();
	
	int pareto_boundary();
};

std::ostream& operator<<(std::ostream& os, const FLPModel& instance);
// An Example for Benders Decomposition
/************************************************************************************
 * This example on Benders Decomposition is taken from the tutorial on Benders Decomposition.
 * Refer to Benders-Slides.pdf
 * Minimize -2X1 + 3X2 - 4X3 - 5Y1 + 2Y2 - 9Y3 
 * s.t.:
 * -2X1 - 3X2 -6X3 - 5Y1 + 3Y2 - 7Y3 >= 2
 * -3X1 + 1X2 -3X3 - 4Y1 - 2Y2 - 4Y3 >= -10
 * X1, X2, X3 >= 0, Y integer in {0,5}
 ***********************************************************************************/

 //////USING IBD//////

#pragma warning(disable : 4996) //For Visual Studio 2012
#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <iosfwd>
#include <string>
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector> //for vectors
#include <math.h>

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

using namespace std;

ILOSTLBEGIN

int main(int argc, char** argv)
{
	IloEnv env;
	try
	{
		/// DECISION VARIABLES///
		IloNumVarArray X(env, 3, 0, IloInfinity, ILOFLOAT);
		// cout << X << endl;
		IloNumVarArray Y(env, 3, 0, 5, ILOINT);


		//////////DEVELOP GENERIC MODEL //////////////////////////

		//////SET MASTER PROBLEM///////////////////////////
		IloModel model_master(env);
		IloExpr Objective_master(env);
		Objective_master = -2 * X[0] + 3 * X[1] - 4 * X[2] - 5 * Y[0] + 2 * Y[1] - 9 * Y[2];

		// Objective Function//
		model_master.add(IloMinimize(env, Objective_master));

		//Constraints//
		IloExpr constraint_1(env);
		constraint_1 = -2 * X[0] - 3 * X[1] - 6 * X[2] - 5 * Y[0] + 3 * Y[1] - 7 * Y[2];
		model_master.add(constraint_1 >= 2);


		IloExpr constraint_2(env);
		constraint_2 = -3 * X[0] + 1 * X[1] - 3 * X[2] - 4 * Y[0] - 2 * Y[1] - 4 * Y[2];
		model_master.add(constraint_2 >= -10);

		IloCplex cplex_master(model_master);
		Objective_master.end();
		constraint_1.end();
		constraint_2.end();
		cplex_master.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		cplex_master.setWarning(env.getNullStream()); //This is to supress warning messages on screen


		cplex_master.setParam(IloCplex::Param::Benders::Strategy, IloCplex::BendersFull);//For Cplex Internal Benders

		///SOLVING///
		cout << "SOLVING MASTER PROBLEM" << endl;

		if (!cplex_master.solve())
		{
			cout << "Failed" << endl;
			throw(-1);
		}

		cout << "Master Problem Solution Status: " << cplex_master.getCplexStatus() << endl;

		double obj_value = cplex_master.getObjValue();
		cout << "Objective Value: " << obj_value << endl;

		double Y1_value = cplex_master.getValue(Y[0]);
		double Y2_value = cplex_master.getValue(Y[1]);
		double Y3_value = cplex_master.getValue(Y[2]);

		double X1_value = cplex_master.getValue(X[0]);
		double X2_value = cplex_master.getValue(X[1]);
		double X3_value = cplex_master.getValue(X[2]);

		cout << "Y1 : " << Y1_value << endl;
		cout << "Y2 : " << Y2_value << endl;
		cout << "Y3 : " << Y3_value << endl;

		cout << "X1 : " << X1_value << endl;
		cout << "X2 : " << X2_value << endl;
		cout << "X3 : " << X3_value << endl;



		model_master.end();

		cplex_master.end();
	}//try
	catch (IloException& e)
	{
		env.out() << "ERROR: " << e << endl;
	}
	catch (...)
	{
		env.out() << "Unknown exception" << endl;
	}
	env.end();
	return 0;
}
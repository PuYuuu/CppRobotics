#include <iostream>

#include <Eigen/Dense>
#include <OsqpEigen/OsqpEigen.h>

int main(int argc, char** argv)
{
    // Problem description
    // min  0.5 * X^T * H * X + Q * X
    // s.t. L <= A * X <= U
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;

    hessian.resize(3,3);
    hessian.insert(0,0) = 1;
    hessian.insert(1,0) = -1;
    hessian.insert(2,0) = 1;
    hessian.insert(0,1) = -1;
    hessian.insert(1,1) = 2;
    hessian.insert(2,1) = -2;
    hessian.insert(0,2) = 1;
    hessian.insert(1,2) = -2;
    hessian.insert(2,2) = 4;
    std::cout << "hessian:" << std::endl << hessian << std::endl;

    gradient.resize(3);
    gradient << 2, -3, 1;
    std::cout << "gradient:" << std::endl << gradient << std::endl;

    linearMatrix.resize(4,3);
    linearMatrix.insert(0,0) = 1;
    linearMatrix.insert(1,0) = 0;
    linearMatrix.insert(2,0) = 0;
    linearMatrix.insert(3,0) = 1;
    linearMatrix.insert(0,1) = 0;
    linearMatrix.insert(1,1) = 1;
    linearMatrix.insert(2,1) = 0;
    linearMatrix.insert(3,1) = 1;
    linearMatrix.insert(0,2) = 0;
    linearMatrix.insert(1,2) = 0;
    linearMatrix.insert(2,2) = 1;
    linearMatrix.insert(3,2) = 1;
    std::cout << "linearMatrix:" << std::endl << linearMatrix << std::endl;

    lowerBound.resize(4);
    lowerBound << 0, 0, 0, 0.5;
    std::cout << "lowerBound:" << std::endl << lowerBound << std::endl;

    upperBound.resize(4);
    upperBound << 1, 1, 1, 0.5;
    std::cout << "upperBound:" << std::endl << upperBound << std::endl;

    int NumberOfVariables = 3;
    int NumberOfConstraints = 4;

    // instantiate the solver
    OsqpEigen::Solver solver;

    // settings
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(NumberOfVariables);
    solver.data()->setNumberOfConstraints(NumberOfConstraints);
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound);
    solver.data()->setUpperBound(upperBound);

    if (!solver.initSolver()) {
        return 1;
    }

    if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
        return 1;
    }

    Eigen::VectorXd QPSolution = solver.getSolution();
    std::cout << "QPSolution:" << std::endl << QPSolution << std::endl;

    return 0;
}

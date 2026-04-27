# Unconstrained Multivariable Optimization

This repository contains a MATLAB implementation and analysis of classical **multivariable optimization methods** for unconstrained function minimization.

The project focuses on both the **algorithmic development** and the **comparative evaluation** of first- and second-order methods, highlighting how different step size strategies affect convergence behavior.

It represents represents the second assignment of an *Optimization Techniques* course and serves as a continuation of a previous study on one-dimensional optimization.
See also: [One-Dimensional Optimization](https://github.com/thodorisvlachakis/one-dimensional-optimization)

---

## 📌 Overview

This repository provides a structured and comprehensive study of multivariable optimization techniques with experimental evaluation through numerical simulations and visualizations, highlighting the importance of step size selection and the trade-offs between convergence speed and computational cost.

The objective of this project is to investigate the performance of optimization algorithms when applied to a nonlinear two-variable function.

The project focuses on: 

* Comparison between **Steepest Descent**, **Newton**, and **Levenberg–Marquardt** methods
* Investigation of **step size selection strategies**
* Analysis of **convergence speed, stability, and accuracy** for each method
* Evaluation under **different initial conditions**

All methods are tested on the same function and starting points to ensure a fair and structured comparison.

---

## 🧠 Problem Description

We consider the minimization of the function:

f(x, y) = x³ · exp(-x² - y⁴)

The exact minimizer is computed analytically by solving:

∇f(x, y) = 0

and is used as a reference point for evaluating the numerical methods.

---

## ⚙️ Methods Implemented

### 🔹 Steepest Descent Method

A first-order method using the negative gradient direction.

### 🔹 Newton’s Method

A second-order method utilizing the Hessian matrix for faster convergence near the solution.

### 🔹 Levenberg–Marquardt Method

A hybrid approach combining gradient descent and Newton’s method, ensuring positive definiteness of the update matrix.

---

## 📏 Step Size Strategies

Each method is implemented with three different step size selection techniques:

### 🔸 Constant Step Size

A fixed step size for all iterations.

* Simple implementation
* Requires careful tuning
* May lead to slow convergence or instability

---

### 🔸 Armijo Rule (Backtracking Line Search)

An adaptive strategy that ensures sufficient decrease of the objective function.

* Dynamically adjusts the step size
* Improves robustness
* Widely used in practice

---

### 🔸 Inner Optimization (Exact Line Search)

The step size is computed by solving an auxiliary optimization problem:

min f(xₖ + g·dₖ)
over g ≥ 0 (condition for descent step)

This is implemented numerically using the **Golden Section Method**.

* Provides more accurate step selection
* Although, it has higher computational cost per iteration

---

## 🧪 Experimental Setup

The algorithms are tested using:

* **Three different starting points**:

  * (0, 0)
  * (-1, -1)
  * (1, 1)

* **Termination criterion**:
  ‖∇f(x, y)‖ < ε, with ε = 10⁻³

For each configuration, the following are recorded:

* Number of iterations
* Final approximation of the minimizer
* Function value progression
* Convergence plots

---

## 🎯 Key Observations

* Newton-based methods converge significantly faster when close to the solution.
* Newton’s method, despite its theoretically fast convergence, may fail to converge to the correct minimizer if the Hessian is not positive definite.
In this problem, the Hessian of the objective function is not positive definite even from the first iteration for all tested initial points, leading to non-convergent or unstable behavior. Find more details in the report analysis.
* Steepest Descent is more stable but typically slower.
* Levenberg–Marquardt provides a good balance between robustness and speed.
* Armijo rule improves stability compared to constant step size.
* Exact line search (inner optimization) yields accurate steps but increases computational cost.

---

## 📁 Project Structure

```
unconstrained-multivariable-optimization
│
├── src/
│   ├── SecondLaboratoryExerciseCode.m              # Main script
│   │
│   ├── methods/                                    # Optimization algorithms
│   │   ├── steepest_descent/
│   │   │   ├── SteepestDescentMethodWithConstantDescentStep.m
│   │   │   ├── SteepestDescentMethodDescentStepByInnerOptimization.m
│   │   │   └── SteepestDescentMethodArmijoStepSizeRule.m
│   │   │
│   │   ├── newton/
│   │   │   ├── NewtonMethodWithConstantDescentStep.m
│   │   │   ├── NewtonMethodDescentStepByInnerOptimization.m
│   │   │   ├── NewtonMethodArmijoStepSizeRule.m
│   │   │   ├── NewtonMethodWithConstantDescentStepAmendment.m
│   │   │   ├── NewtonMethodDescentStepByInnerOptimizationAmendment.m
│   │   │   └── NewtonMethodArmijoStepSizeRuleAmendment.m
│   │   │
│   │   └── levenberg_marquardt/
│   │       ├── LevenbergMarquadtMethodWithConstantDescentStep.m
│   │       ├── LevenbergMarquadtMethodDescentStepByInnerOptimization.m
│   │       └── LevenbergMarquadtMethodArmijoStepSizeRule.m
│   │
│   ├── step_size_strategies/                       # Step size selection strategies
│   │   ├── ArmijoStepSizeRule.m
│   │   ├── DescentStepByInnerOptimization.m
│   │   └── GoldenSectionMethod.m
│   │
│   └── utils/                                      # Helper functions
│       └── HessianOfFunctionAtSpecificPoint.m
│
├── docs/                                           # Statement & Report
│   ├── lab02.pdf
│   └── report_lab02.pdf
│
├── README.md
└── .gitignore
```

---

## ▶️ How to Run

### 🔧 Requirements

* MATLAB (any recent version)

### 🚀 Execution

1. Open MATLAB
2. Navigate to the `src/` directory
3. Run the main script:

```matlab
SecondLaboratoryExerciseCode.m
```

The script will:

* Execute all implemented methods
* Generate convergence plots
* Display useful results in the console, used in analysis in the report

---

## 📝 Notes

* The project emphasizes **comparative understanding**, not just implementation.
* Figures are generated dynamically by running the code.
* A detailed analysis, with plots, interpretations, comparative analysis and final conclusions regarding the convergence and efficiency of the algorithms studied, is provided in report_lab02.pdf file
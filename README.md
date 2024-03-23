# Problem Set 4 (PS4): Refactor the _Klebsiella oxytoca_ Simulation
This problem set will familiarize students with [refactoring code](https://en.wikipedia.org/wiki/Code_refactoring) and solve the Kompala Cybernetic Model for the batch growth of _Klebsiella oxytoca_ on sugar mixtures using the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) package. The model is described [here](paper/Kompala-BiotechBioengineering-1986.pdf).

## Prerequisites
Some `secret people` that only I can see and communicate with have told me they like our `solve` method signature __much better__ than [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). However, they have acknowledged that the [DifferentialEquations.jl package](https://diffeq.sciml.ai/stable/) is a powerful and flexible package for solving differential equations in Julia, i.e., it has better performance and more features than the solvers we have implemented in class.

They have also told me that they are interested in using the [DifferentialEquations.jl package](https://diffeq.sciml.ai/stable/) to solve the Kompala Cybernetic Model for the batch growth of _Klebsiella oxytoca_ on sugar mixtures. They have asked me to refactor the codes we have developed in the examples and labs to use the [DifferentialEquations.jl package](https://diffeq.sciml.ai/stable/) but to keep our `solve` method signature the same.

## Tasks
1. Complete the implementation of the `_mysolve(...)` method in the [Solver.jl file](src/Solver.jl) in the `src` directory. 
2. Complete the implementation of the `runme_task_1.jl` script. This script will call the `mysolve(...)` method to solve the model equations using the `RK4` solver from the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) package. 
    * We will simulate `Fig. 6` of the paper. Let glucose be state index `1` and xylose state `2`. The pseudo enzyme for glucose $e_{1}$ will be index `3`. The pseudo enzyme for xylose $e_{2}$ will be index `4`, and the biomass $C$ will be state `5`. This is the same case as `Lab-9d`.
    * The `runme_task_1.jl` script will generate a time array `T` and a state array `X` holding the model solutions. This data is saved to a binary `jld2` file using the [JLD2.jl package](https://github.com/JuliaIO/JLD2.jl) in the `mydata` directory.
3. Execute the `testme_task_1.jl` script. This script will load the model solutions from the `jld2` file and compare them to the expected solutions, which are provided in the `testdata` directory. We will compare the sizes of your time array `T` and the state array `X` and the values of `T` and `X` to the expected solutions. The tests will pass if the sizes and values match the expected solutions.

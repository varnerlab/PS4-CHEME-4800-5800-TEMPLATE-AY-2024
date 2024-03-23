function _mysolve(problem::MySimpleProblemModel, solver::MyForwardEulerMethod)::Tuple{Array{Float64,1}, Array{Float64,2}}
    
    # initialize -
    t0, tf, dt = problem.time_span;
    initial_conditions = problem.initial_conditions;
    time_array = range(t0, tf, step=dt) |> collect;
    X = Array{Float64,2}(undef, length(time_array), length(initial_conditions));

    # main loop -
    for i ∈ eachindex(time_array)
       
        if (i == 1)
            for j ∈ eachindex(initial_conditions)
                X[i,j] = initial_conditions[j];
            end
        else
            X[i,:] = X[i-1,:] + dt*problem.model(X[i-1,:], i-1, problem.parameters);
        end
    end

    # return the (T,X) tuple -
    return (time_array, X);
end

function _mysolve(problem::MySimpleProblemModel, solver::MyRungeKuttaMethod)::Tuple{Array{Float64,1}, Array{Float64,2}}
    
    # initialize -
    t0, t1, dt = problem.time_span; # t0 = initial time, t1 = final time, dt = time step
    ic = problem.initial_conditions; # initial conditions

    # TODO: get the parameters from the problem object, construct an ODEProblem object, and solve it using the Runge-Kutta method.
    # TODO: use the saveat keyword argument on ODEProblem to save at time steps of dt.
    # TODO: set the reltol and abstol keyword arguments to 1e-8 on the solve function.
    # TODO: the solve function should use the RK4() method, and store the solution in the variable soln.
    throw(ArgumentError("Not implemented yet for Runge-Kutta method."));

   
    # --- DO NOT MODIFY BELOW THIS LINE ------------------------------------------------------ %
    T = soln.t
    tmp = soln.u
    number_of_time_steps = length(T)
    number_of_dynamic_states = length(ic);
    X = Array{Float64,2}(undef, number_of_time_steps,  number_of_dynamic_states);
    for i ∈ 1:number_of_time_steps
        soln_vector = tmp[i]
        for j ∈ 1:number_of_dynamic_states
            X[i,j] = soln_vector[j]
        end
    end
    # --- DO NOT MODIFY ABOVE THIS LINE ------------------------------------------------------ %

    # return - 
    return (T, X)
end

"""
    solve(balances::Function, tspans::Tuple{Float64,Float64,Float64}, initial_conditions::Array{Float64,1}, parameters::Dict{String,Array{Float64,1}}; 
        solver::AbstractIVPSolverType = MyForwardEulerMethod())::Tuple{Array{Float64,1}, Array{Float64,2}}

Solve the system of ODEs defined by the `balances` function, with the given `initial_conditions` and `parameters` over the time span `tspans`. 
The `solver` keyword argument is used to specify the method to solve the system of ODEs. The default solver is the `MyForwardEulerMethod`.

### Arguments
- `balances::Function`: the function that defines the system of ODEs.
- `tspans::Tuple{Float64,Float64,Float64}`: the time span of the simulation.
- `initial_conditions::Array{Float64,1}`: the initial conditions of the system of ODEs.
- `parameters::Dict{String,Array{Float64,1}}`: the parameters of the system of ODEs.

### Optional keyword arguments
- `solver::AbstractIVPSolverType = MyForwardEulerMethod()`: the method to solve the system of ODEs.

### Returns
- `Tuple{Array{Float64,1}, Array{Float64,2}}`: a tuple containing the time array and the solution array.
"""
function mysolve(balances::Function, tspan::Tuple{Float64,Float64,Float64}, initial::Array{Float64,1}, parameters::Dict{String, Any}; 
    solver::AbstractIVPSolverType = MyForwardEulerMethod())::Tuple{Array{Float64,1}, Array{Float64,2}}
    
    # create the problem, object -
    problem = build(MySimpleProblemModel, (
        parameters = parameters,
        initial_conditions = initial,
        time_span = tspan,
        model = balances
    ));

    # solve the problem using the appropriate solver -
    return _mysolve(problem, solver)
end
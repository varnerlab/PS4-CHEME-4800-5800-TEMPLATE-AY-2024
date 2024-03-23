# include -
include("Include.jl");

# TODO: Specify the model and simulation parameters
# ...
throw(ErrorException("The runme_task_1.jl script is not implemented yet!"));


# --- DO NOT MODIFY BELOW THIS LINE ---------------------------------------------------------- #
# call our solver function
(T,X) = mysolve(balances, (0.0, 10.0, 0.1), xₒ, parameters, solver = MyRungeKuttaMethod());

# save the results to a file -
save(joinpath(_PATH_TO_MY_DATA, "mysimulationdata.jld2"), Dict("T" => T, "X" => X));
# --- DO NOT MODIFY ABOVE THIS LINE ---------------------------------------------------------- #
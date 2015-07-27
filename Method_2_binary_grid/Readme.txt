This folder contains the code for Method 2: value iteration based on discretized states. 

Discretize0.m and Discretize1.m returns the discretized grids for a continuous states given certain problem, where Discretize0.m follows the zeroth-order discretization and Discretize1.m follows the first-order discretization.

zero_order_grid.m and first_order_grid.m returns the value function after performing value iteration on discretized states, where the grids are discretized by zeroth-order or first-order discretization. 

Sequential.m returns the lower and upper bound for each time step of the threshold-based decision policy for a sequential problem specified by the parameters. There is a test case included as well. 

policy_eval.m is doing the policy evaluation using value iteration, where given a threshold, the infinite horizon value function is plotted.
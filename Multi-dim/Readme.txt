This folder contains the codes for multi-dimensional discretization. 

npermutek.m and mydistance.m are basic function that do specific computation. 

multidimension.m returns the discretized grids given any dimension. 

discretization3.m shows the grids for three dimension with 20 grids on each axis.

Discretize0.m, Discretize1.m and Discretize_infinite.m returns the value function and index by doing value iteration either with finite horizon or infinite horizon, and the grids are either discretized using zeroth-order or first-order discretization.

first_order.m and zero_order.m shows the optimal threshold for a specific three dimensional problem. 

test.m checks whether Discretize0.m and Discretize1.m are working correctly. 

This folder contains the code for orthogonal search, the idea is to iteratively find the optimal threshold for one decision maker by fixing the threshold of the other decision maker. 

iterations.m is the main function, it calls find_bound.m to find the optimal threshold for one decision maker using value iteration and calls absorption.m to find the conditionally probability for one decision maker using the absorption probability of Markov chain. 

iterations_norm.m is the modified iterations.m to solve the problem that has normally distributed probability matrix.
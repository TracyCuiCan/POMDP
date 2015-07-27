This folder contains the codes for direct search, where we uses non-convex optimization  function fminsearchbnd.m to directly search for optimal threshold and computes the conditional probability using absorption probability of Markov chain.  

fminsearchbnd.m is the Matlab function to use.

fmin.m and fmin_des.m are the functions where one set parameters and call fminsearchbnd.m to find the optimal thresholds. fmin.m solves centralized problem while fmin_des.m solves decentralized problem.

The functions that fminsearchbnd.m evaluate on are absorption1.m and absorption2.m, which uses absorption probability to find the conditional probability for one decision maker and two decision makers.

fmin_des_norm.m is the modified fmin_des.m for problem that has normally distributed probability.
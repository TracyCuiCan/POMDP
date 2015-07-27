This folder contains the code for Method 3: non-convex optimization combined with Monte Carlo sampling and asymptotic expression. 

fminsearchbnd.m is the Matlab function to use.

fmin.m and fmin_des.m are the functions where one set parameters and call fminsearchbnd.m to find the optimal thresholds. fmin.m solves centralized problem while fmin_des.m solves decentralized problem.

The functions that fminsearchbnd.m evaluate on are asymp_centrl.m, asymp_decentrl.m, MC_centrl.m and MC_decentrl.m, which uses asymptotic expression or Monte Carlo sampling and solve centralized or decentralized problem. 

MonteCarlo2.m is the Monte Carlo sampling code for two decision makers, it returns the number of samplings and the type1 and type2 error.
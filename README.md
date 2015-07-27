# POMDP
This work contains most of my Master's research project. It's about finding optimal thresholds for sequential hypothesis testing.

Sequential hypothesis testing is a partially observable Markov decision problem. In a sequential test, there are two kinds of errors. 
We may reject the null hypothesis when it is true (also called missed detection), or we may accept the null hypothesis when some 
alternative hypothesis is true (also called false alarm). There are cost incurred when we make a wrong decision and cost incurred 
when we make an additional observation. The goal is to design an optimal stopping rule that minimize the total cost. 

I used different ways to compute the optimal thresholds, including Sondik's enumeration method (see Method_1_alpha_vector); 
value iteration based on discretizing continuous belief state (see Method_2_binary_grid); non-convex optimization combined with 
Monte-Carlo sampling and asymptotic expression (see Method_3_MC&asymp); non-convex optimization combined with computing absorption 
probability of Markov chain (see Method_4_OS&DS). Discretization of multi-dimensional belief state is also included (see Multi-dim).

All codes are written in Matlab m. files. I hope this helps anyone who is interested in similar research in this field.

%%% EBS-SCA --An Enhanced Brain Storm Sine Cosine Algorithm for Global Optimization Problems

% This software is copyrighted by CHUNQUAN LI, ZU LUO, ZHENSHOU SONG, FENG YANG, JINGHUI FAN, and PETER X. LIU.  
%
% Permission is granted to copy and use the software for scientific, 
% noncommercial purposes, provided this copyright notice is retained and the 
% origin of the code is cited. The software is provided "as is" and without 
% any warranties, express or implied.

This package contains three parts:

1.CEC 2013 benchmark functions: 

	cec13_func.cpp: The function file for the minimization problems.

	input_data: The input data files of the 28 benchmark functions.

	CEC 2013 benchmark functions.pdf: The introduction document of CEC 2013 benchmark functions.

2. EBS_SCA.m : The main file for EBS-SCA algorithm.



3. TEST.m is an example test code with EBS-SCA algorithm.


	% max_iteration = 60000              Maximum number of iterations
	% max_FES = 300000                   Maximum number of fitness evaluations
	% N = 50                             Population size
	% D = 30                             Dimension
	% M = 10                             Number of clusters
	% lb = -100                          Lower bound of a problem
	% ub = 100                           Upper bound of a problem

Note, operating environment configuration requires:

	MatlabR2015a or a newer version;
	
	The computer is configured with a C++ compiler.
	


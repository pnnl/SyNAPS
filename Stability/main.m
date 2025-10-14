%% Testing parametric set identification for stability

clc; clear all; 
% close all;

addpath('./functions/')

% load example
model_opt.sysidx = 10;
model_opt.MOR   = 'yes';        % OPTIONS: {'yes','no'}
model_opt.xfrm  = 'no';        % OPTIONS: {'yes','no'}

model_opt.T = [1 2; -0.5 1];
% model_opt.T = [1 -1 0; 1 1 0; 0 0 1];%diag([1 1 1]); %[1 0 1; 0 1 0; -1 0 1];

[q,q0,xfrm_q,xfrm_q0,Aq,Ai] = extract_model_info(model_opt);

% identify stable region
solver_opt.mode = 'symmetric';      % OPTIONS: {'symmetric','positive','negative'}
solver_opt.norm = 'fro';            % OPTIONS: {'fro',Inf,2,1}
solver_opt.algo = 'oneshot';        % OPTIONS: {'oneshot','combine'}
solver_opt.iter = 'yes';            % OPTIONS: {'yes','no'}
solver_opt.test = 'yes';            % OPTIONS: {'yes','no'}
[optsol,constr_func] = solve_optimization(xfrm_q,xfrm_q0,Ai,solver_opt); 

% now plot and verify
plot_opt.figtyp = 'verify';         % OPTIONS: {'verify','demo'}           
plot_opt.valtyp = 'rel';            % OPTIONS: {'abs','rel'}
range_of_intrst = 1e0*[-1.3 1.3 -0.1 0.1 -3 3];
% range_of_intrst = [0 20 0 10];
plot_opt.npoint = 5e1;
plot_opt.algo   = solver_opt.algo;
plot_opt.tol    = 1e-6;
for ivar = 1:length(q)-1
    for jvar = ivar+1:length(q)
        fprintf('\n*** Plotting variables: %s vs. %s ***\n',char(q(ivar)),char(q(jvar)))
        plot_opt.varidx = [ivar jvar]; 
        plot_opt.ranges = range_of_intrst([2*ivar-1 2*ivar 2*jvar-1 2*jvar]);
        [sample_set, assess_stability] = plot_results(q,q0,Aq,constr_func,plot_opt);
    end
end
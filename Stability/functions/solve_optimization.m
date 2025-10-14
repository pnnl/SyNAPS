%% optimization formulation and solution

function [optsol,constr_func] = solve_optimization(q,q0,Ai,sol_opt)

nvar = size(Ai{1},1);  % number of state variables
pvar = length(q);   % number of design parameters

% normalize for better performance
A_norm = norm(Ai{1});
for iA = 1:length(Ai)
    Ai{iA} = Ai{iA}/A_norm;
end

if strcmp(sol_opt.iter,'yes')
    % optimize the variable transformation
    cmap = optimal_var_transform(Ai(2:end),eye(nvar),sol_opt.norm);    
elseif strcmp(sol_opt.iter,'no')
    cmap = eye(pvar);
else
    error('solver_opt.iter must be either ''yes'' or ''no''!')
end

% define flags to include specific constraints
incl_constr.pos = 1;
incl_constr.neg = 1;
if contains(sol_opt.mode,'pos')
    incl_constr.neg = 0;
elseif contains(sol_opt.mode,'neg')
    incl_constr.pos = 0;
elseif ~contains(sol_opt.mode,'sym')
    error('unknown optmization mode selected!')
end

switch sol_opt.algo
    case 'oneshot'
        % define objective function weighing factors
        Wobj = ones(pvar,1);
        
        if strcmp(sol_opt.iter,'yes')
            old_c = zeros(pvar);
            niter = 1;
            fprintf('\n*** STARTING iterative transformation loop ***\n')
            while norm(full(old_c-cmap)) > 1e-3
                old_c = cmap;

                % solve optimization
                [P,epsvar] = optimization_problem(Ai,old_c,sol_opt,Wobj,incl_constr);

                % optimize the variable transformation
                cmap = optimal_var_transform(Ai(2:end),P,sol_opt.norm);
                
                fprintf('\nITER %2d: %f',niter,norm(full(old_c-cmap)))                
                niter = niter + 1;
            end
            fprintf('\n*** COMPLETED iterative transformation loop ***\n')
            cmap
        elseif strcmp(sol_opt.iter,'no')
            % solve optimization
            [P,epsvar] = optimization_problem(Ai,cmap,sol_opt,Wobj,incl_constr);
        else
            error('solver_opt.iter must be either ''yes'' or ''no''!')
        end

        % rescale epsvar
        epsvar = rescale_epsvar(P,epsvar,Ai{1});

        qvec = inv(cmap) * transpose(q-q0);

        optsol.P{1}        = P;
        optsol.cmap{1}     = cmap;
        optsol.epsvar{1}   = epsvar;

        % finally: retrieve optimal solution, f_constr>=0
        constr_func = sym(zeros(nvar+2*pvar,1));
        
        constr_func(1:nvar) = vpa(simplify(1-optsol.epsvar{1}*abs(qvec)));
        
        big_M = max(1,norm(q0))*1e9;        
        constr_func(nvar+1:nvar+pvar) = vpa(qvec) + incl_constr.neg * big_M;  % q>=q0
        constr_func(nvar+pvar+1:nvar+2*pvar) = vpa(-qvec) + incl_constr.pos * big_M;  % q<=q0        

    case 'combine'
        num_of_sets    = 1+pvar;

        optsol.P       = cell(num_of_sets,1);
        optsol.cmap    = cell(num_of_sets,1);
        optsol.epsvar  = cell(num_of_sets,1);
        
        for iq = 1:pvar
            % define objective function weighing factors
            Wobj_i      = 1e-2*ones(size(q));
            Wobj_i(iq)  = 1;

            % solve optimization
            [optsol.P{iq},optsol.epsvar{iq}] = optimization_problem(Ai,cmap,sol_opt,Wobj_i,incl_constr);
        end
        
        % run one final time with Wobj=1
        [optsol.P{end},optsol.epsvar{end}] = optimization_problem(Ai,cmap,sol_opt,ones(size(q)),incl_constr);  
        
        qvec = inv(cmap) * transpose(q-q0);

        % finally: retrieve optimal solution, f_constr>=0
        asymvar = sym('a',[1 num_of_sets]);        
        asymvar(end) = 1-sum(asymvar(1:end-1));     % since SUM_i a_i = 1
        epsvar  = asymvar(1)*optsol.epsvar{1};
        for iSet = 2:num_of_sets
            epsvar  = simplify(epsvar + asymvar(iSet) * optsol.epsvar{iSet});
        end
        
        constr_func = sym(zeros(nvar+2*pvar+num_of_sets,1));
        
        constr_func(1:nvar) = vpa(simplify(1-epsvar*abs(qvec)));
        
        big_M = max(1,norm(q0))*1e9;        
        constr_func(nvar+1:nvar+pvar) = vpa(qvec) + incl_constr.neg * big_M;  % q>=q0
        constr_func(nvar+pvar+1:nvar+2*pvar) = vpa(-qvec) + incl_constr.pos * big_M;  % q<=q0  
        
        for iSet = 1:num_of_sets
            constr_func(nvar+2*pvar+iSet) = asymvar(iSet);      % a_i>=0
        end
    otherwise
        error('unknown algorithm selected!')
end


end


%% function to adjust epsvar by identifying eps0>=1
function epsvar = rescale_epsvar(P,epsvar,A0)

nvar = size(A0,1);

cvx_begin quiet
variable eps0(nvar) nonnegative

subject to
eps0 >= 1
-diag(eps0) - P*A0 - A0'*P == semidefinite(nvar)

minimize(-sum(eps0))

cvx_end

epsvar = epsvar./eps0;

end
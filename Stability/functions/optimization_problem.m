% Optimization problem
function [P,epsvar] = optimization_problem(Ai,cmap,sol_opt,Wobj,incl_constr)

nvar = size(Ai{1},1);   % number of state variables
pvar = length(Wobj);    % number of design parameters

cvx_begin quiet

variable epsvar(nvar,pvar) nonnegative
expression prob_objfun

sdp_constr = cell(1+2*pvar,1);

subject to

% first: nominal system
A1_block = split_block_diagonals(Ai{1});

PblockStr = 'P = blkdiag(';
for iBlock = 1:length(A1_block)
    eval(sprintf('variable P%d(%d,%d) semidefinite', iBlock, ...
        length(A1_block{iBlock}), length(A1_block{iBlock})))
    PblockStr = strcat(PblockStr, sprintf('P%d,',iBlock));
end
eval(strcat(PblockStr(1:end-1),');'))

sdp_constr{1} = -eye(nvar)-(P*Ai{1}+transpose(Ai{1})*P);
sdp_constr{1} == semidefinite(nvar);

Ai(2:end) = transform_Ai(Ai(2:end),cmap);

% next: for the parametric parts
for iq = 1:pvar
    sdp_constr{2*iq}    = diag(epsvar(:,iq)) - incl_constr.pos*(P*Ai{1+iq}+transpose(Ai{1+iq})*P);
    sdp_constr{1+2*iq}  = diag(epsvar(:,iq)) + incl_constr.neg*(P*Ai{1+iq}+transpose(Ai{1+iq})*P); 
    sdp_constr{2*iq}    == semidefinite(nvar);
    sdp_constr{1+2*iq}  == semidefinite(nvar);
end

prob_objfun = norm(epsvar*diag(Wobj),sol_opt.norm);

minimize(prob_objfun)

cvx_end

flag_infeas = 0;

if ~contains(cvx_status,'Solved')
    if contains(cvx_status,'Infeas')
        flag_infeas = 1;
    else
        cvx_status = verify_feasiblity('Prob-I',cvx_status,sdp_constr,{},nvar);
        if ~contains(cvx_status,'Solved')
            flag_infeas = 1;
        end        
    end
end

if flag_infeas
    error('I cannot proceed with this solution')
end

end


%% function to decompose a block-diagonal (BD) matrix into BD components
function A_block = split_block_diagonals(Amat)

% check if block_diagonal
eig_A   = eig(Amat);
diag_A  = diag(Amat);
n_eigs  = length(eig_A);

if max(abs(sort(real(eig_A))-sort(diag_A)))<1e-6
    % perform the decomposition
    A_block = cell(length(uniquetol(diag_A)),1);

    iRow = 1;
    iBlock = 1;
    while iBlock<=length(A_block)
        if iRow==n_eigs || abs(diag_A(iRow)-diag_A(iRow+1))>1e-6 
            % real eigenvalue (non-repeated)
            A_block{iBlock} = Amat(iRow, iRow);

            % advance one row
            iRow = iRow + 1;     
        else    
            % complex eigvalues
            A_block{iBlock} = Amat(iRow:iRow+1, iRow:iRow+1);
            
            % advance two rows
            iRow = iRow + 2;   
            
        end
        % move on to the next block
        iBlock = iBlock + 1;    
    end
else
    A_block{1} = Amat;
end

% verify
if 0
    recov_str = 'A_recov = blkdiag(';
    for iBlock = 1:length(A_block)
        recov_str = strcat(recov_str, sprintf('A_block{%d},',iBlock));
    end
    eval(strcat(recov_str(1:end-1),');'))

    if norm(A_recov-Amat)>1e-9
        error('incorrect block diagonal form!')
    end
end

end
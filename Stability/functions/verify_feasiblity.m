function cvx_status = verify_feasiblity(myprob, cvx_status, sdp_constr, lin_constr, nvar)

fprintf('\n\n[%s] original problem: %s',myprob,cvx_status)

cvx_begin quiet

subject to

for ic = 1:length(sdp_constr)
    sdp_constr{ic} == semidefinite(nvar);
end

for ic = 1:length(lin_constr)
    lin_constr{ic} >= 0;
end

cvx_end

fprintf('\n[%s] verification prob: %s\n',myprob,cvx_status)
pause(0.1)

if ~strcmp(cvx_status,'Solved')
    for ic = 1:length(sdp_constr)
        ic 
        min(eig(sdp_constr{ic}))
    end

    for ic = 1:length(lin_constr)
        ic 
        lin_constr{ic}
    end        
end
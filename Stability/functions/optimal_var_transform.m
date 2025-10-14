%% function to identify optimal variable transformation

function cmap = optimal_var_transform(Ai,P,opt_norm)

nvar = size(Ai{1},1);
pvar = length(Ai);

for ivar = pvar:-1:1
    norm_Ai(ivar)   = norm(full(Ai{ivar}));
    Anorm{ivar}     = Ai{ivar}/norm_Ai(ivar);
end

cvx_begin quiet
variable c(pvar,pvar)
variable epsvar(nvar,pvar) nonnegative

subject to

Ac = transform_Ai(Anorm,c);

for ivar = pvar:-1:1
    c(ivar,ivar) == 1
    for jvar = 1:ivar-1
        c(jvar,ivar) == 0
    end

    Ad{ivar} = P*Ac{ivar} + transpose(Ac{ivar})*P;

    diag(epsvar(:,ivar)) - Ad{ivar} == semidefinite(nvar)
    Ad{ivar} + diag(epsvar(:,ivar)) == semidefinite(nvar)
end

prob_objfun = norm(epsvar,opt_norm);

minimize(prob_objfun)

cvx_end

% define the inverse map
cmap = diag(1./norm_Ai) * c;

if 0 % verify accuracy
    c
    Arec = transform_Ai(Ai,cmap);
    for ivar = 1:pvar
        if norm(full(Arec{ivar}-Ac{ivar}))>1e-9
            error('incorrect variable transformation')
        end
    end

    fprintf('\n *** variable transformation verified! ***\n')
    pause(0.1)
end

end
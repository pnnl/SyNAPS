%% example systems

function [q,q0,new_q,new_q0,Aq,Ai] = extract_model_info(model_opt)

% read model from file
[q,q0,Aq] = example_system(model_opt.sysidx);

% transform, as needed
[temp_q,new_q,new_q0,new_Aq] = transform_parameters(q,q0,Aq,model_opt);

% drop inconsequential modes, if requested
new_Aq = reduce_model(new_Aq,temp_q,new_q0,model_opt.MOR);

Ai = split_matrix(temp_q,new_q0,new_Aq);

end


%% transform to a different set of parameters
function [temp_p,new_q,new_q0,new_Aq] = transform_parameters(q,q0,Aq,model_opt)

% normalize first
abs_q0 = abs(q0);
abs_q0(abs_q0<1e-6) = 1;
Tnorm = pinv(diag(abs_q0));

temp_p = q;
for iq = 1:length(q)
    eval(sprintf('syms new_%s',char(q(iq))))
    eval(sprintf('temp_p(iq) = new_%s;',char(q(iq))))
end

if strcmp(model_opt.xfrm,'no')
    T = Tnorm;
elseif strcmp(model_opt.xfrm,'yes')
    T = model_opt.T * Tnorm;
else
    error('please check the model_opt.xfrm flag!')
end

% define the transformed state-space
new_q   = q * T';
new_q0  = q0 * T';

% now split-up Aq as per the transformed state-space
old_q = temp_p*inv(T)';
new_Aq = vpa(simplify(subs(Aq,q,old_q)));

end

%% transform to a different set of parameters
function Ai = split_matrix(q,q0,Aq)

Ai      = cell(1+length(q),1);
Ai{1}   = double(subs(Aq,q,q0));
for iq = 1:length(q)
    Ai{iq+1}    = double(subs(diff(Aq,q(iq)),q,q0));
end

if max(real(eig(Ai{1})))>=0
    warning('nominal system not stable!')
end

% verify
recov_A    = Ai{1};
for iq = 1:length(q)
    recov_A = recov_A + Ai{iq+1} * (q(iq)-q0(iq));
end

if norm(double(subs(Aq-recov_A,q,q0)))>1e-9
    % sanity check
    Aq
    recov_A
    norm(double(subs(Aq-recov_A,q,q0)))
    pause
end

end


%% reduce dimension of A by discarding well-damped eigenvalues
function Ared = reduce_model(Aq,q,q0,opt_MOR)

A0 = double(subs(Aq,q,q0));
[V,D] = eig(A0);
eig_A = diag(D);

% identify pairs of complex eigenvalues
idx_complex_eig = find(imag(eig_A)~=0);

VV = V;
for iEig = 1:length(idx_complex_eig)/2
    idx1 = idx_complex_eig(2*iEig-1);
    idx2 = idx_complex_eig(2*iEig);
    VV(:,idx1) = real(V(:,idx1));
    VV(:,idx2) = imag(V(:,idx2));
end

% now block-diagonalize
T_mat = inv(VV);       % such that: D_Block = T_mat * A0 * inv(T_mat)
diag_Aq = vpa(T_mat * Aq * VV);

if 1
    % sanity check: compare eigvalues before/after diagonalization
    eig1 = real(eig_A);
    eig2 = real(diag(double(subs(diag_Aq,q,q0))));

    if max(abs(eig1-eig2))>1e-6
        error('incorrect diagonalization')
    end

end

if strcmp(opt_MOR,'yes')
    max_eig = max(real(eig_A));

    if max_eig>=0
        error('cannot apply MOR to unstable matrix')
    end

    % identify the dominant eigenvalues
    dom_idx = find(real(eig_A)>=40*max_eig); % use: 20

    if 0
        % further reduce based on eigenvalue sensitivity
        max_sen = find_eigenvectors(Aq,q,q0,A0,eig_A);

        % identify the sensitive eigenvalues
        sen_idx = find(max_sen>=20*max(max_sen)); % use: ?
    else
        sen_idx = [];
    end

    % identify the eigenvalues to keep
    keep_idx = union(dom_idx,sen_idx);
    drop_idx = setdiff(1:length(eig_A),keep_idx);

    A11 = diag_Aq(keep_idx,keep_idx);
    A12 = diag_Aq(keep_idx,drop_idx);
    A21 = diag_Aq(drop_idx,keep_idx);
    A22 = double(subs(diag_Aq(drop_idx,drop_idx),q,q0));

    Ared = A11 - A12 * (A22\A21);

    if 1
        % sanity check
        eig1 = sort(eig_A(keep_idx));
        eig2 = sort(eig(double(subs(Ared,q,q0))));

        if max(abs(eig1-eig2))>1e-6
            error('invalid model reduction')
        end
    end

elseif strcmp(opt_MOR,'no')
    Ared = diag_Aq;

else    
    error('the model_opt.MOR flag must be either ''yes'' or ''no''!')
end

end

% find right and left eigenvectors
function max_sen = find_eigenvectors(Aq,q,q0,A0,eig_A)

[RV,RD] = eig(A0);
[LV,LD] = eig(transpose(A0));

dA = Aq - A0;

RE = diag(RD);
LE = diag(LD);
eig_sen = sym(NaN(size(A0,1),length(q)));

for iE = 1:length(eig_sen)
    lv = transpose(LV(:,abs(LE-eig_A(iE))<1e-9));
    rv = RV(:,abs(RE-eig_A(iE))<1e-9);
    eig_sen(iE,:) = vpa(simplify(jacobian( (lv*dA*rv)/(lv*rv), q)));
end

eig_sen = double(subs(eig_sen,q,q0));

rel_sen = abs(eig_sen)./kron(abs(real(eig_A)),ones(1,length(q)));
max_sen = max(rel_sen,[],2);

end
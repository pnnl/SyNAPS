clear all; close all; clc;

q1_sym  = sym('k4','real');
q2_sym  = sym('k5','real');
q3_sym  = sym('kpdc','real');
q4_sym  = sym('kidc','real');

q   = [q1_sym q2_sym q3_sym q4_sym];

% load MTDC network system matrices
load('A_ACDC.mat','A')

q1_val = 1;
q2_val = 1;
q3_val = 0.1;
q4_val = 10;

q0 = [q1_val q2_val q3_val q4_val];

orig_Aq = sym(A);

% add q1
Avec_q1 = [0; A(2,4); 0; -0.201; 0; A(6,4); zeros(4,1)]*(q1_sym-q1_val);
orig_Aq(:,4) = orig_Aq(:,4) + Avec_q1;

% add q2
Avec_q2 = [0; A(2,6); 0; A(4,6); 0; 0.0475; zeros(4,1)]*(q2_sym-q2_val);
orig_Aq(:,6) = orig_Aq(:,6) + Avec_q2;

% add q3
Avec_q3 = [zeros(8,1); -0.01; 0]*(q3_sym-q3_val);
orig_Aq(:,9) = orig_Aq(:,9) + Avec_q3;

% add q4
Avec_q4 = [zeros(8,1); A(9,10); 0]*(q4_sym/q4_val-1);
orig_Aq(:,10)= orig_Aq(:,10) + Avec_q4;

orig_A0 = double(subs(orig_Aq,q,q0));
eig(orig_A0)

% now do transformation to express in relative angles
ref_idx = 1;
rel_idx = [3 5];

T = eye(length(orig_A0));
T(rel_idx,ref_idx) = -1;
AqT = T*orig_Aq*inv(T);

% discard the reference angle state
Aq = AqT(setdiff(1:length(AqT),ref_idx),setdiff(1:length(AqT),ref_idx));

% verify
A0 = double(subs(Aq,q,q0));

eig(A0)

% save data
save('example_acdc5bus.mat','Aq','q','q0')
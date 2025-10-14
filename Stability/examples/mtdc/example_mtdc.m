clear all; close all; clc;

q1_sym  = sym('q1','real');
q2_sym  = sym('q2','real');
q3_sym  = sym('rscale','real');

q   = [q1_sym q2_sym q3_sym];

% load MTDC network system matrices
load('A1.mat','A1_mat')
load('A2.mat','A2_mat')

q1_val = abs(A2_mat(end-1,end-1));
q2_val = abs(A2_mat(end-1,end));
q3_val = 1;

q0 = [q1_val q2_val q3_val];

A2_sym = sym(A2_mat);
A2_sym(end-1,end-1) = q1_sym * sign(A2_mat(end-1,end-1));
A2_sym(end-1,end)   = q2_sym * sign(A2_mat(end-1,end));

Aq = A1_mat * q3_sym + A2_sym;

A = double(subs(Aq,q,q0));

eig(A)

% save data
save('example_mtdc.mat','Aq','q','q0')
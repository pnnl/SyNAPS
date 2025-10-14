close all; clear all; clc

mod_opt = 'noload';   % OPTIONS: {'load','noload'}

% R/X ratio
% ratio = [a b c];

ratio = 1.3;

% Drhoop Gains (%)
kpval = 1;          % 1.1881;
kqval = 0.1;        % kpval/0.3;

syms kp1 kq1 kp2 kq2 kp3 kq3 kp4 kq4
qvar = [kp1 kq1 kp2 kq2 kp3 kq3 kp4 kq4];%[kp, kq];

qval = [kpval kqval kpval kqval kpval kqval kpval kqval];

% choose a subset
qvar(3:end) = qval(3:end);

% Base Peak Phase Voltage (V)
U_b = 230;
% Base Inverter Apparent Power
S_b = 10e3/3;
% Nominal Frequency (rad)
omega_0 = 2*pi*50;
% Filter Cut off Frequency (rad)
omega_c = 10*pi;
% Line Impedance (Ohm/km)
Z = 0.1587*(ratio+1i);
% Z_base (Ohm)
Z_b = U_b^2/S_b;
% Line Impedance (pu/km)
Z = Z/Z_b;

% Kundur system line distances(km)
l = [6; 100; 3];
Z_line = Z.*l;

switch mod_opt
    case 'load'     % include load
        % Load impedances (Ohm)
        Z_load = [20+1j; 25+1j;20+4.75j;40+12.58j];
        Z_load = conj(Z_load)*1j;

        % Load impedsances (p.u.)
        Z_load = Z_load./Z_b;
    
    case 'noload'   % ignore load    
        Z_load = [];

    otherwise
        error('unknown model option')
end

% Combine all impedances into one set
Z_set = [Z_line;Z_load];

Y_n = [1/Z_set(1) -1/Z_set(1) 0 0;
    0 1/Z_set(2) -1/Z_set(2) 0
    0 0 1/Z_set(3) -1/Z_set(3)];
Y_l = diag(1./Z_load);

Y = [Y_n; Y_l];

% Susceptance and Conductance matrices
B = -imag(Y);
G = real(Y);

% Incidence matrix
tol = 1e-6;
del = transpose((B>tol) - (B<-tol));

[n,m] = size(del);

% Droop Gains
K_p = qvar(1:2:end);%*ones(n,1);
K_q = qvar(2:2:end);%*ones(n,1);

k_p = omega_0*diag(K_p./100);
k_q = diag(K_q./100);

% The state space matrix for Electro-Magnetic model
% s*x = Asym_all*x, where x = (theta, omega, V, I_d, I_q)
Asym_all   = sym(zeros(3*n+2*m));
% set of indixes for five state variables x = (theta, omega, V, I_d, I_q)
theta   = 1:n;
omega   = n+1:2*n;
V       = 2*n+1:3*n;
I_d     = 3*n+1:3*n+m;
I_q     = 3*n+m+1:3*n+2*m;
% The state matrix construction
% trivial equation for s*theta = omega
Asym_all(theta, omega) = eye(n);
% equation of the Low Pass filter s*omega = omega_c*omega
Asym_all(omega, omega) = -omega_c*eye(n);
Asym_all(omega, I_d) = -omega_c*k_p*del;
% equation of the Low Pass filter s*V = omega_c*V
Asym_all(V, V) = -omega_c*eye(n);
Asym_all(V, I_q) = omega_c*k_q*del;
Asym_all(I_d, V) = omega_0*(ratio^2+1)*B;
Asym_all(I_d, I_d) = -omega_0*ratio*eye(m);
Asym_all(I_d, I_q) = omega_0*eye(m);
Asym_all(I_q, theta) = omega_0*(ratio^2+1)*B;
Asym_all(I_q, I_d) = -omega_0*eye(m);
Asym_all(I_q, I_q) = -omega_0*ratio*eye(m);


% reduced A (removing line and load dynamics)
Asym_red = Asym_all(1:3*n,1:3*n);
Asym_red(omega, theta) = Asym_red(omega, theta) - omega_c*k_p*del*B;
Asym_red(omega, V)     = Asym_red(omega, V) - omega_c*ratio*k_p*del*B;
Asym_red(V, theta)     = Asym_red(V, theta) + omega_c*ratio*k_q*del*B;
Asym_red(V, V)         = Asym_red(V, V) - omega_c*k_q*del*B;

%% Do a similarity transformation
Tall = eye(3*n+2*m);
Tall(2:n,1) = -1;

Tred = eye(3*n);
Tred(2:n,1) = -1;

ATsym_all   = Tall * (Asym_all/Tall) ;
ATsym_red   = Tred * (Asym_red/Tred) ;

ATsym_all   = ATsym_all(2:end,2:end);
ATsym_red   = ATsym_red(2:end,2:end);

eig1 = eig(double(subs(Asym_all,qvar,qval)));
eig2 = eig(double(subs(Asym_red,qvar,qval)));
eig3 = eig(double(subs(ATsym_all,qvar,qval)));
eig4 = eig(double(subs(ATsym_red,qvar,qval)));

fprintf('\n:===: PRINT DOMINANT EIGENVALUE (REAL PART) :===:\n')
fprintf('\n full\t\t\t\t: %f',max(real(eig1)))
fprintf('\n reduced\t\t\t: %f\n',max(real(eig2)))
fprintf('\n full effective (trans)\t\t: %f',max(real(eig3)))
fprintf('\n reduced effective (trans)\t: %f\n',max(real(eig4)))
    

figure
plot(real(eig1),imag(eig1),'+','markersize',8,'linewidth',1.5)
hold on
plot(real(eig3),imag(eig3),'+','markersize',8,'linewidth',1.5)
plot(real(eig2),imag(eig2),'+','markersize',8,'linewidth',1.5)
plot(real(eig4),imag(eig4),'+','markersize',8,'linewidth',1.5)
legend('full','transformed','reduced','reduced tranformed')


save(sprintf('example_vorobev_%s.mat',mod_opt),'Asym_all','Asym_red', ...
                                    'ATsym_all','ATsym_red','qvar');
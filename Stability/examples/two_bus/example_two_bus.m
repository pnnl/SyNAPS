clear all; close all; clc;

mp_sym      = sym('mp','positive');
mq_sym      = sym('mq','positive');
kpv_sym     = sym('kpv','positive');
kiv_sym     = sym('kiv','positive');
load_sym    = sym('load','positive');

mp_val      = 0.01;
mq_val      = 0.05;
kpv_val     = 1;
kiv_val     = 5.86;
load_val    = 2;

% choose which are symbolic vs. numerical variables
mp_var      = mp_val;
mq_var      = mq_sym;
kpv_var     = kpv_val;
kiv_var     = kiv_val;
load_var    = load_sym;

q   = [load_sym mq_sym];
q0  = [load_val mq_val];

% two-bus GFM-load network

base_params.w = 120*pi;

inv_params = struct('tau',  0.01,   ...
                    'mp',   mp_var,   ...
                    'mq',   mq_var,   ...
                    'R',    0.01,   ...
                    'X',    0.1,    ...
                    'kpv',  kpv_var,  ...   % 1
                    'kiv',  kiv_var);       % 5.86

net_params = struct('Zload', load_var*(1+0.1i), ...
                    'Iload', -2*(1+0.1i));

% ------------ P-f droop ------------
% d_dot = w-w0
% w_dot = 1/tau * [ w0-w  + m_p  * (P_set - P_grid) ]

% ------------ Q-V droop (PI-Loop) ------------
% Ve_dot = 1/tau * [ V_set - Ve - V + m_q * (Q_set - Q_grid) ]
% E_dot  = 1/tau * [ k_pv * (V_set - V) - (k_pv - tau * k_iv) * Ve + k_pv * m_q * (Q_set - Q_grid) ]

d = sym('d','real');
w = sym('w','real');
e = sym('e','real');
v = sym('v','real');

xsym = [d; w; e; v];

% solve the network flow
i_load = conj(net_params.Iload);
z_load = conj(1/net_params.Zload);
z_conn = inv_params.R + inv_params.X*1i;   

factor.a = z_load/(z_conn+z_load);
factor.b = z_conn*factor.a*i_load;

vnet_real = v * (cos(d)*real(factor.a) - sin(d)*imag(factor.a)) - real(factor.b);
vnet_imag = v * (cos(d)*imag(factor.a) + sin(d)*real(factor.a)) - imag(factor.b);
vnet = sqrt( vnet_real^2 + vnet_imag^2 );
Pnet = vnet^2 * real(net_params.Zload) + vnet_real*real(i_load) + vnet_imag*imag(i_load);
Qnet = vnet^2 * imag(net_params.Zload) + vnet_real*imag(i_load) - vnet_imag*real(i_load);

% equilibrium value
vnom = 1 + z_conn * (1/z_load + i_load);

set_points = struct('w0',   0, ...
                    'Vs',   1, ...
                    'Ps',   real(net_params.Zload)+real(net_params.Iload), ...
                    'Qs',   imag(net_params.Zload)+imag(net_params.Iload));

% d_dot
fsym(1) = (w - set_points.w0) * base_params.w;

% w_dot
fsym(2) = 1/inv_params.tau * ( set_points.w0-w  + inv_params.mp  * (set_points.Ps - Pnet) );

% e_dot
fsym(3) = 1/inv_params.tau * ( set_points.Vs - vnet - e + inv_params.mq * (set_points.Qs - Qnet) );

% v_dot
fsym(4) = inv_params.kpv * fsym(3) + inv_params.kiv * e;


% analyze
xnom_sym = [ angle(vnom) set_points.w0 0 abs(vnom) ]';
Aq = vpa(subs(jacobian(fsym,xsym),xsym,xnom_sym));

A = double(subs(Aq,q,q0));

eig(A)

% save data
save('example_two_bus.mat','Aq','q','q0')
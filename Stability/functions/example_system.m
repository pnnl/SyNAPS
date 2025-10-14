%% example systems

function [q,q0,Aq] = example_system(sysidx)

switch sysidx
    case 1
        syms q1 q2
        q       = [q1 q2];
        q0      = [-1 -2];%[1 1.5]; %zeros(size(q));   % [-3 -3];
        Aq      = [-2+q1      0     -1+q1; 
                    0       -3+q2     0; 
                    -1+q1   -1+q2   -4+q1];
    case 2
        syms q1 q2 q3 q4
        q       = [q1 q2 q3 q4];
        q0      = [-3 -3 3 -3]; %[-2 -0.035 2.7 -0.002];
        Aq      = [ q1      0       0;
                    0       q2      q3;
                    0   -0.7115     q4];   

    case 3
        syms q1 q2
        q       = [q1 q2];
        q0      = [-5 -5];%[-1.25 -3.5];
        Aq      = [ q1   -12.06 -0.06 0;
                    -0.25 -0.03 1.00 0.5;
                    0.25 -4.00 -1.03 0;
                    0 0.50 0 q2];

    case 4
        load('./examples/example_vorobev_noload.mat','ATsym_all','qvar');
        q       = symvar(qvar);
        q0      = [3 3];%[1 0.1];
        Aq      = ATsym_all;

    case 5
        syms q1 q2
        q       = [q1 q2];
        q0      = [-1 2];
        Aq      = [q1 0.1; 0 1-q2];

    case 6
        load('./examples/example_two_bus.mat','Aq','q','q0');
        q0 = [1.5 0.03];

    case 7
        load('./examples/example_oliviera.mat','Aq','q','q0');

    case 8
        load('./examples/example_mtdc.mat','Aq','q','q0'); 

    case 9 
        load('./examples/example_acdc5bus.mat','Aq','q','q0'); 
        % q0(1:2) = [0.01 0.01];

    case 10
        load('./examples/example_acdc14bus.mat','Aq','q','q0'); 
        % q0(1:2) = [2 0.01];

    otherwise
        error('example not available')
end

end
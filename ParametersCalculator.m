%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                Parameters Calculator                        %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
%
% This code claculate the parameters related to the family of-
% - CNe methods. To better Understand those parameters and their-
% - physical meaning, please check the publications of Dr.Endre-
% -Kovacs from University of Miskolc
%******************************************************************
% INPUt

% M1                  % Matrix System
% N                   % Number of Nodes
% dt                  % time Step Size


% OUTPUT

% ee, ee1, b: to understand those parameters, one should read the formula-
% of the methods in the published paper.

%******************************************************************

function [b, ee, emid, ee1, emid1]  = ParametersCalculator(M1, N, dt)

%=========================================================
% pre-allocation:

b     = zeros(N,1);
ee    = zeros(N,1);
emid  = zeros(N,1);    % mid = midpoint, at h/2
ee1   = zeros(N,1);
emid1 = zeros(N,1);    % mid = midpoint, at h/2

%=========================================================
% calculations:

for i=1:N
    b(i)    = -M1(i,i);
    bh      = b(i)*dt;
    ee(i)   = exp(-bh);
    emid(i) = exp(-bh/2);

    % Now we define 1-ee:
    ee1(i)   = 1-ee(i);
    emid1(i) = 1-emid(i);
    % This is good if bh is not too small
    
    % but if bh is small:
    if(bh<0.01)
        ee1(i)   = 0;
        emid1(i) = 0;
        for s=1:8
            ee1(i)   = ee1(i)-(-1)^s*bh^s/factorial(s); %% power series of 1-exp(-bh)
            emid1(i) = emid1(i)-(-0.5)^s*bh^s/factorial(s); %% the same but at the middle of the time step
        end
    end
end
end
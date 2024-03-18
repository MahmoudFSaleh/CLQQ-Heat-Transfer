%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                 Stiffness Calculator                        %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
%
% This function does the following: 1- Calculate the Stiffness-
%                                      -Ratio of the System
%                                   2- Calculate the Max Time-Step
%******************************************************************
% INPUt
% M1                  % The matrix of the ODE system

% OUTPUT

% StiffRatio             % The Stiffness Ratio of the System
% dtMAX                  % Maximum time-step for Explicit Euler

%******************************************************************

function [StiffRatio, dtMAX] = StiffnessCalculator(M1, N)

E    = eig(M1);  %%Eigenvalues
Emax = 0;
E0   = -1;
Emin = -2;

for i=1:N
    if E(i)<Emax
        Emax = E(i);
    end
    if E(i)>E0
        E0 = E(i);
    end
end
for i=1:N
    if E(i)>Emin
        if E(i)<E0
            Emin = E(i);
        end
    end
end
dtMAX      = -2/Emax;  %% max timestep for Explicit Euler
StiffRatio = Emax/Emin;
end
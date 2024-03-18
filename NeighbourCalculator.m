%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                Neighbour Calculator                        %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:s
%
% This Code will calculate the four cell neighbours of each cell,-
% - which are the left, right, upper and lower neighbours. Those-
% neighbours are calculated by the matrix system M1. Actually, one-
% -can perform the calculation of the numerical method even without-
%- using those neighbour vectors, but those neighbour vectors simpl-
%-ify the calculations. It is similar to that one in SPH method

%******************************************************************
% INPUt

% M1                  % Matrix System
% Nx                  % Number of Nodes in x Direction
% Nz                  % Number of Nodes in z Direction

% OUTPUT

% Neb1:               % The left neighbour
% Neb2:               % The right neighbour
% Neb3:               % The Upper neighbour
% Neb1:               % The Lower neighbour


%******************************************************************
% The Function:

function [Neb1, Neb2, Neb3, Neb4] = NeighbourCalculator(M1, Nx, Nz)

%=========================================================
% pre-allocation:

N    = Nx*Nz;
Neb1 = zeros(N, 1);
Neb2 = zeros(N, 1);
Neb3 = zeros(N, 1);
Neb4 = zeros(N, 1);


%=========================================================
% calculations:

for i=1:N
    rx=rem(i,Nx);
    if (rx>1 || rx<1)   % left neighbor
        Neb1(i)=M1(i,i-1);
    end
    if (rx>0)           % right neighbor
        Neb2(i)= M1(i,i+1);
    end
    if (i>Nx)           % upper neighbor
        Neb3(i)=M1(i,i-Nx);
    end
    if (i<=Nx*(Nz-1))    % lower neighbor
        Neb4(i)=M1(i,i+Nx);
    end
end
end
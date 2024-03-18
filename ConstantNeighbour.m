%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                Constant Neighbour Method                    %%%
%%%                      or CNe Method                          %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
%
% This code is the time integrator based on the Constant Neighbour-
% Method. This method is explained in the publication of Dr. Endre-
% Kovacs from University of Miskolc.
%
%******************************************************************
% INPUt

% Neb1:               % The left neighbour
% Neb2:               % The right neighbour
% Neb3:               % The Upper neighbour
% Neb1:               % The Lower neighbour
% Nx                  % Number of Nodes in x Direction
% Nz                  % Number of Nodes in z Direction
% U                   % The temperature at the beginning of the time step
% ee, ee1, b: to understand those parameters, one should read the formula-
% of the methods in the published paper.

% OUTPUT

% UTEMP         % The Temperature calculated by CNe method
% a0            %  to understand this parameter,-
%                 -Check the papers related to CNe method

%******************************************************************

function [UTEMP, a0]  = ConstantNeighbour(Nx, Nz, Neb1, Neb2, Neb3, Neb4, U, ee, ee1, b, t)

N     = Nx*Nz;
a0    = zeros(N,1);
UTEMP = zeros(N,1);

for i=1:N  %%%Neighbours are constans
    rx    = rem(i,Nx);

    if (rx>1 || rx<1)     % left neighbor
        a0(i) = a0(i) + Neb1(i)*U(i-1);
    end  
    if (rx>0)             % right neighbor
        a0(i) = a0(i) + Neb2(i)*U(i+1);
    end 
    if (i>Nx)             % upper neighbor
        a0(i) = a0(i) + Neb3(i)*U(i-Nx);
    end  
    if (i <= Nx*(Nz-1))   % lower neighbor
        a0(i) = a0(i) + Neb4(i)*U(i+Nx);
    end 
    
    UTEMP(i) = U(i)*ee(i) + a0(i)/b(i)*ee1(i);
end


end
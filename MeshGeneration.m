%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                 Generate the Mesh                           %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
% 
% This function does the following: 1- Divide the Space
%                                   2- Set the Properties of the Space
%                                   3- Set the Initial Temperature
%                                   4- Calculate the matrix M1
%******************************************************************
% INPUt
% dx                  % Mesh Size in x Direction
% dz                  % Mesh Size in z Direction
% Nx                  % Number of Nodes in x Direction
% Nz                  % Number of Nodes in z Direction

% OUTPUT

% Node.C              % The Heat Capacity
% Node.Rx             % Resistance in x Direction
% Node.Rz             % Resistance in z Direction
% Node.x              % Position in x Direction
% M1                  % The matrix of the ODE system
% U0                  % Initial Consitions
%******************************************************************

function [Node, M1, U0] = MeshGeneration(dx, dz, Nx, Nz, MaterialOption)

%******************************************************************
% pre-allocation:

N       = Nx*Nz;
Node.C  = zeros(N,1);            % The Heat Capacity
Node.Rx = zeros(N,1);            % Resistance in x Direction
Node.Rz = zeros(N,1);            % Resistance in z Direction
Node.x  = zeros(N,1);            % Position in x Direction
U0      = zeros(N,1);            % Initial Consitions
M1      = zeros(N, N);           % The matrix of the ODE system

%******************************************************************
% Divide the space:
for i=1:1:Nx
    Node.x(i) = (i-1)*dx;
end
%******************************************************************
% Set the Initial Conditions
for i=1:1:N
    U0(i) = 100*(1 - Node.x(i))*Node.x(i);
end

%******************************************************************
% Set the Properties

switch MaterialOption

    case 1
        % Diffusivity is 1/4
        for i=1:1:N
            Node.Rx(i) = 0.5;
            Node.Rz(i) = 0.5;
            Node.C(i)  = (dx^2)/((1/4)*Node.Rx(1));
        end

    case 2
        % The Property has not been defined yet
         warndlg([' The second option of material Property has not been set. It can affect' ...
            ' the analytical solution. Please check function "MeshGeneration" >> switch' ...
            'That value is calculated in MeshGeneration function'],'Warning')
end

Node.C(1)  = (1/2)*Node.C(5);
Node.C(N)  = (1/2)*Node.C(5);  % The capacity of the cells of the Borders 
%                             is half the capacity of inner cells

%******************************************************************
% Construct the Matrix of the ODE system

for i=1:N
    rx=rem(i,Nx);
    if (rx>1 || rx<1)                                % left neighbor
        M1(i,i-1) = 1/Node.Rx(i-1)/Node.C(i);
        M1(i,i)   = -M1(i,i-1);
    end
    if (rx>0)                                        % right neighbor
        M1(i,i+1) = 1/Node.Rx(i)/Node.C(i);
        M1(i,i)   = M1(i,i) - M1(i,i+1);
    end
    if (i>Nx)                                         % upper neighbor
        M1(i,i-Nx) = 1/Node.Rz(i-Nx)/Node.C(i);
        M1(i,i)    = M1(i,i) - M1(i,i-Nx);
    end
    if (i<=Nx*(Nz-1))                                 % lower neighbor
        M1(i,i+Nx) = 1/Node.Rz(i)/Node.C(i);
        M1(i,i)    = M1(i,i) - M1(i,i+Nx);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                Reference Solution                           %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
%
% This function does the following: 1- Produce a Reference 
%                                      Solution Analytically
%******************************************************************
% INPUt
% PositionofNodes           % The Position of the Nodes
% Nx                        % Number of Nodes in x Direction
% tfin                      % Final Time

% OUTPUT
% BenchMarlSolution           % Reference Solution

%******************************************************************

function [BenchMarkSolution] = ReferenceSolution(x, Nx, tfin)
e        = exp (1);
solution = zeros(Nx,1);
UMAT     = zeros(Nx,1);

for j=1:1:10000
    for i=1:1:Nx
        term_1 = (1 + (-1)^j)/(j^2);
        term_2 = e^(-((j*pi)^2)*tfin/(4));
        term_3 = cos(j*pi*x(i));
        solution(i,1) = term_1*term_2*term_3;
    end
    UMAT(:,1) = UMAT(:,1)+ solution(:,1);
end
UMAT = (50/3)-(200/(pi^2))*UMAT;
BenchMarkSolution = UMAT(:,1);

end
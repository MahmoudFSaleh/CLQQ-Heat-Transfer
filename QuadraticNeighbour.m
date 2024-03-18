%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                Quadratic Neighbour Method                   %%%
%%%                      or Q Method                            %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
%
% This code is the time integrator based on the Quadratic Neighbour-
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
% UTEMP               % The temperature at the end of the time step-
%                       -calculated by LNe
% Umid                % The temperature at dt/2 of the time step-
%                     % calcilated by LNe method
% ee, ee1, b: to understand those parameters, one should read the formula-
% of the methods in the published paper.

% OUTPUT

% QUTEMP            % The Temperature calculated by Q method at the end-
%                     -of the time step
% QUmid             %  The Temperature calculated by LNe method at dt/2-
%                   % by Q method


%******************************************************************

function [QUmid, QUTEMP]  = QuadraticNeighbour(Nx, Nz, Neb1, Neb2, Neb3, Neb4, Umid, UTEMP, U, emid, emid1, ee, ee1, a0, b, t, dt)

N      = Nx*Nz;
QUmid  = zeros(N,1);
QUTEMP = zeros(N,1);

for i=1:N   %%%Neighbours are quadratic
    anew = 0;
    amid = 0;
    rx   = rem(i,Nx);
    if (rx>1 || rx<1)              % Left Neighbour
        anew = anew + Neb1(i)*UTEMP(i-1);
        amid = amid + Neb1(i)*Umid(i-1);
    end
    if (rx > 0)                     % Right Neighbour
        anew = anew + Neb2(i)*UTEMP(i+1);
        amid = amid + Neb2(i)*Umid(i+1);
    end
    if (i > Nx)                       % Upper Neighbour
        anew = anew + Neb3(i)*UTEMP(i-Nx);
        amid = amid + Neb3(i)*Umid(i-Nx);
    end
    if (i <= Nx*(Nz-1))              % Lower Neighbour
        anew = anew + Neb4(i)*UTEMP(i+Nx);
        amid = amid + Neb4(i)*Umid(i+Nx);
    end
    a1 = 2*(amid - a0(i))/dt;
    w  = 2*((anew - a0(i))/dt - a1)/dt;
    s  = a1 - w*dt/2;
    QUTEMP(i) = U(i)*ee(i) + (ee1(i)*((2*w/b(i) - s)/b(i) + a0(i)) + dt*(w*(dt - 2/b(i)) + s))/b(i);
    QUmid(i)  = U(i)*emid(i) + (emid1(i)*((2*w/b(i) - s)/b(i) + a0(i)) + dt/2*(w*(dt/2 - 2/b(i)) + s))/b(i);
end


end
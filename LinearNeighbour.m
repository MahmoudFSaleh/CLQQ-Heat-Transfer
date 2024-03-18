%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                Linear Neighbour Method                      %%%
%%%                      or LNe Method                          %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
%
% This code is the time integrator based on the Linear Neighbour-
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
%                       -calculated by CNe
% ee, ee1, b: to understand those parameters, one should read the formula-
% of the methods in the published paper.

% OUTPUT

% LUTEMP            % The Temperature calculated by LNe method at the end-
%                     -of the time step
% LUmid             %  The Temperature calculated by LNe method at dt/2-
%                   % by LNe method


%******************************************************************

function [LUmid, LUTEMP]  = LinearNeighbour(Nx, Nz, Neb1, Neb2, Neb3, Neb4, UTEMP, U, emid, emid1, ee, ee1, a0, b, t, dt)

N      = Nx*Nz;
LUmid   = zeros(N,1);
LUTEMP = zeros(N,1);

   for i=1:N  
         anew = 0; 
         rx   = rem(i,Nx);
         tau  = 1/b(i);

     if (rx>1 || rx<1)                       % left neighbour
         anew = anew + Neb1(i)*UTEMP(i-1); 
     end               
     if (rx > 0)                             % right neighbour
         anew = anew + Neb2(i)*UTEMP(i+1); 
     end            
     if (i > Nx)                             % upper neighbou
         anew = anew + Neb3(i)*UTEMP(i-Nx); 
     end  
     if (i <= Nx*(Nz-1))                     % lower neighbour
         anew = anew + Neb4(i)*UTEMP(i+Nx); 
     end 

     a1       = (anew-a0(i))/dt;
     LUTEMP(i) = U(i)*ee(i) + (a0(i) - a1*tau)*ee1(i)*tau + a1*dt*tau;
     LUmid(i)  = U(i)*emid(i) + (a0(i)- a1*tau)*emid1(i)*tau + a1*dt/2*tau;
   end

end
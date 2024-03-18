%%  Section 0:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%             Analytical and Numerical Solution of            %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
%%%                      U_{t} = (1/4)*U_{xx}                   %%%
%%%  Spatial Domain      x ∈ [0, 1]                             %%%
%%%  Boundry Conditions  U(0, t)_{x} = U(1, t)_{x} = 0          %%%
%%%  Initial Conditions  U(x, 0)     = 100*(1-x)*x              %%%
%******************************************************************
% Description:
%
% This code solves the heat equation numerically using CLQQ method
% and compare the result with the analytical solution
%******************************************************************
clear;
clc;
close all;
%%%% End of Section 0:
%% Section 1: Describe the Domain and Introduce the Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%******************************************************************
% The Data of the Geometry
%
Lx  = 1;                 % plate length on x-axis
Lz  = 0;                 % plate length on z-axis
Nx  = 2001;              % Number of Nodes in x Direction
Nz  = 1;                 % Number of Nodes in z Direction
N   = Nx*Nz;             % Total Number of Nodes
dx  = Lx/(Nx - 1);       % Mesh Size in x Direction
dz  = Lz/(Nz - 1);       % Mesh Size in z Direction


%******************************************************************
% Material Property
%
% You can choose on of the material options. Here I introduce
% only two options. Anyway, you can set as many options you want.
% Go to the function "MeshGeneration" where the property of the-
% material is defined.
% 1- MaterialOption = 1: means that the diffusivity = 1/4
% 1- MaterialOption = 2: means that the diffusivity = 1

MaterialOption = 1;
% MaterialOption = 2;
%******************************************************************
% Time Specification
t0   = 0;              % Initial Time
tfin = 0.03;           % Final Time

%%%% End of Section 1:
%% Section 2: Generate the Domain, Create the Mesh, Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Argument "Node" is a structure variable Containing Information
% about the the Capacity, Resistance...etc
% M1: is the matrix system
% U0: is the initial temperature

[Node, M1, U0] = MeshGeneration(dx, dz, Nx, Nz, MaterialOption);

%%%% End of Section 2:
%% Section 3: Test the Calculations, Calculate the Stiffness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%******************************************************************
% Test

NumericalDiffusivity = (dx^2)/(Node.Rx(5)*Node.C(5));
disp(['The Value of the Numerical Diffusivity',' is ' ...
    ,num2str(NumericalDiffusivity)])
warndlg([' Check if the value of the Numerical Diffusivity is ' ...
    'equal to the theoretical one. ' ...
    'The value is printed in the command window'],'Important')


%******************************************************************
% Calculate the Stiffness of the ODE system
% and the maximum time-step size for explicit Euler

[StiffRatio, dtMAX] = StiffnessCalculator(M1, N);


%%%% End of Section 3:
%% Section 4: Benchmark Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%******************************************************************
% Solve the Equation
x = Node.x;
[BenchMarklSolution] = ReferenceSolution(x, Nx, tfin);


%%%% End of Section 4:
%% Section 5: CLQQ Numerical Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%******************************************************************
% pre-allocation:
Mh       = 16;                     % Number of Time Steps
Axhstep  = zeros(Mh,1);            % Time Step Size
MaxD     = zeros(Mh,1);            % Maximum Error
RunT     = zeros(Mh,1);            % Running Time

%******************************************************************
% Big Loop

for ih=1:1:Mh

    %=========================================================
    % pre-calculations ans initilization:

    dt             = (tfin)/(2^(ih));     % time step size
    Axhstep(ih, 1) = dt;

    [b, ee, emid, ee1, emid1] = ParametersCalculator(M1, N, dt);
    [Neb1, Neb2, Neb3, Neb4]  = NeighbourCalculator(M1, Nx, Nz);


    t           = t0;
    TimeIndex   = 1;             % index to Control the time History

    U               = U0;        % Initial Teperature


    %=========================================================
    % Time Integration:
    tic
    while (t < tfin)

        if ((tfin - t) <= dt)    % this will be executed maximum one time
            dt = min(dt, tfin - t);
            [b, ee, emid, ee1, emid1] = ParametersCalculator(M1, N, dt);
        end
        %-----------------------------------------------------------
        % CLQQ stages:


        [UTEMP, a0] = ...   % CNe stage
            ConstantNeighbour(Nx, Nz, Neb1, Neb2, Neb3, Neb4, U, ee, ee1, b, t);
        [Umid, UTEMP] = ...  % LNe stage
            LinearNeighbour(Nx, Nz, Neb1, Neb2, Neb3, Neb4, UTEMP, U, emid, emid1, ee, ee1, a0, b, t, dt);
        [Umid, UTEMP] = ...  % Q1 stage
            QuadraticNeighbour(Nx, Nz, Neb1, Neb2, Neb3, Neb4, Umid, UTEMP, U, emid, emid1, ee, ee1, a0, b, t, dt);
                [Umid, UTEMP] = ...  % Q2 stage
            QuadraticNeighbour(Nx, Nz, Neb1, Neb2, Neb3, Neb4, Umid, UTEMP, U, emid, emid1, ee, ee1, a0, b, t, dt);

        U(:)      = UTEMP(:);
        t         = t + dt;


    end
    %=========================================================
    % Calculate the Errors:

    RunT(ih) = toc;
    disp([' Iteration Number ',' is ',num2str(ih), ...
       ' out of  ', num2str(Mh),'   Running Time', ' is ', num2str(RunT(ih)) ])

    Kul = zeros(N,1); MaxKul = 0;

    for i=1:N
        Kul(i) = U(i) - BenchMarklSolution(i);
        if MaxKul < abs(Kul(i))
            MaxKul=abs(Kul(i));
        end

    end
    MaxD(ih) = MaxKul;


end
%%%% End of Section 5:
%% Section 00: Save the Data and Post-Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%******************************************************************
% Save the Data

save DataOfSolution

%******************************************************************
% Go to "PlotResult.m" file and run the code to process the data


warndlg([' The Numerical calculations have been achieved. ' ...
    'Please Run the PlotResult Script'],'Important')


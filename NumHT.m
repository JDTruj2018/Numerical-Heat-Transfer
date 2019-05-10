%% 18.0851 Project
% Author      : Jered Dominguez-Trujillo
% Date        : May 2, 2019
% Description : Numerical Solution to Heat Equation

% SCHEME = 0 -> EXPLICIT
% SCHEME = 1 -> IMPLICIT
% SCHEME = 2 -> CRANK_NICOLSON

function U = NumHT(SCHEME, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG, W, SIGMA, QMAX)
    clc;
    
    % Boundary Conditions
    LENGTH = L;                             % Length of Domain
    TIME_RAMP = TR;                         % Time Ramp at X = 0 for BC to Go from 0 to C1
    C1 = BC1;                               % Dirilecht Condition U(x) = C1 at X = 0
    C2 = BC2;                               % Neumann Condition at dU/dX = C2 at X = L
    K = KT;                                 % Thermal Diffusivity
    
    % Default Arguments
    if nargin <= 12
        QMAX = 1;                           % Default Maximum Value of Heat Source Function
        
        if nargin <= 11
            SIGMA = 0.2 .* LENGTH;           % Default Standard Deviation of Heat Source Function
        end
        
        if nargin <= 10
            W = 1;                          % Default Frequency of Heat Source Function (rad/s)
        end
    end
    
    % Spatial Domain
    NODES = NX;                             % Nodes
    DX = LENGTH ./ NODES;                    % DX Calculation
    X = linspace(0, LENGTH + DX, NODES + 2);% X Vector with Ghost Node
    
    % Time Domain
    TMAX = TM;                              % End Time of Simulation
    DT = TMAX ./ NT;                         % DT Calculation
    TIMESTEPS = TMAX ./ DT + 1;              % Number of Time Steps
    T = linspace(0, TMAX, TIMESTEPS);       % Time Vector
    
    % Calculate CFL Number     
    CFL = (DT .* K) ./ (DX .* DX);          % Multiplication Factor K (DT / DX^2)
    
    % Print Out Simulation Info
    if SCHEME == 0
        fprintf('Explicit Method:\n');
    elseif SCHEME == 1
        fprintf('Implicit Method:\n');
    elseif SCHEME == 2
        fprintf('Crank-Nicolson Method:\n');
    end
    
    fprintf('Thermal Diffusivity [K]: %.3f\n', K);
    fprintf('BCs:\n(1) U(x) = %.2f At X = 0\n(2) dU/dX = %.2f at X = L = %.2f\n\n', C1, C2, L);
    fprintf('Length: %.2f\t\tDX: %.5f\t\tNodes: %.0f\n', LENGTH, DX, NODES);
    fprintf('Max Time: %.2f\t\tDT: %.5f\t\t\n\n', TMAX, DT);
    if SOURCE_FLAG == 1
        fprintf('Q(x) = - QMAX * SIN(OMEGA * T) * EXP(-(X - L/2)^2 / SIGMA ^ 2)\n');
        fprintf('QMAX: %.2f\t\tOMEGA: %.2f rad/s\t\tSIGMA: %.2f\t\t L: %.2f\n\n', QMAX, W, SIGMA, LENGTH);
    end
    fprintf('CFL Number: %.5f\n', CFL);
    
    % Initialize Matrices
    U = zeros(TIMESTEPS, NODES + 2);
    Q = zeros(TIMESTEPS, NODES + 2);
    
    % Time Tolerance
    eps = 10^-7;
    
    % If Heat Source Flag is True
    if SOURCE_FLAG == 1
        for ii = 1:TIMESTEPS
            Q(ii, :) = - QMAX .* (sin((W .* T(ii)))) .* exp(-((X - (LENGTH ./ 2)) .^ 2) ./ (SIGMA .^ 2));
        end
    end
      
    % Initialize
    CURRENT_T = 0;
    TIMESTEP = 0;
    
    % Plot Initial Conditions
    fTemperature = figure('Name', 'Temperature History', 'NumberTitle', 'off');
    figure(fTemperature);
    
    % Temperature Profile
    subplot(211);
    plot(X, U(1, :), 'o-');
    xlabel('X', 'FontSize', 18); ylabel('Temperature [u]', 'Fontsize', 14);
    axis([0 LENGTH -2 2]);
    
    % Heat Source Profile
    subplot(212);
    plot(X, Q(1, :), 'o-');
    xlabel('X', 'FontSize', 18); ylabel('Heat Source [Q]', 'Fontsize', 14);
    axis([0 LENGTH -QMAX QMAX]);
    
    suptitle(['Time = ', num2str(CURRENT_T), ' s']);
    % pause();
    
    % Explicit Scheme
    if SCHEME == 0
        
        % Iterate Until TMAX
        while CURRENT_T < TMAX - eps
            CURRENT_T = CURRENT_T + DT;
            TIMESTEP = TIMESTEP + 1;

            n = TIMESTEP;

            % Explicit Euler Method
            for ii = 2:NODES + 1
                U(n + 1, ii) = CFL .* U(n, ii + 1) + (1 - 2 .* CFL) .* U(n, ii) + CFL .* U(n, ii - 1) + Q(n, ii) .* DT;
            end

            % Boundary Conditions
            U(n + 1, 1) = C1 .* (CURRENT_T / TIME_RAMP) .* (CURRENT_T <= TIME_RAMP) + C1 .* (CURRENT_T > TIME_RAMP);
            U(n + 1, NODES + 2) = U(n + 1, NODES) + 2 .* DX .* C2;

            % Plot New Temperature Profile
            figure(fTemperature);

            % Temperature Profile
            subplot(211);
            plot(X, U(n + 1, :), 'o-');
            xlabel('X', 'FontSize', 18); ylabel('Temperature [u]', 'Fontsize', 14);
            axis([0 L -2 2]);

            % Heat Source Profile
            subplot(212);
            plot(X, Q(n + 1, :), 'o-');
            xlabel('X', 'FontSize', 18); ylabel('Heat Source [Q]', 'Fontsize', 14);
            axis([0 LENGTH -QMAX QMAX]);

            suptitle(['Time = ', num2str(CURRENT_T), ' s']);

            pause(0.01);
        end
    
    % Implicit Scheme
    elseif SCHEME == 1
        LOWER = zeros(1, NODES);
        DIAG  = zeros(1, NODES + 1);
        UPPER = zeros(1, NODES);
        
        LOWERI = 2:NODES + 1; LOWERJ = 1:NODES;
        DIAGI = 1:NODES + 1; DIAGJ = 1:NODES + 1;
        UPPERI = 1:NODES; UPPERJ = 2:NODES + 1;
        
        DIAG(1) = 1; UPPER(1) = 0;
        
        for ii = 2:NODES
            LOWER(ii) = -CFL;
            DIAG(ii) = 2 .* CFL + 1;
            UPPER(ii) = -CFL;
        end
        
        DIAG(NODES + 1) = 2 .* CFL + 1;
        
        % Iterate Until TMAX
        while CURRENT_T < TMAX - eps
            CURRENT_T = CURRENT_T + DT;
            TIMESTEP = TIMESTEP + 1;

            n = TIMESTEP;
            
            RHS = zeros(NODES + 1, 1);
            
            % Implicit Euler Method
            for ii = 2:NODES + 1
                RHS(ii) = U(n, ii) + Q(n, ii) .* DT;
            end
            
            % Boundary Conditions
            fBC1 = C1 .* (CURRENT_T / TIME_RAMP) .* (CURRENT_T <= TIME_RAMP) + C1 .* (CURRENT_T > TIME_RAMP);
            LOWER(1) = 0; 
            RHS(1) = fBC1;
            RHS(2) = RHS(2) + CFL .* fBC1;
            
            % Forward Difference Scheme
            % DIAG(NODES + 1) = CFL + 1;
            % RHS(NODES + 1) = RHS(NODES + 1) + C2 .* CFL .* DX;
            
            % Central Difference Scheme
            L(NODES) = -2 .* CFL;
            RHS(NODES + 1) = RHS(NODES + 1) + 2 .* C2 .* CFL .* DX;
            
            MA = sparse([LOWERI, DIAGI, UPPERI], [LOWERJ, DIAGJ, UPPERJ], [LOWER, DIAG, UPPER], NODES + 1, NODES + 1);
            
            U(n + 1, 1:end-1) = MA \ RHS;
            
            % Forward Difference Scheme
            % U(n + 1, end) = U(n + 1, end - 1) + C2 .* DX;
            
            % Central Difference Scheme
            U(n + 1, end) = U(n + 1, end - 2) + 2 .* C2 .* DX;
            
            % Plot New Temperature Profile
            figure(fTemperature);

            % Temperature Profile
            subplot(211);
            plot(X, U(n + 1, :), 'o-');
            xlabel('X', 'FontSize', 18); ylabel('Temperature [u]', 'Fontsize', 14);
            axis([0 LENGTH -2 2]);

            % Heat Source Profile
            subplot(212);
            plot(X, Q(n + 1, :), 'o-');
            xlabel('X', 'FontSize', 18); ylabel('Heat Source [Q]', 'Fontsize', 14);
            axis([0 LENGTH -QMAX QMAX]);

            suptitle(['Time = ', num2str(CURRENT_T), ' s']);
            
            pause(0.01);
        end
        
    % Crank-Nicolson Scheme
    elseif SCHEME == 2
        LOWER = zeros(1, NODES);
        DIAG  = zeros(1, NODES + 1);
        UPPER = zeros(1, NODES);
        
        LOWERI = 2:NODES + 1; LOWERJ = 1:NODES;
        DIAGI = 1:NODES + 1; DIAGJ = 1:NODES + 1;
        UPPERI = 1:NODES; UPPERJ = 2:NODES + 1;
        
        DIAG(1) = 1; UPPER(1) = 0;
        
        for ii = 2:NODES
            LOWER(ii) = -CFL ./ 2;
            DIAG(ii) = CFL + 1;
            UPPER(ii) = -CFL ./ 2;
        end
        
        DIAG(NODES + 1) = CFL + 1;
        
        % Iterate Until TMAX
        while CURRENT_T < TMAX - eps
            CURRENT_T = CURRENT_T + DT;
            TIMESTEP = TIMESTEP + 1;

            n = TIMESTEP;
            
            RHS = zeros(NODES + 1, 1);
            
            % Crank-Nicolson Method
            for ii = 2:NODES + 1
                RHS(ii) = CFL .* U(n, ii - 1) ./ 2 + (1 - CFL) .* U(n, ii) + CFL .* U(n, ii + 1) ./ 2 + Q(n, ii) .* DT;
            end
            
            % Boundary Conditions
            fBC1 = C1 .* (CURRENT_T / TIME_RAMP) .* (CURRENT_T <= TIME_RAMP) + C1 .* (CURRENT_T > TIME_RAMP);
            LOWER(1) = 0; 
            RHS(1) = fBC1;
            RHS(2) = RHS(2) + CFL .* fBC1 ./ 2;
            
            % Forward Difference Scheme
            % DIAG(NODES + 1) = CFL ./ 2 + 1;
            % RHS(NODES + 1) = RHS(NODES + 1) + C2 .* CFL .* DX ./ 2;
            
            % Central Difference Scheme
            L(NODES) = -CFL;
            RHS(NODES + 1) = RHS(NODES + 1) + C2 .* CFL .* DX;
            
            MA = sparse([LOWERI, DIAGI, UPPERI], [LOWERJ, DIAGJ, UPPERJ], [LOWER, DIAG, UPPER], NODES + 1, NODES + 1);
            
            U(n + 1, 1:end-1) = MA \ RHS;
            
            % Forward Difference Scheme
            % U(n + 1, end) = U(n + 1, end - 1) + C2 .* DX;
            
            % Central Difference Scheme
            U(n + 1, end) = U(n + 1, end - 2) + 2 .* C2 .* DX;
            
            % Plot New Temperature Profile
            figure(fTemperature);

            % Temperature Profile
            subplot(211);
            plot(X, U(n + 1, :), 'o-');
            xlabel('X', 'FontSize', 14); ylabel('Temperature [u]', 'Fontsize', 14);
            axis([0 LENGTH -2 2]);

            % Heat Source Profile
            subplot(212);
            plot(X, Q(n + 1, :), 'o-');
            xlabel('X', 'FontSize', 14); ylabel('Heat Source [Q]', 'Fontsize', 14);
            axis([0 LENGTH -QMAX QMAX]);

            suptitle(['Time = ', num2str(CURRENT_T), ' s']);
            
            pause(0.01);
        end
    end
    
    close(fTemperature);
end


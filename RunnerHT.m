%% 18.0851 Project
% Author      : Jered Dominguez-Trujillo
% Date        : May 9, 2019
% Description : Wrapper for NumHT.m

% SCHEME = 0 -> EXPLICIT
% SCHEME = 1 -> IMPLICIT
% SCHEME = 2 -> CRANK_NICOLSON

clear all; close all;
COLORS = get(gca,'colororder');

%% Plot 3 Methods at t = tfinal against each other without Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; MNX = 2^18; TM = 40; MNT = 2^12; TR = 1; SOURCE_FLAG = 0;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, MNX, TM, MNT, TR, SOURCE_FLAG);

% 3 Method Runs
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; NX = 2^6; TM = 40; NT = 2^10; TR = 1; SOURCE_FLAG = 0;
DXX = L ./ NX; XX = linspace(0, L + DXX, NX + 2);
DT = TM ./ NT; TIMESTEPS = TM ./ DT + 1; TT = linspace(0, TM, TIMESTEPS);       
    
AA = NumHT(0, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
BB = NumHT(1, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
CC = NumHT(2, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);

BASELINE = BASELINE(1:MNT/NT:end, 1:MNX/NX:end);
ResA = BASELINE - AA(:, 1:end-1); ErrorA = sqrt(sum(ResA' .* ResA') ./ NX);
ResB = BASELINE - BB(:, 1:end-1); ErrorB = sqrt(sum(ResB' .* ResB') ./ NX);
ResC = BASELINE - CC(:, 1:end-1); ErrorC = sqrt(sum(ResC' .* ResC') ./ NX);

fErrorTime1 = figure('Name', 'Error Evolution with Time', 'NumberTitle', 'off');
figure(fErrorTime1); hold on;

plot(TT, ErrorA, '-', 'LineWidth', 2, 'DisplayName', 'Explicit Euler');
plot(TT, ErrorB, '-', 'LineWidth', 2, 'DisplayName', 'Implicit Euler');
plot(TT, ErrorC, '-', 'Linewidth', 2, 'DisplayName', 'Crank-Nicolson');

xlabel('Time', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
title('Error Evolution', 'FontSize', 18); legend('show');

saveas(fErrorTime1, 'Figures/MATLAB/NoSourceErrorTime.png');
saveas(fErrorTime1, 'Figures/MATLAB/Figs/NoSourceErrorTime.fig');

fCompare1 = figure('Name', 'Solution Comparison: T = 40', 'NumberTitle', 'off');
figure(fCompare1); hold on;

subplot(2, 2, 1);
plot(XX(1:end-1), BASELINE(end, :), '-', 'Color', COLORS(4, :), 'LineWidth', 2, 'DisplayName', 'Baseline');
xlabel('X'); ylabel('Temperature [u]');
title('Baseline (Analytical)'); axis([0 L -2 2]);

subplot(2, 2, 2);
plot(XX, AA(end, :), '-', 'Color', COLORS(1, :), 'LineWidth', 2, 'DisplayName', 'Explicit Euler');
xlabel('X'); ylabel('Temperature [u]');
title('Explicit Euler'); axis([0 L -2 2]);



subplot(2, 2, 3);
plot(XX, BB(end, :), '-', 'Color', COLORS(2, :), 'LineWidth', 2, 'DisplayName', 'Implicit Euler');
xlabel('X'); ylabel('Temperature [u]');
title('Implicit Euler'); axis([0 L -2 2]);

subplot(2, 2, 4);
plot(XX, CC(end, :), '-', 'Color', COLORS(3, :), 'LineWidth', 2, 'DisplayName', 'Crank-Nicolson');
xlabel('X'); ylabel('Temperature [u]');
title('Crank-Nicolson'); axis([0 L -2 2]);

suptitle('Time = 40 s');  

saveas(fCompare1, 'Figures/MATLAB/NoSourceCompare.png');
saveas(fCompare1, 'Figures/MATLAB/Figs/NoSourceCompare.fig');

%% Plot 3 Methods at t = tfinal against each other with Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; MNX = 2^18; TM = 40; MNT = 2^12; TR = 1; SOURCE_FLAG = 1;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, MNX, TM, MNT, TR, SOURCE_FLAG);

% 3 Method Runs
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; NX = 2^6; TM = 40; NT = 2^10; TR = 1; SOURCE_FLAG = 0;

DXX = L ./ NX; XX = linspace(0, L + DXX, NX + 2);

DT = TM ./ NT; TIMESTEPS = TM ./ DT + 1; TT = linspace(0, TM, TIMESTEPS); 

AA = NumHT(0, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
BB = NumHT(1, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
CC = NumHT(2, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);

BASELINE = BASELINE(1:MNT/NT:end, 1:MNX/NX:end);
ResA = BASELINE - AA(:, 1:end-1); ErrorA = sqrt(sum(ResA' .* ResA') ./ NX);
ResB = BASELINE - BB(:, 1:end-1); ErrorB = sqrt(sum(ResB' .* ResB') ./ NX);
ResC = BASELINE - CC(:, 1:end-1); ErrorC = sqrt(sum(ResC' .* ResC') ./ NX);

fErrorTime2 = figure('Name', 'Error Evolution with Time', 'NumberTitle', 'off');
figure(fErrorTime2); hold on;

plot(TT, ErrorA, '-', 'LineWidth', 2, 'DisplayName', 'Explicit Euler');
plot(TT, ErrorB, '-', 'LineWidth', 2, 'DisplayName', 'Implicit Euler');
plot(TT, ErrorC, '-', 'Linewidth', 2, 'DisplayName', 'Crank-Nicolson');

xlabel('Time', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
title('Error Evolution', 'FontSize', 18); legend('show');

saveas(fErrorTime2, 'Figures/MATLAB/SourceErrorTime.png');
saveas(fErrorTime2, 'Figures/MATLAB/Figs/SourceErrorTime.fig');

BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; NX = 2^6; TM = 40; NT = 2^10; TR = 1; SOURCE_FLAG = 1;

DXX = L ./ NX; XX = linspace(0, L + DXX, NX + 2);

AA = NumHT(0, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
BB = NumHT(1, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
CC = NumHT(2, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);





fCompare2 = figure('Name', 'Solution Comparison: T = 40', 'NumberTitle', 'off');
figure(fCompare2); hold on;

subplot(2, 2, 1);
plot(XX(1:end-1), BASELINE(end, :), '-', 'Color', COLORS(4, :), 'LineWidth', 2, 'DisplayName', 'Baseline');
xlabel('X'); ylabel('Temperature [u]');
title('Baseline (Analytical)'); axis([0 L -2 2]);

subplot(2, 2, 2);
plot(XX, AA(end, :), '-', 'Color', COLORS(1, :), 'LineWidth', 2, 'DisplayName', 'Explicit Euler');
xlabel('X'); ylabel('Temperature [u]');
title('Explicit Euler'); axis([0 L -2 2]);

subplot(2, 2, 3);
plot(XX, BB(end, :), '-', 'Color', COLORS(2, :), 'LineWidth', 2, 'DisplayName', 'Implicit Euler');
xlabel('X'); ylabel('Temperature [u]');
title('Implicit Euler'); axis([0 L -2 2]);

subplot(2, 2, 4);
plot(XX, CC(end, :), '-', 'Color', COLORS(3, :), 'LineWidth', 2, 'DisplayName', 'Crank-Nicolson');
xlabel('X'); ylabel('Temperature [u]');
title('Crank-Nicolson'); axis([0 L -2 2]);

suptitle('Time = 40 s');  

saveas(fCompare2, 'Figures/MATLAB/SourceCompare.png');
saveas(fCompare2, 'Figures/MATLAB/Figs/SourceCompare.fig');


































%% Hold DX Constant with Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi;NX = 2^6; TM = 40; NT = 2^14; TR = 1; SOURCE_FLAG = 1;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
BLAST = BASELINE(end, :);

sch = {'Explicit Euler', 'Implicit Euler', 'Crank-Nicolson'};
    
% Run Through DT
MAXNT = 2^11; MINNT = 2^3;

NT = round(logspace(log(MINNT)/log(10), log(MAXNT)/log(10), 20), 0); DT = TM ./ NT;
RVX = zeros(3, length(NX));

fErrorDT1 = figure('Name', 'Numerical Error - Constant DX', 'NumberTitle', 'off');
figure(fErrorDT1);

for SCHEME = 0:2
    for ii = 1:length(NT)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX, TM, NT(ii), TR, SOURCE_FLAG);
        U = U(:, 1:end-1);
        ULAST = U(end, :);
        
        B = BLAST(1:end-1);
        
        RES = ULAST - B; 
        RVX(SCHEME + 1, ii) = sqrt(sum(RES .* RES) ./ NX);
        
        clf;
        hold off;
        for jj = 0:SCHEME
            if jj == 0
                inx = find(RVX(jj + 1, :) < 10^2, 1);
                loglog(DT(inx:ii), RVX(jj + 1, inx:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            else
                loglog(DT(1:ii), RVX(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            end
            hold on;
        end
        
        loglog(DT, DT, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order');
        loglog(DT, DT .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order');

        xlabel('DT', 'Fontsize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Time Step', 'FontSize', 18);
        ylim([0, max(DT) .^ 2]); legend('show');
    end
end

saveas(fErrorDT1, 'Figures/MATLAB/ConstantDXSource.png');
saveas(fErrorDT1, 'Figures/MATLAB/Figs/ConstantDXSource.fig');







%% Hold DT Constant with Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; MNX = 2^18; TM = 40; NT = 2^12; TR = 1; SOURCE_FLAG = 1;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, MNX, TM, NT, TR, SOURCE_FLAG);
BLAST = BASELINE(end, :);

sch = {'Explicit Euler', 'Implicit Euler', 'Crank-Nicolson'};
    
% Run through DX
MAXNX = 2^8; MINNX = 2^3;

NX = round(logspace(log(MINNX)/log(10), log(MAXNX)/log(10), 20), 0); DX = L ./ NX;

RVT = zeros(3, length(NX));
fErrorDX1 = figure('Name', 'Numerical Error - Constant DT', 'NumberTitle', 'off');
figure(fErrorDX1);

for SCHEME = 0:2
    for ii = 1:length(NX)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX(ii), TM, NT, TR, SOURCE_FLAG);
        U = U(:, 1:end-1);
        ULAST = U(end, :);
        
        B = BLAST(1:MNX/NX(ii):end-1);
        
        RES = ULAST - B; 
        RVT(SCHEME + 1, ii) = sqrt(sum(RES .* RES) ./ NX(ii));
        
        clf;
        hold off;
        for jj = 0:SCHEME
            if jj == 0
                inx = find(RVT(jj + 1, :) < 10^2, 1, 'last');
                loglog(DX(1:inx), RVT(jj + 1, 1:inx), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            else
                loglog(DX(1:ii), RVT(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            end
            hold on;
        end
        
        loglog(DX, DX, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order');
        loglog(DX, DX .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order');

        xlabel('DX', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Spatial Step', 'FontSize', 18);
        ylim([0, max(DX) .^ 2]); legend('show'); 
    end
end

saveas(fErrorDX1, 'Figures/MATLAB/ConstantDTSource.png');
saveas(fErrorDX1, 'Figures/MATLAB/Figs/ConstantDTSource.fig');







%% Hold DX Constant without Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; NX = 2^6; TM = 40; NT = 2^14; TR = 1; SOURCE_FLAG = 0;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, NX, TM, NT, TR, SOURCE_FLAG);
BLAST = BASELINE(end, :);

sch = {'Explicit Euler', 'Implicit Euler', 'Crank-Nicolson'};
    
% Run Through DT
MAXNT = 2^11; MINNT = 2^3;

NT = round(logspace(log(MINNT)/log(10), log(MAXNT)/log(10), 20), 0);
DT = TM ./ NT;

RVX = zeros(3, length(NX));

fErrorDT2 = figure('Name', 'Numerical Error - Constant DX', 'NumberTitle', 'off');
figure(fErrorDT2);

for SCHEME = 0:2
    for ii = 1:length(NT)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX, TM, NT(ii), TR, SOURCE_FLAG);
        U = U(:, 1:end-1);
        ULAST = U(end, :);
        
        B = BLAST(1:end-1);
        
        RES = ULAST - B; 
        RVX(SCHEME + 1, ii) = sqrt(sum(RES .* RES) ./ NX);
        
        clf;
        hold off;
        for jj = 0:SCHEME
            if jj == 0
                inx = find(RVX(jj + 1, :) < 10^2, 1);
                loglog(DT(inx:ii), RVX(jj + 1, inx:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            else
                loglog(DT(1:ii), RVX(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            end
            hold on;
        end
        
        loglog(DT, DT, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order');
        loglog(DT, DT .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order');

        xlabel('DT', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Time Step', 'FontSize', 18);
        ylim([0, max(DT) .^ 2]); legend('show');
    end
end

saveas(fErrorDT2, 'Figures/MATLAB/ConstantDXNoSource.png');
saveas(fErrorDT2, 'Figures/MATLAB/Figs/ConstantDXNoSource.fig');





%% Hold DT Constant without Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; MNX = 2^18; TM = 40; NT = 2^12; TR = 1; SOURCE_FLAG = 0;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, MNX, TM, NT, TR, SOURCE_FLAG);
BLAST = BASELINE(end, :);

sch = {'Explicit Euler', 'Implicit Euler', 'Crank-Nicolson'};
    
% Run through DX
MAXNX = 2^8; MINNX = 2^3;

NX = round(logspace(log(MINNX)/log(10), log(MAXNX)/log(10), 20), 0); DX = L ./ NX;

RVT = zeros(3, length(NX));
fErrorDX2 = figure('Name', 'Numerical Error - Constant DT', 'NumberTitle', 'off');
figure(fErrorDX2);

for SCHEME = 0:2
    for ii = 1:length(NX)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX(ii), TM, NT, TR, SOURCE_FLAG);
        U = U(:, 1:end-1);
        ULAST = U(end, :);
        
        B = BLAST(1:MNX/NX(ii):end-1);
        
        RES = ULAST - B; 
        RVT(SCHEME + 1, ii) = sqrt(sum(RES .* RES) ./ NX(ii));
        
        clf;
        hold off;
        for jj = 0:SCHEME
            if jj == 0
                inx = find(RVT(jj + 1, :) < 10^2, 1, 'last');
                loglog(DX(1:inx), RVT(jj + 1, 1:inx), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            else
                loglog(DX(1:ii), RVT(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            end
            hold on;
        end
        
        loglog(DX, DX, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order');
        loglog(DX, DX .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order');

        xlabel('DX', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Spatial Step', 'FontSize', 18);
        ylim([0, max(DX) .^ 2]); legend('show'); 
    end
end

saveas(fErrorDX2, 'Figures/MATLAB/ConstantDTNoSource.png');
saveas(fErrorDX2, 'Figures/MATLAB/Figs/ConstantDTNoSource.fig');







%% Hold CFL Constant with Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; MNX = 2^16; TM = 40; NT = 2^14; TR = 1; SOURCE_FLAG = 1;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, MNX, TM, NT, TR, SOURCE_FLAG);
BLAST = BASELINE(end, :);

sch = {'Explicit Euler', 'Implicit Euler', 'Crank-Nicolson'};
    
% Run through DX
MAXNT = 2^11; MINNT = 2^3;
NT = round(logspace(log(MINNT)/log(10), log(MAXNT)/log(10), 9), 0); DT = TM ./ NT;

CFL_HOLD = 0.4;
NX = round(L ./ sqrt(KT .* DT ./ CFL_HOLD), 0); DX = L ./ NX;

RVCFL = zeros(3, length(NX));
fErrorCFL1 = figure('Name', 'Numerical Error - Constant CFL', 'NumberTitle', 'off');
figure(fErrorCFL1);

for SCHEME = 0:2
    for ii = 1:length(NX)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX(ii), TM, NT(ii), TR, SOURCE_FLAG);
        U = U(:, 1:end-1); ULAST = U(end, :); B = BLAST(1:MNX/NX(ii):end-1);
        RES = ULAST - B; RVCFL(SCHEME + 1, ii) = sqrt(sum(RES .* RES) ./ NX(ii));
        
        clf; hold off;
        subplot(2, 1, 1); hold off;
        for jj = 0:SCHEME
            loglog(DX(1:ii), RVCFL(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            hold on;
        end
        loglog(DX, DX, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order Space');
        loglog(DX, DX .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order Space');
        xlabel('DX', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Spatial Step', 'FontSize', 18);
        ylim([0, max(DX) .^ 2]); legend('show'); 
        
        subplot(2, 1, 2); hold off;
        for jj = 0:SCHEME
            loglog(DT(1:ii), RVCFL(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            hold on;
        end
        loglog(DT, DT, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order Time');
        loglog(DT, DT .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order Time');
        xlabel('DT', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Temporal Step', 'FontSize', 18);
        ylim([0, max(DT) .^ 2]); legend('show'); 
    end
end

saveas(fErrorCFL1, 'Figures/MATLAB/ConstantCFLSource.png');
saveas(fErrorCFL1, 'Figures/MATLAB/Figs/ConstantCFLSource.fig');



%% Hold CFL Constant without Source
% Baseline
SCHEME = 2;                 % Crank-Nicolson
BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi; MNX = 2^16; TM = 40; NT = 2^14; TR = 1; SOURCE_FLAG = 0;
BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, MNX, TM, NT, TR, SOURCE_FLAG);
BLAST = BASELINE(end, :);

sch = {'Explicit Euler', 'Implicit Euler', 'Crank-Nicolson'};
    
% Run through DX
MAXNT = 2^11; MINNT = 2^3;
NT = round(logspace(log(MINNT)/log(10), log(MAXNT)/log(10), 9), 0); DT = TM ./ NT;

CFL_HOLD = 0.4;
NX = round(L ./ sqrt(KT .* DT ./ CFL_HOLD), 0); DX = L ./ NX;

RVCFL = zeros(3, length(NX));
fErrorCFL2 = figure('Name', 'Numerical Error - Constant CFL', 'NumberTitle', 'off');
figure(fErrorCFL2);

for SCHEME = 0:2
    for ii = 1:length(NX)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX(ii), TM, NT(ii), TR, SOURCE_FLAG);
        U = U(:, 1:end-1); ULAST = U(end, :); B = BLAST(1:MNX/NX(ii):end-1);
        RES = ULAST - B; RVCFL(SCHEME + 1, ii) = sqrt(sum(RES .* RES) ./ NX(ii));
        
        clf; hold off;
        subplot(2, 1, 1); hold off;
        for jj = 0:SCHEME
            loglog(DX(1:ii), RVCFL(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            hold on;
        end
        loglog(DX, DX, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order Space');
        loglog(DX, DX .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order Space');
        xlabel('DX', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Spatial Step', 'FontSize', 18);
        ylim([0, max(DX) .^ 2]); legend('show'); 
        
        subplot(2, 1, 2); hold off;
        for jj = 0:SCHEME
            loglog(DT(1:ii), RVCFL(jj + 1, 1:ii), '-o', 'Color', COLORS(jj + 1, :), 'LineWidth', 2, 'DisplayName', sch{jj+1});
            hold on;
        end
        loglog(DT, DT, '--', 'Color', COLORS(5, :), 'LineWidth', 2, 'DisplayName', 'First-Order Time');
        loglog(DT, DT .^ 2, '--', 'Color', COLORS(6, :), 'LineWidth', 2, 'DisplayName', 'Second-Order Time');
        xlabel('DT', 'FontSize', 14); ylabel('Error', 'FontSize', 14);
        title('Error as a Function of Temporal Step', 'FontSize', 18);
        ylim([0, max(DT) .^ 2]); legend('show'); 
    end
end

saveas(fErrorCFL2, 'Figures/MATLAB/ConstantCFLNoSource.png');
saveas(fErrorCFL2, 'Figures/MATLAB/Figs/ConstantCFLNoSource.fig');
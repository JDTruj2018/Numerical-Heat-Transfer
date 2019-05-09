%% 18.0851 Project
% Author      : Jered Dominguez-Trujillo
% Date        : May 9, 2019
% Description : Wrapper for NumHT.m

% SCHEME = 0 -> EXPLICIT
% SCHEME = 1 -> IMPLICIT
% SCHEME = 2 -> CRANK_NICOLSON

clear all; close all;

%% Baseline
SCHEME = 2;                 % Crank-Nicolson

BC1 = 1; BC2 = -0.2; KT = 0.1; L = 2*pi;
MNX = 2^16; TM = 40; NT = 50; TR = 1; SOURCE_FLAG = 0;

BASELINE = NumHT(SCHEME, BC1, BC2, KT, L, MNX, TM, NT, TR, SOURCE_FLAG);
BLAST = BASELINE(end, :);

%% Hold DX Constant
MAXNT = 2^12; MINNT = 2^8;

NT = round(logspace(log(MINNT)/log(10), log(MAXNT)/log(10), 9), 0);
NX = 2^6;
DT = TM ./ NT;

RVX = zeros(3, length(NX));
    
for SCHEME = 0:2

    for ii = 1:length(NT)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX, TM, NT(ii), TR, SOURCE_FLAG);
        U = U(:, 1:end-1);
        ULAST = U(end, :);
        
        B = BLAST(1:MNX/NX:end-1);
        
        RES = ULAST - B; 
        RVX(SCHEME + 1, ii) = sqrt(sum(RES .* RES));
    end
end

fErrorDT = figure('Name', 'Numerical Error - Constant DX', 'NumberTitle', 'off');
figure(fErrorDT);
hold on;

for SCHEME = 0:2
    if SCHEME == 0
        sch = 'Explicit Euler';
    elseif SCHEME == 1
        sch = 'Implicit Euler';
    elseif SCHEME == 2
        sch = 'Crank-Nicolson';
    end
    
    loglog(DT, RVX(SCHEME + 1, :), 'LineWidth', 2, 'DisplayName', sch);
end

xlabel('DT'); ylabel('Error');
ylim([0, 1]);
legend('show');

%% Hold DT Constant
MAXNX = 2^10; MINNX = 2^3;

NX = round(logspace(log(MINNX)/log(10), log(MAXNX)/log(10), 9), 0);
NT = 2^10;
DX = L ./ NX;

RVT = zeros(3, length(NX));
    
for SCHEME = 0:2

    for ii = 1:length(NX)
        U = NumHT(SCHEME, BC1, BC2, KT, L, NX(ii), TM, NT, TR, SOURCE_FLAG);
        U = U(:, 1:end-1);
        ULAST = U(end, :);
        
        B = BLAST(1:MNX/NX(ii):end-1);
        
        RES = ULAST - B; 
        RVT(SCHEME + 1, ii) = sqrt(sum(RES .* RES));
    end
end

fErrorDX = figure('Name', 'Numerical Error - Constant DT', 'NumberTitle', 'off');
figure(fErrorDX);
hold on;

for SCHEME = 0:2
    if SCHEME == 0
        sch = 'Explicit Euler';
    elseif SCHEME == 1
        sch = 'Implicit Euler';
    elseif SCHEME == 2
        sch = 'Crank-Nicolson';
    end
    
    loglog(DX, RVT(SCHEME + 1, :), 'LineWidth', 2, 'DisplayName', sch);
end

xlabel('DX'); ylabel('Error');
ylim([0, 1]);
legend('show');
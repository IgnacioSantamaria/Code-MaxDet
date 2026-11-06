% This script shows the performance of the Max-Det BD-RIS solution
% described in [A]. It allows to reproduce the results of Fig. 5 of [A], showing 
% the achievable rate vs SNR for a given MxM BD-RIS. 
%
% We compare the following methods:
% a) Closed-form Max-Det symmetric and unitary solution (method in [A])
% b) Unitary but not symmetric BD-RIS maximizing capacity
% c) Manifold optimization (MO) algorithm (this is the method in [B]): Obtains a full-rank
% unitary and symmetric BD-RIS maximizing capacity.
% d) Diagonal RIS maximizing capacity.
% 
% [A] I. Santamaria, M. Soleymani, J. Gutierrez, E. Jorswieck, "On the optimal 
% symmetric BD-RIS maximizing the determinant of a MIMO link", IEEE
% Transactions on Signal Processing, 2026.
%
% [B] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutiérrez and C. Beltrán, 
% "Riemannian Optimization on the Manifold of Unitary and Symmetric
% Matrices with Application to BD-RIS-Assisted Systems",
% IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 
% Barcelona, Spain, 2026, pp. 531–535
%
% I. Santamaria, UC Nov. 2026

format compact
clc; clear;
%% Parameters
Nt = 4;                 % Number of transmit antennas
Nr = 4;                 % Number of receive antennas
r = min(Nt,Nr);         % DoF
M = [16,64,256];        % Number of BD-RIS elements  
SNR = -10:10;                 % Received SNR in dBs
PmaxdBm = 0;                  % Pmax (in dBm)
Pmax = 10.^(PmaxdBm/10);      % Pmax
Qiso = (Pmax/Nt)*eye(Nt);     % Isotropic covariance matrix (fixed)
NsimMC = 200;                 % Number of Monte Carlo simulations (you should use at least 25 to get a reasonable smooth curve)

%% MO algorithm parameters (ICASSP'26 algorithm)
% Parameters for the method in [B]

opt_paramsBDRIS_MO = struct();
opt_paramsBDRIS_MO.maxiter = 100;        % Maximum number of iterations
opt_paramsBDRIS_MO.threshold = 1e-3;     % To check convergence

%% Parameters optimization RIS
% It implements the algorithm in S. Zhang, R. Zhang, "Capacity
% characterization for IRS aided MIMO Communications", JSAC, vol. 38, no.
% 8, pp. 1823-1838, 2020.
opt_paramsRIS = struct();
opt_paramsRIS.niterRIS = 100;
opt_paramsRIS.thresholdRIS = 1e-3;

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 6;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Channel parameters
channelparams = struct;
channelparams.blocked = 1;         % Set to 1 if direct channel is blocked (Note: The Max-Det solution assumes a blocked direct link!)
channelparams.RiceRIS = 3;         % Rician factor for the channel btw RIS and Tx/Rx (if Inf-> pure LoS channels)
channelparams.RiceDirect = 0;      % Rician factor for the direct link (if 0 -> Rayleigh fading)
channelparams.pl_0 = -28;          % Path loss at a reference distance (d_0)
channelparams.alpha_RIS = 2;       % Path loss exponent for the RIS links
channelparams.alpha_direct = 4;    % Path loss exponent for the direct links
channelparams.ray_fading = 0;      % Set to 1 if all channels Rayleigh

% ===== Position of the Tx/Rx/RIS (units in meters) ======
sqr_size = 50;                  % square of size sqr_size
PosTx_XYZ = [0 0 1.5];          % Position Tx
PosRx_XYZ = [sqr_size 0 1.5];   % Position Rx
PosRIS_XYZ = [5, 3, 3];         % (x,y,z) coordinates for the RIS position

%% Variables to store the rates (b/s/Hz)
RateMaxDet = zeros(length(M),length(SNR));    % BD-RIS (closed-form Max-Det)
RateMO = zeros(length(M),length(SNR));        % BD-RIS (unitary+ symmetric MO)
RateUnit = zeros(length(M),length(SNR));      % BD-RIS (Unitary)
RateRIS = zeros(length(M),length(SNR));       % BD-RIS (Unitary)

for mm = 1:NsimMC  
    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(M)
        disp(['RIS elements:', num2str(M(dd))])

        %% Generate channels        
        [Hd,G,F] = ChannelsMIMO(M(dd),Nr,Nt,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  
       
        %% Max-Det BDRIS (Solution in [A])
        BDRISCF = BDRIS_MaxDet(F,G,channelparams.blocked);
        Heq_CF = Hd + F*BDRISCF*G';
        for hh = 1:length(SNR)
            sigma2n = Pmax/(Nt*Nr)*real(norm(F*G','fro')^2)*10.^(-SNR(hh)/10);  % noise variance
            RateMaxDet(dd,hh) =  RateMaxDet(dd,hh) + log2(real(det(eye(Nr)+ (Heq_CF*Qiso*Heq_CF')./sigma2n)));
        end
        %% Unitary solution
        [UF,SigmaF,VF] = svd(F,'econ');
        [UG,SigmaG,VG] = svd(G,'econ');
        BDRISU = VF(:,1:r)*VG(:,1:r)';  % we consider a low-rank version even though the full-rank version gets the same results
        HeqU = Hd + F*BDRISU*G';
        for hh =1:length(SNR)
            sigma2n = Pmax/(Nt*Nr)*real(norm(F*G','fro')^2)*10.^(-SNR(hh)/10); % noise variance
            RateUnit(dd,hh) =  RateUnit(dd,hh) + log2(real(det(eye(Nr)+ (HeqU*Qiso*HeqU')./sigma2n)));
        end
        %% Max-Cap BDRIS (MO algorithm)
        % Initialization (unitary+symmetric matrix)
        [Uf,Df,Vf] = svd(F);
        [Ug,Dg,Vg] = svd(G');
        Q = Vf*Ug';
        Thetaini = (Q + Q.')/2;       % symmetric projection
        [F_SVD,~,G_SVD] = svd(Thetaini);
        dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);   
        Thetaini = F_SVD*F_SVD.';
        for hh =1:length(SNR)
            sigma2n = Pmax/(Nt*Nr)*real(norm(F*G','fro')^2)*10.^(-SNR(hh)/10); % noise variance
            [CfinalMO, CtotalMO, BDRIS_MO] = OptimizeBDRIS_MOUs_FullRank(Hd,F,G,Thetaini,Qiso,sigma2n,opt_paramsBDRIS_MO);
            Heq_MO = Hd + F*BDRIS_MO*G';
            RateMO(dd,hh) =  RateMO(dd,hh) + log2(real(det(eye(Nr)+ (Heq_MO*Qiso*Heq_MO')./sigma2n)));
        end

        %% RIS
        for hh =1:length(SNR)
            sigma2n = Pmax/(Nt*Nr)*real(norm(F*G','fro')^2)*10.^(-SNR(hh)/10); % noise variance
            RIS = diag(exp(1i*2*pi*rand(M(dd),1)));  % initial RIS
            [CfinalRIS, CtotalRIS, RIS] = OptimizeRIS(Hd,F,G,RIS,Qiso,sigma2n,opt_paramsRIS);
            Heq_RIS = Hd + F*RIS*G';
            RateRIS(dd,hh) =  RateRIS(dd,hh) + log2(real(det(eye(Nr)+ (Heq_RIS*Qiso*Heq_RIS')./sigma2n)));
        end
    end
end

RateMaxDet = RateMaxDet/NsimMC;
RateMO = RateMO/NsimMC;
RateUnit = RateUnit/NsimMC;
RateRIS = RateRIS/NsimMC;

%% Plot results
figure(30);clf; plot(SNR,RateUnit(1,:), 'k--','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(SNR,RateMaxDet(1,:),'b-s','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMO(1,:),'r-.','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateRIS(1,:),'g-+','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateUnit(3,:), 'k--','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMaxDet(3,:),'b-s','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMO(3,:),'r-.','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateRIS(3,:),'g-+','MarkerSize',ms,'LineWidth',lw);

legend('Max-Rate Unit.','Max-Det Unit.+Symm', 'Max-Rate Unit.+Symm. (MO)',...
    'Max-Rate RIS', 'Location','best');
ylabel('Rate (b/s/Hz)');
xlabel('SNR (dB)');
grid on;
hold off


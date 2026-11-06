% This script allows to reproduce the results of Fig. 4 of [A], showing 
% the achievable rate vs SNR for a given MxM BD-RIS. This script considers the otpimization of the Tx 
% covariance matrix through an alternating optimization algorithm. 
% Only one step of the AO is enough to get most of the improvement, so we
% only optimize the Tx cov. matrix after optimizing the BD-RIS
%
% We compare the following methods:
% a) Closed-form Max-Det solution (method in [A]) + Watefilling
% b) Closed-form Max-Det solution (method in [A]) + Isotropic covariance
% c) Manifold optimization (MO) algorithm (method in [B]) + Waterfilling
% d) Manifold optimization (MO) algorithm (method in [B]) + Isotropic
% 
% 
% [A] I. Santamaria, M. Soleymani, J. Gutierrez, E. Jorswieck, "On the optimal 
% symmetric BD-RIS maximizing the determinant of a MIMO link", IEEE
% Transactions on Signal Processing, 2026.
%
% [B] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutiérrez and C. Beltrán, 
% "Riemannian Optimization on the Manifold of Unitary and Symmetric
% Matrices with Application to BD-RIS-Assisted Systems",
% IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 
% Barcelona, Spain, 2026, pp. 531–535.
% I. Santamaria, UC Jan. 2026


format compact
clc; clear;
%% Parameters
Nt = 8;                 % Number of transmit antennas
Nr = 4;                 % Number of receive antennas
r = min(Nt,Nr);         % DoF
M = [16,256];           % Number of BD-RIS elements  
SNR = -15:15;                 % Received SNR in dBs
rho = 10.^(SNR/10);           % rho in Eq (5).
PmaxdBm = 0;                  % Pmax (in dBm)
Pmax = 10.^(PmaxdBm/10);      % Pmax
Rxxiso = (Pmax/Nt)*eye(Nt);   % Isotropic covariance matrix (fixed)
NsimMC = 200;              % Number of Monte Carlo simulations (you should use at least 25 to get a reasonable smooth curve)

%% MO algorithm parameters (ICASSP'26 algorithm)
% Parameters for the method in [B]

opt_paramsBDRIS_MO = struct();
opt_paramsBDRIS_MO.maxiter = 100;        % Maximum number of iterations
opt_paramsBDRIS_MO.threshold = 1e-3;     % To check convergence

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
RateMaxDet_iso = zeros(length(M),length(SNR));    % BD-RIS (closed-form Max-Det)
RateMO_iso = zeros(length(M),length(SNR));        % BD-RIS (unitary+ symmetric MO)
RateMaxDet_Rxxopt = zeros(length(M),length(SNR));    % BD-RIS (closed-form Max-Det)
RateMO_Rxxopt = zeros(length(M),length(SNR));        % BD-RIS (unitary+ symmetric MO)

for mm = 1:NsimMC  
    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(M)
        disp(['RIS elements:', num2str(M(dd))])

        %% Generate channels        
        [Hd,G,F] = ChannelsMIMO(M(dd),Nr,Nt,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  
       
        %% Max-Det BDRIS isotropic (Solution in [A])
        BDRISCF = BDRIS_MaxDet(F,G,channelparams.blocked);
        Heq_CF = Hd + F*BDRISCF*G';
        for hh = 1:length(SNR)
            sigma2n = Pmax/(Nt*Nr)*real(norm(F*G','fro')^2)*10.^(-SNR(hh)/10);  % noise variance
            [Ueq,Deq,Veq] = svd(Heq_CF); % svd of eq channel
            Rxxiso_opt = (Pmax/Nr)*(Veq(:,1:r)*Veq(:,1:r)');   
            RateMaxDet_iso(dd,hh) =  RateMaxDet_iso(dd,hh) + log2(real(det(eye(Nr)+ (Heq_CF*Rxxiso_opt *Heq_CF')./sigma2n)));
            % Waterfilling
            [QUdet,VU,pU] = OptTransmitCovMatrix(Heq_CF,sigma2n*eye(Nr),Pmax);
            RateMaxDet_Rxxopt(dd,hh) =  RateMaxDet_Rxxopt(dd,hh) + log2(real(det(eye(Nr)+ (Heq_CF*QUdet*Heq_CF')./sigma2n)));
        end
        
        %% Max-Cap BDRIS isotropic (MO algorithm)
        % Initialization (unitary+symmetric matrix)
        [Uf,Df,Vf] = svd(F);
        [Ug,Dg,Vg] = svd(G');
        Q = Vf*Ug';
        Thetaini = (Q + Q.')/2;  % symmetric projection
        [F_SVD,~,G_SVD] = svd(Thetaini);
        dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);   % THIS IS GENERAL
        Thetaini = F_SVD*F_SVD.';
        for hh =1:length(SNR)
            sigma2n = Pmax/(Nt*Nr)*real(norm(F*G','fro')^2)*10.^(-SNR(hh)/10); % noise variance
            [CfinalMO, CtotalMO, BDRIS_MO] = OptimizeBDRIS_MOUs_FullRank(Hd,F,G,Thetaini,Rxxiso,sigma2n,opt_paramsBDRIS_MO);
            Heq_MO = Hd + F*BDRIS_MO*G';
            [Ueq,Deq,Veq] = svd(Heq_MO); % svd of eq channel
            Rxxiso_opt = (Pmax/Nr)*(Veq(:,1:r)*Veq(:,1:r)');   
            
            RateMO_iso(dd,hh) =  RateMO_iso(dd,hh) + log2(real(det(eye(Nr)+ (Heq_MO*Rxxiso_opt*Heq_MO')./sigma2n)));
            % Waterfilling
            [QUmo,VU,pU] = OptTransmitCovMatrix(Heq_MO,sigma2n*eye(Nr),Pmax);
            RateMO_Rxxopt(dd,hh) =  RateMO_Rxxopt(dd,hh) + log2(real(det(eye(Nr)+ (Heq_MO*QUmo*Heq_MO')./sigma2n)));
        end
    end
end

RateMaxDet_iso = RateMaxDet_iso/NsimMC;
RateMO_iso = RateMO_iso/NsimMC;
RateMaxDet_Rxxopt = RateMaxDet_Rxxopt/NsimMC;
RateMO_Rxxopt = RateMO_Rxxopt/NsimMC;
%% Plot results
figure(30);clf; plot(SNR,RateMO_Rxxopt(1,:), 'r-o','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(SNR,RateMO_iso(1,:),'r--','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMaxDet_Rxxopt(1,:),'b-s','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMaxDet_iso(1,:),'b--','MarkerSize',ms,'LineWidth',lw);

plot(SNR,RateMO_Rxxopt(2,:), 'r-o','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMO_iso(2,:),'r--','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMaxDet_Rxxopt(2,:),'b-s','MarkerSize',ms,'LineWidth',lw);
plot(SNR,RateMaxDet_iso(2,:),'b--','MarkerSize',ms,'LineWidth',lw);

ylabel('Rate (b/s/Hz)');
xlabel('SNR (dB)');
legend('Max-Rate (MO) + waterfilling','Max-Rate (MO), isotropic',...
    'Max-Det (MO) + waterfilling','Max-Det (MO), isotropic' , 'Location','best');
grid on;
hold off


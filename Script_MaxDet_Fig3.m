
% This script shows the performance of the Max-Det BD-RIS solution
% described in [A]. It allows to reproduce the results in Fig. 3 of [A], showing
% the achievable rate vs number of BD-RIS elements. In this scenario, there is a direct link
% and we evaluate the impact of the phase correction term (solution of Eq.
% (17))
%
% We compare the following methods:
% a) Closed-form Max-Det symmetric and unitary solution without phase
% correction
% b) Closed-form Max-Det symmetric and unitary solution wit phase
% correction
% c) Low-Cost app solution in T. Fang and Y. Mao "Low Comlexity beamforming design for
% BDRIS aided multi-user networks," IEEE Comm. Letters, 2024.
% d) Manifold optimization (MO) algorithm (this is the method in [B]): Obtains a full-rank
% unitary and symmetric BD-RIS maximizing capacity.
% e) Random BDRIS
% f) No RIS
%
% [A] I. Santamaria, M. Soleymani, J. Gutierrez, E. Jorswieck, "On the optimal
% symmetric BD-RIS maximizing capacity in a MIMO link" submitted to IEEE
% Transactions on Signal Processing, 2026.

% I. Santamaria, UC Jan. 2026


format compact
clc; clear;
%% Parameters
Nt = 4;                       % Number of transmit antennas
Nr = 4;                       % Number of receive antennas
r = min(Nt,Nr);               % DoF
M = 16;                       % Number of BD-RIS elements
SNR = 10;                     % Received SNR in dBs
a = [1e-3 1e-2 0.1:0.01:20];  % scaling factor of direct link
PmaxdBm = 10;                 % Pmax (in dBm)
Pmax = 10.^(PmaxdBm/10);      % Pmax
Qiso = (Pmax/Nt)*eye(Nt);     % Isotropic covariance matrix (fixed)
NsimMC = 25;                  % Number of Monte Carlo simulations

%% MO algorithm parameters (ICASSP'26 algorithm)
% Parameters for the method in I. Santamaria, M. Soleymani, E. Jorswieck,
% J. Gutierrez, C, Beltran, "Riemannian Optimization on the manifold of
% unitary and symmetric matrices with application to BD-RIS assisted
% systems" ICASSP 2026.

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
channelparams.blocked = 0;         % Set to 1 if direct channel is blocked
% For this scenario we consider a direct
% link of varying strength so this param
% is set to 0
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

B = 1;    % Bandwidth MHz
NF = 0;   % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance


%% Variables to store the rates (b/s/Hz)
RateMaxDet = zeros(1,length(a));         % BD-RIS (closed-form MaxDet)
RateMO = zeros(1,length(a));             % BD-RIS (unitary+ symmetric MO)
RateMaxDetPhase = zeros(1,length(a));    % BD-RIS (closed-form MaxDet + phase correction)
RateLowCost = zeros(1,length(a));        % BD-RIS (Low-Cost)
RateNoRIS = zeros(1,length(a));          % No RIS
Raternd = zeros(1,length(a));            % Rand BD-RIS

for mm = 1:NsimMC

    disp(['Simulation:', num2str(mm)])

    %% Generate channels
    [Hd,G,F] = ChannelsMIMO(M,Nr,Nt,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  % With this version the RiceanFactor also affects the direct links
    
    %% Rnd BDRIS (we calculate the rnd BD-RIS here)
    Q = orth(randn(M,M)+1i*randn(M,M));
    BDRISrnd = Q*Q.';
    for hh = 1:length(a)
        Hdaux = a(hh)*Hd;   % Scale the direct link
        %% Max-Det BDRIS (Solution in [A] w/o phase correction)
        BDRISMaxDet =  BDRIS_MaxDet(F,G,1); % without phase correction (we ste the blocked param to 1)
        HeqMaxDet = Hdaux + F*BDRISMaxDet*G';
        RateMaxDet(hh) =  RateMaxDet(hh) + log2(real(det(eye(Nr)+ (HeqMaxDet*Qiso*HeqMaxDet')./sigma2n)));

        %% Max-Det BDRIS (Solution in [A] with phase correction)
        params = struct;
        params.sigma2n = sigma2n;
        params.Hd = Hdaux;
        params.Rxx = Qiso;
        BDRISMaxDetPhase = BDRIS_MaxDet(F,G,channelparams.blocked,params); % with phase correction
        HeqMaxDetPhase = Hdaux + F*BDRISMaxDetPhase*G';
        RateMaxDetPhase(hh) =  RateMaxDetPhase(hh) + log2(real(det(eye(Nr)+ (HeqMaxDetPhase*Qiso*HeqMaxDetPhase')./sigma2n)));

        %% Max-Cap BDRIS (MO algorithm)
        % Initialization (unitary+symmetric matrix)
        [Uf,Df,Vf] = svd(F);
        [Ug,Dg,Vg] = svd(G');
        Q = Vf*Ug';
        Thetaini = (Q + Q.')/2;  % symmetric projection
        [F_SVD,~,G_SVD] = svd(Thetaini);
        dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);   % THIS IS GENERAL
        Thetaini = F_SVD*F_SVD.';

        [CfinalMO, CtotalMO, BDRIS_MO] = OptimizeBDRIS_MOUs_FullRank(Hdaux,F,G,Thetaini,Qiso,sigma2n,opt_paramsBDRIS_MO);
        Heq_MO = Hdaux + F*BDRIS_MO*G';
        RateMO(hh) =  RateMO(hh) + log2(real(det(eye(Nr)+ (Heq_MO*Qiso*Heq_MO')./sigma2n)));

        %% Max-Cap BDRIS (Low Cost)
        Q = F'*Hdaux*G;
        Thetasymm = (Q+Q.')/2;  % symmetric projection
        [F_SVD,D_SVD,G_SVD] = svd(Thetasymm);
        dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);
        ThetaLowCost = F_SVD*F_SVD.';
        Heq_LowCost = Hdaux + F*ThetaLowCost*G';
        RateLowCost(hh) =  RateLowCost(hh) + log2(real(det(eye(Nr)+ (Heq_LowCost*Qiso*Heq_LowCost')/sigma2n)));

        %% No RIS
        RateNoRIS(hh) =  RateNoRIS(hh) + log2(real(det(eye(Nr)+ (Hdaux*Qiso*Hdaux')/sigma2n)));

        %% Rnd BD-RIS
        Heq_rnd = Hdaux + F*BDRISrnd*G';
        Raternd(hh) = Raternd(hh) + log2(real(det(eye(Nr)+ (Heq_rnd*Qiso*Heq_rnd')/sigma2n)));
    end
end

RateMaxDet = RateMaxDet/NsimMC;
RateMO = RateMO/NsimMC;
RateLowCost = RateLowCost/NsimMC;
RateMaxDetPhase = RateMaxDetPhase/NsimMC;
RateNoRIS = RateNoRIS/NsimMC;
Raternd =  Raternd/NsimMC;

%% Plot results
figure(30);clf; semilogx(a,RateMO, 'r-.','MarkerSize',ms,'LineWidth',lw);
hold on;
semilogx(a,RateMaxDetPhase,'b-','MarkerSize',ms,'LineWidth',lw);
semilogx(a,RateMaxDet,'b--','MarkerSize',ms,'LineWidth',lw);
semilogx(a,RateLowCost,'k-.','MarkerSize',ms,'LineWidth',lw);
semilogx(a,Raternd,'m-','MarkerSize',ms,'LineWidth',lw);
semilogx(a,RateNoRIS,'g--','MarkerSize',ms,'LineWidth',lw);
axis([a(1) a(end) 20 70]);
legend('Max-Rate MO','Max-Det (phase opt.)','Max-Det', 'Low Cost', 'BD-RIS random','No RIS', 'Location','best');
ylabel('Rate (b/s/Hz)');
xlabel('Scale direct link (a)');
grid on;
hold off


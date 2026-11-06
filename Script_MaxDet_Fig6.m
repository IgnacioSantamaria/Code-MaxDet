% This script compares the computational cost
% estimated via de RunTime for the methods in [A] and [B]. It allows to
% reproduce Fig. 6 in [A]
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

% I. Santamaria, UC Nov. 2026


format compact
clc; clear;
%% Parameters
Nt = 8;                 % Number of transmit antennas
Nr = 8;                 % Number of receive antennas
r = min(Nt,Nr);         % DoF
M = 16:16:256;         % Number of BD-RIS elements  
PmaxdBm = 0;                  % Pmax (in dBm)
Pmax = 10.^(PmaxdBm/10);      % Pmax
Qiso = (Pmax/Nt)*eye(Nt);     % Isotropic covariance matrix (fixed)
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

B = 1;    % Bandwidth MHz
NF = 0;   % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance


%% Variables to store the Run Times (secs)
RunTimeMaxDet = zeros(length(M),1);    % BD-RIS (closed-form Max-Det)
RunTimeMO = zeros(length(M),1);        % BD-RIS (unitary+ symmetric MO)

for mm = 1:NsimMC  
    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(M)
        disp(['RIS elements:', num2str(M(dd))])

        %% Generate channels        
        [Hd,G,F] = ChannelsMIMO(M(dd),Nr,Nt,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  
       
        %% Max-Det BDRIS (Solution in [A])
        tic_start = tic;
        BDRISCF = BDRIS_MaxDet(F,G,channelparams.blocked);
        Heq_CF = Hd + F*BDRISCF*G';
        RunTimeMaxDet(dd) =  RunTimeMaxDet(dd) + toc(tic_start);
        
        %% Max-Cap BDRIS (MO algorithm)
        % Initialization (unitary+symmetric matrix)
        tic_start = tic;
        [Uf,Df,Vf] = svd(F);
        [Ug,Dg,Vg] = svd(G');
        Q = Vf*Ug';
        Thetaini = (Q + Q.')/2;  % symmetric projection
        [F_SVD,~,G_SVD] = svd(Thetaini);
        dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);   
        Thetaini = F_SVD*F_SVD.';
        [CfinalMO, CtotalMO, BDRIS_MO] = OptimizeBDRIS_MOUs_FullRank(Hd,F,G,Thetaini,Qiso,sigma2n,opt_paramsBDRIS_MO);
        Heq_MO = Hd + F*BDRIS_MO*G';
        RunTimeMO(dd) =  RunTimeMO(dd) + toc(tic_start);

    end
end

RunTimeMaxDet = RunTimeMaxDet/NsimMC;
RunTimeMO = RunTimeMO/NsimMC;

%% Plot results
figure(30);clf; plot(M,10*log10(RunTimeMO),'r-.','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(M,10*log10(RunTimeMaxDet),'b-s','MarkerSize',ms,'LineWidth',lw);
legend('Max-Rate Unit.+Symm. (MO)', 'Max-Det Unit.+Symm', 'Location','best');
ylabel('10 log_{10} (Run Time [s])');
xlabel('M (Number of elements)');
grid on;
hold off


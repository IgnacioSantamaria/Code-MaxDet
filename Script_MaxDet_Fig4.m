
% This Script illustrates the achievable rate vs number of stems of BD-RIS
% implemente with the q-stem topology with reduced circuit complexity. 
% It allows to reproduce the results of Fig. 4 of [A]. To
% get the B matrix (susceptance matrix) of the q-stem topology we solve a
% linear system of equations. Depending on the rank of the scattering
% matrix this solution can be exavt (if q=2r-1, when the rank is 2r as in
% the Max-Det solution)
%
% I. Santamaria, UC Jan. 2026
%
%
% [A] I. Santamaria, M. Soleymani, J. Gutierrez, E. Jorswieck, "Optimal symmetric low-rank BD-RIS configuration 
% maximizing the determinant of a MIMO link" submitted to IEEE
% Transactions on Signal Processing, 2026.


format compact
clc; clear;
%% Parameters
Nt = 4;                 % Number of transmit antennas
Nr = 4;                 % Number of receive antennas
r = min(Nt,Nr);         % DoF
M = 64;                 % Number of BD-RIS elements  
q = 1:12;               % Number of stems of the q-stem topology
Z0 = 50;                % Impedance matching
PmaxdBm = 0;                  % Pmax (in dBm) of each user
Pmax = 10.^(PmaxdBm/10);      % Pmax
Qiso = (Pmax/Nt)*eye(Nt);     % Isotropic covariance matrix (fixed)
NsimMC = 25;                  % Number of Monte Carlo simulations 

%% MO algorithm parameters (ICASSP'26 algorithm)
% Parameters for the method in I. Santamaria, M. Soleymani, E. Jorswieck,
% J. Gutierrez, C, Beltran, "Riemannian Optimization on the manifold of
% unitary and symmetric matrices with application to BD-RIS assisted
% systems" ICASSP 2026.

opt_paramsBDRIS_MO = struct();
opt_paramsBDRIS_MO.maxiter = 1000;        % Maximum number of iterations
opt_paramsBDRIS_MO.threshold = 1e-5;      % To check convergence

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 8;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Channel parameters
channelparams = struct;
channelparams.blocked = 1;         % Set to 1 if direct channel is blocked
channelparams.RiceRIS = 2;         % Rician factor for the channel btw RIS and Tx/Rx (if Inf-> pure LoS channels)
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
RateMaxDetFC = zeros(size(q));          % BD-RIS (closed-form MaxDet, Fully connected)
RateMaxDetqstem = zeros(size(q));       % BD-RIS (closed-form MaxDet, q-stem connected)
RateMOFC = zeros(size(q));    % BD-RIS (MO, Fully connected)
RateMOqstem = zeros(size(q)); % BD-RIS (MO, q-stem connected)

for mm = 1:NsimMC  % You can also remove the parfor and use for instead

    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(q)  % you may use parfor here
        disp(['q-stem parameter:', num2str(q(dd))])

        %% Generate channels
        [Hd,G,F] = ChannelsMIMO(M,Nr,Nt,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  
        
        %% Max-Det BDRIS (Solution in [A]), fully-connected and q-stem connected architectures
        [ThetaMaxDet, QMaxDet] = BDRIS_MaxDet(F,G,channelparams.blocked);
        BMaxDetqstem = Theta2Bqstem(QMaxDet, Z0, q(dd));  % q-stem susceptance matrix
        ThetaMaxDetqstem =  (eye(M)+1i*Z0*BMaxDetqstem)\(eye(M)-1i*Z0*BMaxDetqstem);  % scattering matrix with the q-stem implementation 
        Heq_FC = Hd + F*ThetaMaxDet*G';
        RateMaxDetFC(dd) =  RateMaxDetFC(dd) + log2(real(det(eye(Nr)+ (Heq_FC*Qiso*Heq_FC')/sigma2n)));
        Heq_qStem = Hd + F*ThetaMaxDetqstem*G';
        RateMaxDetqstem(dd) = RateMaxDetqstem(dd) + log2(real(det(eye(Nr)+ (Heq_qStem*Qiso*Heq_qStem')/sigma2n)));
        
        %% Max-Cap BDRIS (MO Fully-connected & q-stem connected)
        % Initialization (unitary+symmetric matrix)
        [Uf,Df,Vf] = svd(F);
        [Ug,Dg,Vg] = svd(G');
        Qini = Vf*Ug';
        Thetaini = (Qini + Qini.')/2;  % symmetric projection
        [F_SVD,~,G_SVD] = svd(Thetaini);
        dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);   
        Thetaini = F_SVD*F_SVD.';
        
        [CfinalMO, CtotalMO, Theta_MO] = OptimizeBDRIS_MOUs_FullRank(Hd,F,G,Thetaini,Qiso,sigma2n,opt_paramsBDRIS_MO);
        [F_SVD,D_SVD,G_SVD] = svd(Theta_MO);
        dummy = F_SVD'*conj(G_SVD);
        QMO = F_SVD*sqrtm(dummy);
        
        BMOqstem = Theta2Bqstem(QMO, Z0, q(dd));  % q-stem susceptance matrix
        ThetaMOqstem =  (eye(M)+1i*Z0*BMOqstem)\(eye(M)-1i*Z0*BMOqstem);
        
        Heq_MO = Hd + F*Theta_MO*G';
        Heq_MOStemLS = Hd + F*ThetaMOqstem*G';
        RateMOFC(dd) =  RateMOFC(dd) + log2(real(det(eye(Nr)+ (Heq_MO*Qiso*Heq_MO')/sigma2n)));
        RateMOqstem(dd) =   RateMOqstem(dd) + log2(real(det(eye(Nr)+ (Heq_MOStemLS*Qiso*Heq_MOStemLS')/sigma2n)));

    end
end

RateMaxDetFC = RateMaxDetFC/NsimMC;
RateMaxDetqstem = RateMaxDetqstem/NsimMC;
RateMOFC = RateMOFC/NsimMC;
RateMOqstem = RateMOqstem/NsimMC;

%% Plot results
figure(30);clf; plot(q,RateMOFC, 'r--*','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(q,RateMaxDetFC,'b:o','MarkerSize',ms,'LineWidth',lw);
plot(q,RateMaxDetqstem,'b-s','MarkerSize',ms,'LineWidth',lw);
plot(q,RateMOqstem,'r-d','MarkerSize',ms,'LineWidth',lw);
legend('Max-Rate MO (fully-connected)','Max-Det (fully-connected)',...
    'Max-Det (q-stem)','Max-Rate MO (q-stem)', 'Location','best');
ylabel('Rate (b/s/Hz)');
xlabel('q (number of stems)');
grid on;
hold off



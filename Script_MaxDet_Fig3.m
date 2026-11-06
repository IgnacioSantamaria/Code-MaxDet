% This script shows the performance of the Max-Det BD-RIS solution
% described in [A]. It allows to reproduce the results of Fig. 3 of [A], showing 
% the percentage of the achievable rate obtained by the Max-Det solution
% with respect to the max. capacity solution obtained with the MO algoritm
% in [B].
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

% I. Santamaria, UC Aug. 2026


format compact
clc; clear;
%% Parameters
Nt = 6;                 % Number of transmit antennas
Nr = 3;                 % Number of receive antennas
r = min(Nt,Nr);         % DoF
M = 16:16:256;          % Number of BD-RIS elements  
SNR = [-40,-30,-20];    % Received SNR in dBs
rho = 10.^(SNR/10);     % rho in Eq (5) of [A].

Qiso = eye(Nt);         % Isotropic covariance matrix (fixed)
NsimMC = 200;           % Number of Monte Carlo simulations (you should use at least 25 to get a reasonable smooth curve)

%% MO algorithm parameters (ICASSP'26 algorithm)
% Parameters for the method in [B]

opt_paramsBDRIS_MO = struct();
opt_paramsBDRIS_MO.maxiter = 500;        % Maximum number of iterations
opt_paramsBDRIS_MO.threshold = 1e-5;     % To check convergence

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 6;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Variables to store the rates (b/s/Hz)
RateMaxDet = zeros(length(SNR),length(M));    % BD-RIS (closed-form Max-Det)
RateMO = zeros(length(SNR),length(M));        % BD-RIS (unitary+ symmetric MO)

for mm = 1:NsimMC
    disp(['Simulation:', num2str(mm)])
    for hh = 1:length(SNR)
        disp(['SNR:', num2str(SNR(hh))])
        for dd = 1:length(M)
            disp(['RIS elements:', num2str(M(dd))])

            %% Generate channels, we consider nRayleigh channel for this example
            F = (randn(Nr,M(dd))+1i*randn(Nr,M(dd)))/sqrt(2); % iid CN(0,1);
            G = (randn(Nt,M(dd))+1i*randn(Nt,M(dd)))/sqrt(2); % iid CN(0,1);
            Hd = zeros(Nr,Nt);  % direct channel blocked
            %% Max-Det BDRIS (Solution in [A])
            BDRISCF = BDRIS_MaxDet(F,G,1);   % blocked = 1 (no direct link)
            Heq_CF = Hd + F*BDRISCF*G';
            RateMaxDet(hh,dd) =  RateMaxDet(hh,dd) + log2(real(det(eye(Nr)+ rho(hh)*(Heq_CF*Qiso*Heq_CF'))));

            %% Max-Cap BDRIS (MO algorithm in [B])
            % Initialization (unitary+symmetric matrix)
            [Uf,Df,Vf] = svd(F);
            [Ug,Dg,Vg] = svd(G');
            Q = Vf*Ug';
            Thetaini = (Q + Q.')/2;  % symmetric projection
            [F_SVD,~,G_SVD] = svd(Thetaini);
            dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
            F_SVD = F_SVD*sqrtm(dummy);   
            Thetaini = F_SVD*F_SVD.';
            [CfinalMO, CtotalMO, BDRIS_MO] = OptimizeBDRIS_MOUs_FullRank(Hd,F,G,Thetaini,Qiso,1/rho(hh),opt_paramsBDRIS_MO);
            
            Heq_MO = Hd + F*BDRIS_MO*G';
            RateMO(hh,dd) =  RateMO(hh,dd) + log2(real(det(eye(Nr)+ rho(hh)*(Heq_MO*Qiso*Heq_MO'))));

        end
    end
end

RateMaxDet = RateMaxDet/NsimMC;
RateMO = RateMO/NsimMC;

%% Plot results
figure(40);clf; plot(M,100*RateMaxDet(1,:)./RateMO(1,:),'b-s','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(M,100*RateMaxDet(2,:)./RateMO(2,:),'r-*','MarkerSize',ms,'LineWidth',lw);legend('Max-Det Unit.+Symm', 'Max-Rate Unit.+Symm. (MO)', 'Location','best');
plot(M,100*RateMaxDet(3,:)./RateMO(3,:),'k-d','MarkerSize',ms,'LineWidth',lw);legend('Max-Det Unit.+Symm', 'Max-Rate Unit.+Symm. (MO)', 'Location','best');
ylabel('100*R_{det}/R_{max} (%)');
xlabel('M (BD-RIS elements)');
leg = cell(size(rho));

for k = 1:numel(rho)
    exponente = round(log10(rho(k)));
    leg{k} = sprintf('$\\rho = 10^{%d}$', exponente);
end

legend(leg, 'Interpreter', 'latex');
grid on;
hold off



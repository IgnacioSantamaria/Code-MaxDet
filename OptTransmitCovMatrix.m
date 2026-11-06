function [Q,V,p] = OptTransmitCovMatrix(H,Si,Pt,G,Theta,sigma2v,PRIS)

% Description: this function optimizes the transmit covariance matrix for a 
% given MIMO link assisted by a fixed (active or passive) RIS. 
% The capacity formula is C = log det(I + Si^{-1}(H*Q*H')), where Si
% denotes the noise+interference covariance matrix and H is the equivalent
% MIMO channel matrix (including the RIS, i.e., H = Hd + F*Theta*G^H)
% It can be used with active RIS and Passive fixed RISs matrices.
% It can be used with diagonal or beyond-diagonal RIS.
%
% With a passive RIS it performs SVD + power allocation via waterfiling. In
% this case we should use only 3 input parameters: H,Si,Pt.
%
% With an active RIS we have an extra convex contraint so we solve the resulting
% convex problem via cvx. In this case we should use all 7 input parameters.

% Input parameters:
% H: MIMO Channel
% Si: Noise+Interference Covariance matrix,
% Pt: transmit power
% G: channel from Tx to RIS
% Theta: RIS matrix
% sigma2v, PRIS: transmit noise variance + RIS power for active RIS
%
% Output parameters:
% Q,V,p: transmit covariance matrix such that Q = V*diag(p)*V'
%
% Ignacio Santamaria UC Dec. 2023

[Nrx,Ntx] = size(H);
Sitsqrt = sqrtm(inv(Si));
Htint = Sitsqrt*H;
if nargin > 3       % Active RIS case
    if isscalar(H)  % check if the channel is SISO
        %pin = real(trace((Theta*Theta')*(G'*G)));
        %pact = real(trace(sigma2v*(Theta*Theta')));
        %Q = min(Pt, (PRIS-pact/pin));
        Q = Pt;
    else            % MIMO/MISO/SIMO cases
        B = G*Theta';
        A = (B*B');
        b = real(sigma2v*trace(Theta*Theta'));
        cvx_begin sdp quiet
            variable Q(Ntx,Ntx) Hermitian semidefinite
            maximize(log_det(eye(Nrx) + Htint*Q*Htint'))
            subject to
                real(trace(Q)) <= Pt;           % Tx power constraint
                real(trace(Q*A) + b) <= PRIS;   % Active RIS constraint
        cvx_end
    end
else               % passive RIS case
    if isscalar(H) % check if the channel is SISO
        Q = Pt; 
    else           % MIMO/MISO/SIMO cases
        [~,D,V] = svd(Htint,'econ');
        if Ntx>1
            sinrs = (diag(D).^2);
        else
            sinrs = (D.^2);
        end
        sinrs = sinrs(:)';
        [p, ~] = waterfilling(1./sinrs,Pt);
        
        indp = find(p~=0);  % eliminate zeros. Indices with zero allocated power
        p = p(indp);
        V = V(:,indp);
        Q = V*diag(p)*V';
    end
end

%% Output values
[V,D,~] = svd(Q);
p = diag(D);
indp = find(p~=0);  % eliminate zeros. Indices with zero allocated power
p = p(indp);
V = V(:,indp);
Q = V*diag(p)*V';


function [p, mu] = waterfilling(b,Pt)

% Waterfilling algorithm
% Adapted from Palomar & Fonollosa TSP paper "Practical algorithms for a
% family of waterfilling problems"

% Pt = Total power (linear scale)
% b = row vector with the inverse SNRs sigma^2/lambda_i^2, i= 1,..,L

% p = optimal power vector allocated to the L channels
% mu = waterlevel


L = length(b);
a = ones(1,L);

%Reorder substreams with increasing bi
[~,index] = sort(b);
b = b(index);  

%Loop to obtain the optimal waterlevel
b(L+1) = +Inf; 
L_ = L;
while L_>=1
    mu = b(L_);
    if mu < b(L_+1)  &&  mu < ( Pt + sum(b(1:L_)) )/( sum(a(1:L_)) )
        break; %accept hypothesis
    else
        L_= L_-1; %reject hypothesis and form a new one
    end
end

%Compute the definite waterlevel and the power allocation
mu = ( Pt + sum(b(1:L_)) ) / sum(a(1:L_));
p = max(0, mu*a(1:L) - b(1:L));

%Recover the original ordering
p(index)=p;

% Water-level
mu=1/mu;
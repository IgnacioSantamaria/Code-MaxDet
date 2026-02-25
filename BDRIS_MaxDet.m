function [Theta, Q] = BDRIS_MaxDet(F,G,blocked,varargin)

% Description: This function obtains the closed-form solution in [A] that
% maximizes the determinant of the BD-RIS-assisted MIMO channel. 
%
% If the direct link is not blocked (blocked == 0) we include a
% phase term to "align" the through-RIS channel with the direct link. 
% This last step is only effective if there is a direct link, otherwise it's irrelevant,
%
% Input parameters:
% Hd,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% blocked : if == 0 there is a direct link, otherwise the direct link is
% blocked or too weak. If blocked ==0 we should have the following input
% parameters:
% varargin{1}.sigma2n = sigma2n, noise variance 
% varargin{1}.Hd = Hd, direct link 
% varargin{1}.Rxx = Rxx
%
% Output parameters:
% Theta: MxM (low-rank) symmetric BD-RIS
% Q: Takagi factor

% I. Santamaria, UC, Jan 2026
%
% [A] I. Santamaria, M. Soleymani, J. Gutierrez, E. Jorswieck, "Optimal symmetric low-rank BD-RIS configuration 
% maximizing the determinant of a MIMO link" submitted to IEEE
% Transactions on Signal Processing, 2026.

if nargin < 3
    error(message('TooFewInputs'));
elseif nargin == 4
    sigma2n = varargin{1}.sigma2n;
    Hd = varargin{1}.Hd;
    Rxx = varargin{1}.Rxx;
elseif nargin > 4
    error(message('TooManyInputs'));
end

%% Max-Det solution
[Nr,~] = size(F);   % Matrix F is Nr x M
[~,Nt] = size(G);   % Matrix G is Nt x M
r = min(Nt,Nr);
[~,~,VF] = svd(F);
[~,~,VG] = svd(G);
VF1 = VF(:,1:r); 
VG1 = VG(:,1:r);
A = [VF1 conj(VG1)];
[U,~,~] = svd(A,'econ');
U1 = U(:,1:r);
U2 = U(:,r+1:2*r);
Theta =  U1*U1.'- U2*U2.';            % This is the rank-2r Max-Det solution
Bldiag = blkdiag(eye(r), -1i*eye(r));
Q = [U1 U2]*Bldiag;                   % This is the Q factor (Takagi factor)

if blocked == 0     % If direct link is NOT blocked
   
    %% Phase term
    HeqBDRIS = F*Theta *G';
    A = eye(Nr) + (Hd*Rxx*Hd'+ HeqBDRIS*Rxx*HeqBDRIS')/sigma2n;
    B = Hd*Rxx*HeqBDRIS'/sigma2n;
    fun = @(x)-log2(real(det(A + exp(-1i*x)*B + exp(1i*x)*B')));
    x = fminbnd(fun,0,2*pi);  % Optimize phase
    Theta = exp(1i*x)*Theta;  % Rotated BD-RIS   
    Q = exp(1i*x/2)*Q;        % New Q factor

end

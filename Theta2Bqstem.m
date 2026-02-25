function B = Theta2Bqstem(Q, Z0, q)

% Algorithm to find a structured B (q-stem)
% approximating Theta = Q * Q' in LS sense.
% Inputs: Q (M x r complex, Q' * Q = eye(r)), Z0 (scalar >0), q (1 <= q <= M)
% Output: B (M x M real symmetric, q-stem structured)
%
% I. Santamaria Jan 2026


M = size(Q,1);
r = size(Q,2);
if r > M || q > M || q <1
    error('Invalid dimensions: r <= M, 1 <= q <= M');
end

Q_R = Q;
Q_Rs = conj(Q_R);

% Equation (32) with sign flip for correctness 
C = 1i * Z0 * (Q_Rs + Q_R);
D = -(Q_R - Q_Rs);  % Equivalent to BC = -D with original D

ImC = imag(C);
ImD = imag(D);

% Build R (M^2 x num_ind): maps independent params to vec(B)
m = M - q;
num_ind = q*(q+1)/2 + q*m + m;
R_mat = zeros(M*M, num_ind);

param = 0;
% Upper tri B(1:q,1:q)
for i=1:q
    for j=i:q
        param = param + 1;
        idx_ij = (j-1)*M + i;
        R_mat(idx_ij, param) = 1;
        if i ~= j
            idx_ji = (i-1)*M + j;
            R_mat(idx_ji, param) = 1;
        end
    end
end
% Entries B(1:q, q+1:M)
for kk=1:m
    for l=1:q
        param = param + 1;
        j = q + kk;
        idx_lj = (j-1)*M + l;
        idx_jl = (l-1)*M + j;
        R_mat(idx_lj, param) = 1;
        R_mat(idx_jl, param) = 1;
    end
end
% Diag B(q+1:M, q+1:M)
for kk=1:m
    param = param + 1;
    j = q + kk;
    idx_jj = (j-1)*M + j;
    R_mat(idx_jj, param) = 1;
end

A = kron(ImC.', eye(M)) * R_mat;
z = ImD(:);
b = pinv(A) * z;  % Robust for over/under-determined
vecB = R_mat * b;
B = reshape(vecB, M, M);
% Symmetry check (should be symmetric real)
B = real(0.5 * (B + B.'));


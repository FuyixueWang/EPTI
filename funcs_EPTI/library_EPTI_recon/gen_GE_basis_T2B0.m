function [U, X0, S] = gen_GE_basis_T2B0(N, ETL, T0, dt, T2vals,B0vals)

% Generate subspace basis for GE EPTI data
% modified based on Jon Tamir's  basis generation function
% Zijing Dong, 2020, MGH 

% Add simulating B0 phase in the subpsace bases for B0 updating
% EPTI recosntruction
% Zijing Dong, 2021, MGH

% Fuyixue Wang, updated for C2P, 2023

if length(T2vals) > N
    idx = randperm(length(T2vals));
    T2vals = T2vals(idx(1:N));
end

TEs=(dt:dt:ETL*dt)+T0;

LT1 = length(T2vals);
LT2 = length(B0vals);

X0 = zeros(ETL, LT1, LT2);

for ii=1:LT1
    R2 = 1/T2vals(ii);
    for jj=1:LT2
    X0(:,ii,jj)=exp(-TEs(:)*R2).*exp(1i*2*pi*B0vals(jj)*TEs(:));
    end
end
X0=reshape(X0,ETL,[]);
[U, S, ~] = svd(X0, 'econ');
S = diag(S);
end
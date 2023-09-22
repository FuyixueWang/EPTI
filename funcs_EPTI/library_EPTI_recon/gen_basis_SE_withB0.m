function [U, X0, S] = gen_basis_SE_withB0(N1,N2, TEs_SE, TE_SE, T2svals,T2vals,B0vals)
% Generate subspace basis for SE EPTI data
% modified based on Jon Tamir's  basis generation function
% Zijing Dong, 2020, MGH 

% Add simulating B0 phase in the subpsace bases for B0 updating
% EPTI recosntruction
% Zijing Dong, 2021, MGH

% Fuyixue Wang, updated for C2P, 2023, MGH

if length(T2svals) > N1
    idx = randperm(length(T2svals));
    T2svals = T2svals(idx(1:N1));
end
if length(T2vals) > N2
    idx = randperm(length(T2vals));
    T2vals = T2vals(idx(1:N1));
end

np_SE=size(TEs_SE,1);
np_SE1=find(TEs_SE>TE_SE,1)-1;
nt=np_SE;

LT1 = length(T2svals);
LT2 = length(T2vals);
LT3 = length(B0vals);
                     
X0 = zeros(nt, LT1*LT2*LT3);
num=1;
for jj=1:LT2
    R2 = 1/T2vals(jj);
    for kk=1:LT1
        R2s = 1/T2svals(kk);
        if (R2s>=R2) && (R2s/R2<4)
            for b0 = 1:LT3
                R2_p=R2s-R2;
                X0(1:np_SE1,num)=exp(-TE_SE*R2_p).*exp(-TEs_SE(1:np_SE1)*(R2-R2_p));
                X0(np_SE1+1:end,num)=exp(TE_SE*R2_p).*exp(-TEs_SE(np_SE1+1:end)*R2s);
                X0(:,num) = X0(:,num).*exp(1i*2*pi*B0vals(b0)*(TEs_SE(:)-TE_SE));
                num=num+1;
            end
        end
    end
end

X0=X0(:,1:num-1);

[U, S, ~] = svd(X0, 'econ');
S = diag(S);

end
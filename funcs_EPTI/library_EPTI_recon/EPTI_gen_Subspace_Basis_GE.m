function [Phi,TEs_GE,S] = EPTI_gen_Subspace_Basis_GE(parameters,K,show_fig)
% Function for subspace basis generation
% used in subspace EPTI reconstruction
% Zijing Dong 8/05/2019

% GE EPTI sequence and parameter calculation
% Fuyixue Wang, 2021

% clean up and optimized for C2P, 2023 
% Fuyixue Wang, Zijing Dong, 2023, MGH
%%
if nargin<3
    show_fig = 0;
end

T2s=[10:1:50,52:2:100,102:5:200,220:20:300];
% T2s=[10:1:50,52:2:100,102:5:200];
% T2s=[2:1:50,52:2:100,102:5:200,220:20:300];
B0vals = 0;

%% Simulate and generate basis
dt = parameters.iEffectiveEpiEchoSpacing*1e-6; % echo spacing
nechoGE = parameters.nechoGE; % number of time points
necho_select = nechoGE;
if (parameters.pf_echo==0) 
    t0 = parameters.alTE(1)*1e-6 - (parameters.lPhaseEncodingLines/2)*dt;
% elseif (length(parameters.alTE)==1)
%     t0 = parameters.alTE(1)*1e-6 - (parameters.nechoGE/2)*dt;
else
    t0 = parameters.alTE(1)*1e-6 - (floor(floor(parameters.lPhaseEncodingLines*parameters.pf_echo)/2))*dt;
end
TEs_GE=(dt:dt:nechoGE*dt)+t0;

N = 256; % maximum number of unique T2 values for training
T2s=T2s/1000;
TEs_GE = TEs_GE(:);
TEs_GE = TEs_GE(1:necho_select);
[U, X, S] = gen_GE_basis_T2B0(N, length(TEs_GE), 0, dt, T2s, B0vals);
Phi = U(:,1:K);
%%
Z = Phi*Phi'*X;
err = norm(X(:) - Z(:)) / norm(X(:));
fprintf('Relative norm of error of the generated subspace basis: %.6f\n', err);
if show_fig==1
    figure;
    plot(TEs_GE(end-size(X,1)+1:end)*1000, real(X(end-size(X,1)+1:end,:)), 'linewidth', 2);

    % xlim([dt*1000, dt*1000*nechoGE]);
    xlabel('Virtual echo time (ms)');
    ylabel('Signal value');
    ftitle('Signal evolutions for distribution of T2 values', 24)
    faxis;

    figure;
    subplot(1,2,1); plot(real(Phi), 'linewidth', 3);
    ftitle('Subspace curves Real Part', 24)
    faxis;
    subplot(1,2,2); plot(imag(Phi), 'linewidth', 3);
    ftitle('Subspace curves Real Part', 24)
    faxis;
end

end


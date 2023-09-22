function [Phi,TEs_SE,S] = EPTI_gen_Subspace_Basis_SE(parameters,K,show_fig)
% Function for subspace basis generation
% used in subspace EPTI reconstruction
% Zijing Dong 8/05/2019

% Basis generation for SE EPTI sequence
% Fuyixue Wang, 2020, MGH

% Add simulating B0 phase in the subpsace bases for B0 updating
% EPTI recosntruction
% Zijing Dong, 2021, MGH

% clean up and optimized for C2P, 2023 
% Fuyixue Wang, 2023, MGH
%%
if nargin<3
    show_fig = 0;
end

T2=[5:1:50,52:2:100,102:5:200,220:20:400];
T2s = T2;
B0vals= 0; % unit Hz%% Simulate and generate basis
dt = parameters.iEffectiveEpiEchoSpacing*1e-6; % echo spacing
nechoSE = parameters.nechoSE; % number of time points
SE_echoindex = floor(parameters.lPhaseEncodingLines/2)+1;
t0 = parameters.alTE(2)*1e-6 - SE_echoindex*dt; % % time for SE echo
tes_SE=(dt:dt:nechoSE*dt)+t0;
TE_SE=parameters.alTE(2)*1e-6 ;

N1 = 256; % maximum number of unique T2 values for training
N2 = 256; % maximum number of unique T2 values for training
T2s=T2s/1000;
T2=T2/1000;

[U, X, S] = gen_basis_SE_withB0(N1, N2, tes_SE(:), TE_SE, T2s,T2,B0vals);
Phi = U(:,1:K);
TEs_SE = tes_SE(:)-TE_SE;
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


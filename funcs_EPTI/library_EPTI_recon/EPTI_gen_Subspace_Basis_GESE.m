function [Phi,tes_SE] = EPTI_gen_Subspace_Basis_GESE(parameters,K,show_fig)
% Function for subspace basis generation
% used in subspace EPTI reconstruction
% Zijing Dong 8/05/2019

% Basis generation for GE-SE EPTI sequence
% Fuyixue Wang, 2020, MGH
%%
if nargin<3
    show_fig = 0;
end

T2=[5:1:50,52:2:100,102:5:200,220:20:400];
T2s = T2;
dt = parameters.iEffectiveEpiEchoSpacing*1e-6; % echo spacing
nechoSE = parameters.nechoSE; % number of time points
SE_echoindex = floor(nechoSE/2)+1;
t0 = parameters.alTE(2)*1e-6 - SE_echoindex*dt; % % time for SE echo
TEs_SE=(dt:dt:nechoSE*dt)+t0;
TE_SE=parameters.alTE(2)*1e-6 ;

nechoGE = parameters.nechoGE; % number of time points
t0 = parameters.alTE(1)*1e-6 - (parameters.lPhaseEncodingLines/2)*dt; % % time for SE echo
TEs_GE=(dt:dt:nechoGE*dt)+t0;
%% Simulate and generate basis
K = 4; % subspace size
N1 = 256; % maximum number of unique T2 values for training
N2 = 256; % maximum number of unique T2 values for training
deltaS0 = [0.8,0.9,1,1.1,1.2];
T2s=T2s/1000;
T2=T2/1000;

[U, X] = gen_GE_basis_GESE(N1, N2, TEs_GRE(:), TEs_SE(:), TE_SE, T2s,T2,deltaS0);
Phi = U(:,1:K);
%%
Z = Phi*Phi'*X;
err = norm(X(:) - Z(:)) / norm(X(:));
fprintf('Relative norm of error: %.6f\n', err);
show_fig=0;
if show_fig==1
    figure;
    plot(TEs_GRE(:)*1000, real(X(1:size(TEs_GRE,2),:)), 'linewidth', 2);
    hold on;
    plot(TEs_SE(:)*1000, real(X(size(TEs_GRE,2)+1:end,:)), 'linewidth', 2);

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
    %% Project the signal evolutions onto the subspace
    figure;
    plot(TEs_GRE(:)*1000, real(Z(1:size(TEs_GRE,2),6000)), 'linewidth', 2);
    hold on;
    plot(TEs_SE(:)*1000, real(Z(size(TEs_GRE,2)+1:end,6000)), 'linewidth', 2);
    % xlim([dt*1000, dt*1000*nechoGE]);
    xlabel('Virtual echo time (ms)');
    ylabel('Signal value');
    ftitle('Projected signal evolutions', 24)
    faxis;
end

end


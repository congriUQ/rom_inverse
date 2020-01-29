%Script that runs the inverse problem

nIterations = 10000;
testSample = 3404;
uf_true = rom.trainingDataMatfile.uf(:, testSample);
xi_true = rom.trainingDataMatfile.xi(testSample, :)

%Variational Inference parameters
%First elements correspond to log lambda_c, last elements to xi
dim_VI = rom.coarseMesh.nEl + rom.nBochnerBasis;
VI_params.mu = zeros(1, dim_VI);
VI_params.sigma = .1*ones(1, dim_VI);
VI_params_vec = [VI_params.mu, -2*log(VI_params.sigma)];
log_emp_dist = @(x) log_empirical_distribution(x(1:rom.coarseMesh.nEl)', x((rom.coarseMesh.nEl + 1):end), rom, uf_true);

%Plot stuff
f = figure('units','normalized','outerposition',[0 0 1 1]);
ax = subplot(1, 4, 1, 'Parent', f);
[Xf, Yf] = meshgrid(linspace(0, 1, 257));
s = surf(Xf, Yf, reshape(uf_true, 257, 257), 'linestyle', 'none', 'Parent', ax);
axis(ax, 'square');
axis(ax, 'tight');
ax.XTickLabel = {};
ax.YTickLabel = {};
ax.Title.String = "Solution $u_f$ and reconstruction";
ax.View = [-120, -5];
hold(ax, 'on');

ax(2) = subplot(1, 4, 2, 'Parent', f);
im = imagesc(reshape(rom.getDiscretizedConductivityField(xi_true), 256, 256));
axis(ax(2), 'square');
ax(2).GridLineStyle = 'none';
ax(2).XTick = {};
ax(2).YTick = {};
ax(2).Title.String = "Conductivity from $\xi$";
drawnow;




log_lambda_c_opt = 0*mean(rom.XMean, 2);
xi_opt = zeros(1, rom.nBochnerBasis);
% xi_opt = normrnd(0, 1, 1, rom.nBochnerBasis);
% xi_opt = xi_true;

options = optimoptions(@fminunc, 'Display', 'off', 'MaxIterations', 5);

rom.theta_cf.sumLogS = sum(log(rom.theta_cf.S));
for i = 1:nIterations
%     neg_log_p_log_lambda_c_opt = @(log_lambda_c) neg_log_p_lambda_c(rom, xi_opt, log_lambda_c, uf_true);
%     log_lambda_c_opt = fminunc(neg_log_p_log_lambda_c_opt, log_lambda_c_opt, options)
%     
%     neg_log_p_xi_opt = @(xi) neg_log_p_xi(rom, xi, log_lambda_c_opt);
%     xi_opt = fminunc(neg_log_p_xi_opt, xi_opt, options)
    
    
    [VI_params, VI_params_vec] =...
        efficientStochOpt(VI_params_vec, log_emp_dist, 'diagonalGauss', 1e-3, dim_VI);
    fprintf("Current mu_lambda_c == ")
    disp(VI_params.mu(1:rom.coarseMesh.nEl))
    fprintf("Current mu_xi == ")
    disp(VI_params.mu((rom.coarseMesh.nEl + 1):end))
    fprintf("Current sigma_xi == %f\n")
    disp(VI_params.sigma((rom.coarseMesh.nEl + 1):end))
    % Hack for plotting
    xi_opt = VI_params.mu((rom.coarseMesh.nEl + 1):end);
    log_lambda_c_opt = VI_params.mu(1:rom.coarseMesh.nEl);
    
    
    
    
    %Live plot conductivity field
    ax(3) = subplot(1, 4, 3, 'Parent', f);
    im = imagesc(reshape(rom.getDiscretizedConductivityField(xi_opt), 256, 256));
    axis(ax(3), 'square');
    ax(3).GridLineStyle = 'none';
    ax(3).XTick = {};
    ax(3).YTick = {};
    ax(3).Title.String = "Reconstructed conductivity";
    
    ax(4) = subplot(1, 4, 4, 'Parent', f);
    im = imagesc(reshape(exp(log_lambda_c_opt), rom.coarseMesh.nElX, rom.coarseMesh.nElY));
    axis(ax(4), 'square');
    ax(4).GridLineStyle = 'none';
    ax(4).XTick = {};
    ax(4).YTick = {};
    ax(4).Title.String = "Effective conductivity";
    
    %exp(log_lambda_c_opt) needs to be a column vector
    FEMout = heat2d(rom.coarseMesh, exp(log_lambda_c_opt)');
    XXc = meshgrid(linspace(0, 1, rom.coarseMesh.nElX + 1));
    [~, YYc] = meshgrid(linspace(0, 1, rom.coarseMesh.nElY + 1));
    if i == 1
        s(2) = surf(XXc, YYc, reshape(FEMout.u, rom.coarseMesh.nElX + 1, rom.coarseMesh.nElY + 1), 'Parent', ax(1));
        s(2).LineStyle = 'none';
        s(2).FaceColor = 'b';
    else
        s(2).ZData = reshape(FEMout.u, rom.coarseMesh.nElX + 1, rom.coarseMesh.nElY + 1);
    end
    drawnow;
end



function [log_p, dlog_p_dxi] = neg_log_p_xi(rom, xi, log_lambda_c)
    %gradient
    [designMatrix, d_designMatrix] = rom.computeDesignMatrix("inverse", false, xi);
    conductivity = rom.getDiscretizedConductivityField(xi);
    dconductivity_dxi = (conductivity - rom.lowerConductivity).*rom.randFieldGenMatrix;
    d_designMatrix_dxi = d_designMatrix*dconductivity_dxi;
    %reshape for vectorization
    d_designMatrix_dxi = reshape(d_designMatrix_dxi, numel(rom.theta_c.theta), rom.nBochnerBasis*rom.coarseMesh.nEl)';
    d_designMatrix_dxi_times_theta_c = reshape(d_designMatrix_dxi*rom.theta_c.theta, rom.coarseMesh.nEl, rom.nBochnerBasis);
    dlog_p_dxi = -(rom.theta_c.Sigma\(log_lambda_c - designMatrix{1}*rom.theta_c.theta))'*...
        d_designMatrix_dxi_times_theta_c + xi;
    
    %function
    log_p = .5*sum(((log_lambda_c - designMatrix{1}*rom.theta_c.theta).^2)./(diag(rom.theta_c.Sigma)));
end


function [log_p, dlog_p_dlog_lambda_c] = neg_log_p_lambda_c(rom, xi, log_lambda_c, uf_true)
    designMatrix = rom.computeDesignMatrix("inverse", false, xi);
    [log_p, dlog_p_dlog_lambda_c] = log_q_i(log_lambda_c, uf_true, rom.theta_cf, rom.theta_c, designMatrix{1},...
        rom.coarseMesh, rom.conductivityTransformation, false);
    log_p = -log_p;
    dlog_p_dlog_lambda_c = -dlog_p_dlog_lambda_c;
end

function [log_p, d_log_p] = log_empirical_distribution(log_lambda_c, xi, rom, uf_true)
    [designMatrix, d_designMatrix] = rom.computeDesignMatrix("inverse", false, xi);
    conductivity = rom.getDiscretizedConductivityField(xi);
    dconductivity_dxi = (conductivity - rom.lowerConductivity).*rom.randFieldGenMatrix;
    d_designMatrix_dxi = d_designMatrix*dconductivity_dxi;
    %reshape for vectorization
    d_designMatrix_dxi = reshape(d_designMatrix_dxi, numel(rom.theta_c.theta), rom.nBochnerBasis*rom.coarseMesh.nEl)';
    d_designMatrix_dxi_times_theta_c = reshape(d_designMatrix_dxi*rom.theta_c.theta, rom.coarseMesh.nEl, rom.nBochnerBasis);
    dlog_p_dxi = -(rom.theta_c.Sigma\(log_lambda_c - designMatrix{1}*rom.theta_c.theta))'*...
        d_designMatrix_dxi_times_theta_c + xi;
    log_p_xi = .5*sum(((log_lambda_c - designMatrix{1}*rom.theta_c.theta).^2)./(diag(rom.theta_c.Sigma)));
    [log_p_lambda_c, dlog_p_dlog_lambda_c] = log_q_i(log_lambda_c, uf_true, rom.theta_cf, rom.theta_c, designMatrix{1},...
        rom.coarseMesh, rom.conductivityTransformation, false);
    
%     [log_p_lambda_c, dlog_p_dlog_lambda_c] = neg_log_p_lambda_c(rom, xi, log_lambda_c, uf_true);
%     [log_p_xi, dlog_p_dxi] = neg_log_p_xi(rom, xi, log_lambda_c);
    log_p = -log_p_lambda_c + log_p_xi;
    d_log_p = [dlog_p_dlog_lambda_c', -dlog_p_dxi];
end





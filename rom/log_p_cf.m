function [log_p, d_log_p, Tc] = log_p_cf(Tf_i_minus_mu, coarseMesh, Xi, theta_cf, condTransOpts, onlyGrad)
%Coarse-to-fine map
%ignore constant prefactor
%log_p = -.5*logdet(S, 'chol') - .5*(Tf - mu)'*(S\(Tf - mu));
%diagonal S

conductivity = conductivityBackTransform(Xi, condTransOpts);
FEMout = heat2d(coarseMesh, conductivity);

Tc = FEMout.u;
Tf_i_minus_mu_minus_WTc = Tf_i_minus_mu - theta_cf.W*Tc;
if onlyGrad
    %to avoid computation overhead if only d_log_p is needed
    log_p = [];
else
    %only for diagonal S!
    log_p = -.5*(theta_cf.sumLogS + (Tf_i_minus_mu_minus_WTc)'*(theta_cf.Sinv_vec.*(Tf_i_minus_mu_minus_WTc)));
end

if nargout > 1
    %Gradient of FEM equation system w.r.t. conductivities

%     d_r = FEMgrad(FEMout, domainc, conductivity);
    d_r = FEMgrad2(FEMout.naturalTemperatures, coarseMesh);
%     d_rx = d_r;
    if strcmp(condTransOpts.type, 'log')
        %We need gradient of r w.r.t. log conductivities X, multiply 
        %each row with resp. conductivity
        d_rx(1:coarseMesh.nEl, :) = conductivity.*d_r(1:coarseMesh.nEl, :);
    elseif strcmp(condTransOpts.type, 'logit')
        %We need gradient w.r.t. x, where x is
        %- log((lambda_up - lambda_lo)/(lambda - lambda_lo) - 1)
        X = conductivityTransform(conductivity, condTransOpts);
        dLambda_dX = (condTransOpts.limits(2) - condTransOpts.limits(1))./(exp(X) + 2 + exp(-X));
        d_rx = dLambda_dX.*d_r;
    elseif strcmp(condTransOpts.type, 'log_lower_bound')
        %transformation is X = log(Lambda - lambda_lo)
        dLambda_dX = conductivity - condTransOpts.limits(1);
        d_rx(1:coarseMesh.nEl, :) = diag(dLambda_dX)*d_r(1:coarseMesh.nEl, :);
    elseif strcmp(condTransOpts.type, 'square')
        %We need gradient of r w.r.t. sqrt conductivities X. d/dX = 2X d/dlambda
        d_rx(1:coarseMesh.nEl, :) = diag(2*Xi(1:coarseMesh.nEl))*d_r(1:coarseMesh.nEl, :);
    else
        error('Unknown conductivity transformation')
    end
    adjoints = get_adjoints(FEMout.globalStiffness, theta_cf, coarseMesh, Tf_i_minus_mu_minus_WTc);
    d_log_p = - d_rx*adjoints;

    
    %Finite difference gradient check
    FDcheck = false;
    if FDcheck
        disp('Gradient check log p_cf')
        d = 1e-8;
        FDgrad = zeros(coarseMesh.nEl, 1);
        for e = 1:coarseMesh.nEl
            conductivityFD = conductivity;
            conductivityFD(e) = conductivityFD(e) + d;
            
            DFD = zeros(2, 2, coarseMesh.nEl);
            for j = 1:coarseMesh.nEl
                DFD(:, :, j) =  conductivityFD(j)*eye(2);
            end
            FEMoutFD = heat2d(coarseMesh, DFD);
            checkLocalStiffness = false;
            if checkLocalStiffness
                k = FEMout.diffusionStiffness(:, :, e);
                kFD = FEMoutFD.diffusionStiffness(:, :, e);
                d_k = FEMout.diffusionStiffness(:, :, e)/conductivity(e);
                d_kFD = (FEMoutFD.diffusionStiffness(:, :, e) - FEMout.diffusionStiffness(:, :, e))/d;
                relgrad_k = d_k./d_kFD
                pause
            end
            TcFD = FEMoutFD.Tff';
            TcFD = TcFD(:);
            
            WTcFD = theta_cf.W*TcFD;
            log_pFD = -.5*(theta_cf.sumLogS + (Tf_i_minus_mu - WTcFD)'*(theta_cf.Sinv_vec.*(Tf_i_minus_mu - WTcFD)));
            if strcmp(condTransOpts.type, 'log')
                FDgrad(e) = conductivity(e)*(log_pFD - log_p)/d;
            elseif strcmp(condTransOpts.type, 'logit')
                FDgrad(e) = dLambda_dX(e)*(log_pFD - log_p)/d;
            elseif strcmp(condTransOpts.type, 'log_lower_bound')
                FDgrad(e) = dLambda_dX(e)*(log_pFD - log_p)/d;
            else
                error('Unknown conductivity transformation')
            end
        end
        relgrad = FDgrad./d_log_p
        plot(1:numel(FDgrad), FDgrad, 1:numel(FDgrad), d_log_p)
        axis square
        drawnow
        pause(1)
        if(norm(relgrad - 1) > 1e-1)
            log_p
            log_pFD
            d_log_p
            FDgrad
            diff = log_pFD - log_p
        end
    end %FD gradient check
    
    
end

end


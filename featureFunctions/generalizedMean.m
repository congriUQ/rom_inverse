function [m, dm_dx, dm_dz] = generalizedMean(x, mask, z)
%Computes the generalized mean of x to parameter z

x_mask = x(mask);       %e.g., pixels of a coarse cell
if(size(x_mask, 3) > 1)
    error('generalized mean only defined for vector/matrix input')
end

if z == 0
    if any(x_mask == 0)
        m = 0;
    else
        m = exp(mean(mean(log(x_mask))));
    end
else
    sum_x_to_z = sum(sum(x_mask.^z));
    m = ((1/numel(x_mask))*sum_x_to_z)^(1/z);
end

if nargout > 1
    %compute derivative w.r.t. input x
    if z == 0
        if ~any(x_mask)
            error("Geometric mean gradient not yet implemented for zero components")
        else
            dm_dx_mask = m./(n*x_mask(:));
        end
    else
        dm_dx_mask = (m/sum_x_to_z)*(x_mask(:).^(z - 1));
    end
    dm_dx = zeros(size(x));
    dm_dx(mask) = dm_dx_mask;
end

if nargout > 2
    %Compute derivative w.r.t. parameter z
    if z == 0
        if any(x_mask == 0)
            dm_dz = 0;
        else
            dm_dz = .5*m*(1/numel(x_mask))*(sum(sum(log(x_mask).^2)) - (1/numel(x_mask))*(sum(sum(log(x_mask))))^2);
        end
    else
        dm_dz = (m/z)*(sum(sum((x_mask.^z).*log(x_mask)))/sum(sum(x_mask.^z))  - (1/z)*log(mean(mean(x_mask.^z))));
    end
end

end


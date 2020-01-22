function [m, dm_dx, dm_dz] = generalizedMean(x, z)
%Computes the generalized mean of x to parameter z

if(size(x, 3) > 1)
    error('generalized mean only defined for vector/matrix input')
end

if z == 0
    if any(x == 0)
        m = 0;
    else
        m = exp(mean(mean(log(x))));
    end
else
    sum_x_to_z = sum(sum(x.^z));
    m = ((1/numel(x))*sum_x_to_z)^(1/z);
end

if nargout > 1
    %compute derivative w.r.t. input x
    if z == 0
        if any(x) == 0
            error("Geometric mean gradient not yet implemented for zero components")
        else
            dm_dx = m./(n*x);
        end
    else
        dm_dx = (m/sum_x_to_z)*(x.^(z - 1));
    end
end

if nargout > 2
    %Compute derivative w.r.t. parameter z
    if z == 0
        if any(x == 0)
            dm_dz = 0;
        else
            dm_dz = .5*m*(1/numel(x))*(sum(sum(log(x).^2)) - (1/numel(x))*(sum(sum(log(x))))^2);
        end
    else
        dm_dz = (m/z)*(sum(sum((x.^z).*log(x)))/sum(sum(x.^z))  - (1/z)*log(mean(mean(x.^z))));
    end
end

end


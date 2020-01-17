function [bcPressure, bcFlux] = bcCoeffs2Fun(bcCoeffs)
%Uses bc coefficients bcCoeffs to construct handle functions to the
%corresponding pressure and boundary flux fields

bcPressure = @(x) bcCoeffs{i}(1) + bcCoeffs(2)*x(1) + bcCoeffs(3)*x(2) + bcCoeffs(4)*x(1)*x(2);
%lower bound
bcFlux{1} = @(x) -(bcCoeffs(3) + bcCoeffs(4)*x);
%right bound
bcFlux{2} = @(y) (bcCoeffs(2) + bcCoeffs(4)*y);
%upper bound
bcFlux{3} = @(x) (bcCoeffs(3) + bcCoeffs(4)*x);
%left bound
bcFlux{4} = @(y) -(bcCoeffs(2) + bcCoeffs(4)*y);

end


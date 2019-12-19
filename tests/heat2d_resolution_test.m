% performance test for heat2d
clear;
addpath('./heatFEM');
addpath('./FEMgradient');

%Set up boundary condition functions
bc = [100 -200 150 -250];
boundaryTemperature = @(x) bc(1) + bc(2)*x(1) + bc(3)*x(2) + bc(4)*x(1)*x(2);
boundaryHeatFlux{1} = @(x) -(bc(3) + bc(4)*x);      %lower bound
boundaryHeatFlux{2} = @(y) (bc(2) + bc(4)*y);       %right bound
boundaryHeatFlux{3} = @(x) (bc(3) + bc(4)*x);       %upper bound
boundaryHeatFlux{4} = @(y) -(bc(2) + bc(4)*y);      %left bound

c_hi = 1000;

% set up domain object
nRes = 20;
gridVectorX = (1/nRes)*ones(1, nRes);
gridVectorY = (1/nRes)*ones(1, nRes);

mesh = Domain(nRes, nRes, gridVectorX, gridVectorY);
mesh.useConvection = false;

mesh = mesh.setBoundaries([2:(2*nRes + 2*nRes)], boundaryTemperature,...
    boundaryHeatFlux);

%shrink mesh object for memory efficiency
mesh.boundaryNodes = [];
mesh.naturalNodes = [];
mesh.boundaryElements = [];
mesh.naturalBoundaries = [];
mesh.boundaryType = [];
mesh.lx = [];
mesh.ly = [];
mesh.lElX = [];
mesh.lElY = [];
mesh.cum_lElX = [];
mesh.cum_lElY = [];
mesh.AEl = [];
mesh.nEq = [];
mesh.lc = [];
mesh.Bvec =[];
mesh.useConvection = [];
mesh.essentialBoundary = [];
mesh.LocalNode = [];
mesh.fs = [];
mesh.fh = [];

D = ones(nRes^2, 1);
D(1:(.5*nRes^2)) = c_hi;
out1 = heat2d_v2(mesh, D);






% set up domain object
nRes2 = 80;
gridVectorX2 = (1/nRes2)*ones(1, nRes2);
gridVectorY2 = (1/nRes2)*ones(1, nRes2);

mesh2 = Domain(nRes2, nRes2, gridVectorX2, gridVectorY2);
mesh2.useConvection = false;

mesh2 = mesh2.setBoundaries([2:(2*nRes2 + 2*nRes2)], boundaryTemperature,...
    boundaryHeatFlux);

%shrink mesh object for memory efficiency
mesh2.boundaryNodes = [];
mesh2.naturalNodes = [];
mesh2.boundaryElements = [];
mesh2.naturalBoundaries = [];
mesh2.boundaryType = [];
mesh2.lx = [];
mesh2.ly = [];
mesh2.lElX = [];
mesh2.lElY = [];
mesh2.cum_lElX = [];
mesh2.cum_lElY = [];
mesh2.AEl = [];
mesh2.nEq = [];
mesh2.lc = [];
mesh2.Bvec =[];
mesh2.useConvection = [];
mesh2.essentialBoundary = [];
mesh2.LocalNode = [];
mesh2.fs = [];
mesh2.fh = [];

D2 = ones(nRes2^2, 1);
D2(1:(.5*nRes2^2)) = c_hi;
out2 = heat2d_v2(mesh2, D2);


figure
subplot(1,2,1)
imagesc(out1.Tff)
grid off
colorbar
subplot(1,2,2)
imagesc(out2.Tff)
grid off
colorbar























% % analytical solution
% [X, Y] = meshgrid(linspace(0, 1, nRes + 1));
% T_true = 0*X;
% for i = 1:numel(X)
%     x = [X(i), Y(i)];
%     T_true(i) = boundaryTemperature(x);
% end
% 
% figure
% subplot(1,2,1)
% imagesc(out.Tff)
% colorbar
% subplot(1,2,2)
% imagesc(T_true)
% colorbar



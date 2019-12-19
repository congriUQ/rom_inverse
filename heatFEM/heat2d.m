function [Out] = heat2d(mesh, conductivity)
%2D heat conduction main function
%Gives back temperature on point x

%Global stiffness matrix
% Out.globalStiffness = get_glob_stiff2(mesh, Out.diffusionStiffness);
Out.globalStiffness = get_glob_stiff3(mesh.d_glob_stiff_assemble, conductivity, mesh.nEq);


%Finally solving the equation system
Out.naturalTemperatures = Out.globalStiffness\get_glob_force2(mesh, conductivity);

%Temperature field
Out.u(mesh.id, 1) = Out.naturalTemperatures;
Out.u(mesh.essentialNodes, 1) = mesh.essentialTemperatures(mesh.essentialNodes);

end



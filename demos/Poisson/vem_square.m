% Copyright (c) 2019 Adrien Lefieux 
%
% This file is part of mvem.
% 
% mvem is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% mvem is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Foobar.  If not, see <https://www.gnu.org/licenses/>.

% This demo solves the problem : -Laplacian(u)=f, with f=1 on the unit
% square and homogeneuous Dirichlet BC

clear
close all 

degree = 2; % You can check convergence rates
 
% Build Circle with BC
[ Node, Element, Supp] = PolyMesher(@PolySquare, 640, 50);

% Build UVEM polygonal mesh
mesh = Mesh(); % Init
mesh.mesh_from_PolyMesher(Node, Element); 
%clear Node Element % These are not needed anymore
close all

% Build dofs (linear vem)
dofs = Dofs(mesh, degree);

% Operators (Matrices)
operators = Operators(dofs);
operators.build_K();
operators.build_M();

% Home build load (Issue not solved yet: use the L^2 projection)
b = Load();
b.build_constant_load(operators, 1);

% Boundaries
boundaries = Boundaries(operators, b);
boundaries.load_Essential_BCs_connectivity(Supp(:,1)); 
close all

%mesh.build_boundary_alphaShape(0.915); % A bit problematic this boundary story...
% Cooked solution
if degree == 1
  [~, lib] = ismember(Node(Supp(:,1),:), mesh.vertexes,'rows');
  mesh.boundary = unique(lib,'stable');
  for i = 1:length(mesh.boundary)
    exact_sol_BC(i,1) = analytical_solution( ...
      mesh.vertexes(mesh.boundary(i),1), ...
      mesh.vertexes(mesh.boundary(i),2), 100 );
  end
  %boundaries.homogeneous_Essential_BCs(); % Set boundary to zero
  boundaries.nonhomogeneous_Essential_BCs(exact_sol_BC); % Apply given parameters as Dirichlet BC
elseif degree == 2
  % Vertexes BCs
  [~, lib] = ismember(Node(Supp(:,1),:), mesh.vertexes,'rows');
  % Edge BCs
  t = unique(vertcat(dofs.edges_dofs{:}),'stable');
  bc = t(histc(vertcat(dofs.edges_dofs{:}),t)==1); %clear t
  %
  mesh.boundary = union(unique(lib,'stable'),unique(bc,'stable'));
  %boundaries.homogeneous_Essential_BCs();
   mesh.boundary = [unique(lib,'stable') ; unique(bc,'stable')];
  for i = 1:length(unique(lib))
    exact_sol_BC(i,1) = analytical_solution( ...
      mesh.vertexes(mesh.boundary(i),1), ...
      mesh.vertexes(mesh.boundary(i),2), 100 );
  end
  %
  for i = (length(unique(lib))+1):(length(unique(lib))+length(unique(bc)))
    exact_sol_BC(i,1) = analytical_solution( ...
      dofs.edges_coords(mesh.boundary(i)-mesh.nb_vertexes,1), ...
      dofs.edges_coords(mesh.boundary(i)-mesh.nb_vertexes,2), 100 );
  end
  boundaries.nonhomogeneous_Essential_BCs(exact_sol_BC);
end
%boundaries.plot_Boundaries();
pause

% Solve (Preliminary solver system: backslash :) )
x = operators.K\b.b;

% Plot
figure
sol_idx = 1:mesh.nb_vertexes;
mesh.plot_solution(x(sol_idx));

% Plot error
for i = 1:length(mesh.vertexes)
 exact_sol(i,1) = analytical_solution( mesh.vertexes(i,1), mesh.vertexes(i,2), 100 );
end
error = abs(exact_sol - x(sol_idx));

figure
mesh.plot_solution(error);

disp('# Dofs')
disp(mesh.nb_vertexes)

disp('Hmax')
disp(mesh.hmax)

% Compute discrete L^2 error
disp('L^2 error')
error_L2 = sqrt( error'*operators.M(sol_idx,sol_idx)*error );
disp(error_L2)

% Compute discrete sH^1 error
disp('sH^1 error')
error_sH1 = sqrt( error'*operators.K(sol_idx,sol_idx)*error );
disp(error_sH1)

function sol = analytical_solution(x,y,n)

k = 1:2:n;
sol = (1-x^2)/2 - 16/(pi^3)*sum( ( sin(pi*(1+x)/2*k)./(k.^3.*sinh(k*pi) ) )...
  .* ( sinh(pi*(1+y)/2*k) + sinh(pi*(1-y)/2*k ) ) );

end

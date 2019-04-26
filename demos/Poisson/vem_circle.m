% Copyright 2019 Adrien Lefieux
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
% along with mvem.  If not, see <https://www.gnu.org/licenses/>.

% This demo solves the problem : -Laplacian(u)=f, with f=1 on the unit
% circle and homogeneuous Dirichlet BC

clear
close all

degree = 1; % One or two (Two is a patch test)
 
% Build Circle with BC
[ Node, Element, Supp] = PolyMesher(@PolyCircle, 80, 30);

% Build UVEM polygonal mesh
mesh = Mesh(); % Init
mesh.mesh_from_PolyMesher(Node, Element); 
%clear Node Element % These are not needed anymore
%close all
% Build dofs
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
%clear Supp
close all

%mesh.build_boundary_alphaShape(0.915); % A bit problematic this boundary story...
% Cooked solution
if degree == 1 % BC for first order VEM
  [~, lib] = ismember(Node(Supp(:,1),:), mesh.vertexes,'rows');
  mesh.boundary = unique(lib,'stable');

  idx_bc_lin = 1:length(mesh.boundary);
  exact_sol_BC(:,1) = -0.25*( ...
    mesh.vertexes(mesh.boundary(:),1).^2 + ...
    mesh.vertexes(mesh.boundary(:),2).^2 ) + 0.25;

  %boundaries.homogeneous_Essential_BCs(); % set BC to zero
  boundaries.nonhomogeneous_Essential_BCs(exact_sol_BC); % assign given vector
elseif degree == 2 % BC for second order VEM
  % Vertexes BCs
  [~, lib] = ismember(Node(Supp(:,1),:), mesh.vertexes,'rows');
  % Edge BCs
  t = unique(vertcat(dofs.edges_dofs{:}),'stable');
  bc = t(histc(vertcat(dofs.edges_dofs{:}),t)==1); %clear t
  %
  mesh.boundary = union(unique(lib,'stable'),unique(bc,'stable'));
  %boundaries.homogeneous_Essential_BCs(); % set BC to zero
  
  % assign given vector
  mesh.boundary = [unique(lib,'stable') ; unique(bc,'stable')];
  idx_bc_lin = 1:length(unique(lib));
  
  exact_sol_BC(idx_bc_lin) = ...
    -0.25*(mesh.vertexes(mesh.boundary(idx_bc_lin),1).^2 ...
    + mesh.vertexes(mesh.boundary(idx_bc_lin),2).^2) + 0.25;
  % Handle quadratic nodes
  idx_bc_quad = (length(unique(lib))+1):(length(unique(lib))+length(unique(bc)));
  exact_sol_BC( idx_bc_quad ) = -0.25*(...
    dofs.edges_coords(mesh.boundary(idx_bc_quad)-mesh.nb_vertexes,1).^2 + ...
    dofs.edges_coords(mesh.boundary(idx_bc_quad)-mesh.nb_vertexes,2).^2 ) + 0.25;
  
  boundaries.nonhomogeneous_Essential_BCs(exact_sol_BC');
end
pause

% Solve (Preliminary solver system: backslash :) )
x = operators.K\b.b;

% Plot
figure
sol_idx = 1:mesh.nb_vertexes;
mesh.plot_solution(x(sol_idx));

% Plot error
exact_sol = -0.25*(mesh.vertexes(:,1).^2+mesh.vertexes(:,2).^2) + 0.25;
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

%% Conditioning
condK = condest(operators.K);
disp('Conditioning Stiffness:')
disp(condK)

condM = condest(operators.M);
disp('Conditioning Mass:')
disp(condM)

%% Analytical solution 
function u = circle_solution(x,y)

%close all 
u = -0.25*(x.^2+y.^2)+(1/4*(pi/10)^2);

end

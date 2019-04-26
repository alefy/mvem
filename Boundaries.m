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
    
classdef Boundaries
  % The Boundaries class provides all necessary tools for imposing boundary
  % conditions to the BVP problem. 
  %
  %   Boundary supported: 
  %       -> Homogenous Dirichlet BCs
  %       -> nonHomogenous Dirichlet BCs are handled as well, but it is
  %       cumbersome (see demos)
  
  properties
    operators % i.e., the matrices to be modified (handle)
    essential_BC_connectivity
    
    load
  end
  
  methods
    
    % Constructor
    function obj = Boundaries(operators, b)
      
      obj.operators = operators;
      obj.load = b;
    end
    
    function obj = load_Essential_BCs_connectivity(obj, BC_connectivity)
      % This function loads precomputed BCs connectivity matrix. The user
      % is charged to build the essential_BC_connectivity.
      obj.essential_BC_connectivity = BC_connectivity;
    end
    
    function obj = homogeneous_Essential_BCs(obj)
      
      obj.essential_BC_connectivity = obj.operators.dofs.mesh.boundary;
      
      % First stiffness
      obj.operators.K( obj.essential_BC_connectivity, : ) = ...
        sparse( length(obj.essential_BC_connectivity),...
        size( obj.operators.K,1) );
      %
      obj.operators.K( : , obj.essential_BC_connectivity ) = ...
        sparse(size(obj.operators.K,1) ,...
        length(obj.essential_BC_connectivity) );
      
      obj.operators.K(obj.essential_BC_connectivity,...
        obj.essential_BC_connectivity) = ...
        speye( size(obj.essential_BC_connectivity,1),...
        size(obj.essential_BC_connectivity,1));

      % Load
      obj.load.b(obj.essential_BC_connectivity) = ...
        zeros(size(obj.essential_BC_connectivity));
    end
    
    function obj = nonhomogeneous_Essential_BCs(obj, uh)
            
      obj.essential_BC_connectivity = obj.operators.dofs.mesh.boundary;
      
      % Projection
      nonBC = setdiff(1:size(obj.load.b,1),obj.essential_BC_connectivity);
      obj.load.b(nonBC) = obj.load.b(nonBC)...
        -obj.operators.K(nonBC,obj.essential_BC_connectivity)*uh;
      
      % First stiffness (I should build a problem class)
      obj.operators.K( obj.essential_BC_connectivity, : ) = ...
        sparse( length(obj.essential_BC_connectivity),...
        size( obj.operators.K,1) );
      %
      obj.operators.K( : , obj.essential_BC_connectivity ) = ...
        sparse(size(obj.operators.K,1) ,...
        length(obj.essential_BC_connectivity) );
      
      obj.operators.K(obj.essential_BC_connectivity,...
        obj.essential_BC_connectivity) = ...
        speye( size(obj.essential_BC_connectivity,1),...
        size(obj.essential_BC_connectivity,1));
      
                   
      % Load
      obj.load.b(obj.essential_BC_connectivity) = uh;
      
    end
    
    % Plot function
    function plot_Boundaries(obj)
      obj.operators.dofs.mesh.plot_mesh(...
        obj.operators.dofs.mesh.vertexes(...
        obj.operators.dofs.mesh.boundary,:));
    end
    
  end
  
end


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
% along with Foobar.  If not, see <https://www.gnu.org/licenses/>.

classdef Dofs < handle
  %
  % The class Dofs puts altogether local FE Space (VEM, FEM, etc.)
  %
  
  %------------------------------------------------------------------------
  % PROPERTIES
  %------------------------------------------------------------------------
  
  properties
    
    mesh % Of type mesh class
    elements % Actually finite elements
        
    nb_edges_coords = 0;
    nb_edges_dofs = 0;
    edges_dofs % Both for VEM & nodal FEM
    edges_coords % Edge dof coordinates (x,y) in 2D - VEM & nodal FEM
    
    nb_moment_dofs
    
    dofs % All dofs: for matrix assembly

    degree = 1; % We assume uniform order
    
    nb_total_dofs
   
  end
  
  %------------------------------------------------------------------------
  % METHODS
  %------------------------------------------------------------------------
  
  methods
    
    function obj = Dofs(mesh, degree)
      
      if nargin > 1
        obj.degree = degree;
      end
      
      obj.mesh = mesh;
      
      if obj.degree == 1
        for k=1:obj.mesh.nb_cells
          if isa(obj.mesh.cells{k},'Polygon') 
            % Note: this test is not very consistent: not all polygonal
            % methods are VEM. 
            obj.elements{k} = Vem(obj.mesh.cells{k});
            %obj.dofs{k} = obj.mesh.connectivity{k};
          end
        end
        obj.nb_total_dofs = obj.mesh.nb_vertexes;
        
      elseif obj.degree == 2
        
        %obj.edges_dofs = cell(obj.mesh.nb_cells);
        
        for k=1:obj.mesh.nb_cells
          if isa(obj.mesh.cells{k},'Polygon') 
            obj.elements{k} = Vem(obj.mesh.cells{k},2);
            
            % Build connectivity table
            obj.dofs{k} = zeros(obj.elements{k}.nd_total,1);
            
            % Vertices dofs
            elem_nd_V = obj.elements{k}.nd_V;
            obj.dofs{k}(1:elem_nd_V) = obj.mesh.connectivity{k};
            
            % Edges dofs
            elem_nd_E = obj.elements{k}.nd_E;
            obj.edges_dofs{k} = zeros(elem_nd_E,1); % Allocate local memory
            for i = 1:elem_nd_E
              obj.nb_edges_dofs = obj.nb_edges_dofs + 1;
              % Get edge_coord
              edge_coord = obj.elements{k}.edge_points(elem_nd_V+i,:);
              % Add it to the edge table
              position = obj.add_edge_dofs(edge_coord);
              obj.edges_dofs{k}(i) = obj.mesh.nb_vertexes + position;
            end
            obj.dofs{k}( (elem_nd_V+1):(elem_nd_V+elem_nd_E) ) = ...
              obj.edges_dofs{k};
          else
            error('High order polygonal meshes not supported yet.')    
          end

        end
        
        if obj.degree == 2
          for k=1:obj.mesh.nb_cells
            if isa(obj.mesh.cells{k},'Polygon')
              elem_nd_V = obj.elements{k}.nd_V;
              elem_nd_M = obj.elements{k}.nd_M;
              elem_nd_E = obj.elements{k}.nd_E;
              % Virtual dofs (only one for P2)
              virtual_range = ...
                (elem_nd_V+elem_nd_E+1):(elem_nd_V+elem_nd_E+elem_nd_M);
              obj.dofs{k}( virtual_range ) = ...
                obj.mesh.nb_vertexes + obj.nb_edges_coords + k;
            end
          end
          obj.nb_total_dofs = obj.mesh.nb_vertexes + size(obj.edges_coords,1) + ...
        obj.mesh.nb_cells;
        end
      end
      
    end
    
    function dof = get_dofs(obj, k, i)
      % This function returns the i-th dof of element k
      if obj.elements{k}.degree == 1
        dof = obj.mesh.connectivity{k}(i);
      elseif obj.elements{k}.degree == 2
        dof = obj.dofs{k}(i);
      else
        error('High order polygonal elements not supported yet.')
      end
    end
    
    function plot_elem_dofs(obj, l)
      % Plot the mesh with the dofs of element l
      if obj.mesh.dimension == 2
        for k=1:obj.mesh.nb_cells
          h = patch(obj.mesh.vertexes(obj.mesh.connectivity{k},1),...
            obj.mesh.vertexes(obj.mesh.connectivity{k},2),k);
          axis equal; axis tight;
        end
        hold on
        var1 = {'filled', 'r'};
        scatter(obj.mesh.vertexes(obj.mesh.connectivity{l},1),...
                obj.mesh.vertexes(obj.mesh.connectivity{l},2), var1{:});
        var2 = {'d', 'filled' ,'b'};
        
        if obj.degree > 1
          scatter(obj.edges_coords(obj.edges_dofs{l}-obj.mesh.nb_vertexes,1),...
            obj.edges_coords(obj.edges_dofs{l}-obj.mesh.nb_vertexes,2),...
            var2{:});
          uistack(h, 'bottom')
        end
      end
    end
    
  end
  
  methods (Access = private)
    
    function position = add_edge_dofs(obj, edge_coord)
      
      if isempty(obj.edges_coords) % Init
        obj.nb_edges_coords = obj.nb_edges_coords + 1;
        obj.edges_coords(obj.nb_edges_coords,:) = edge_coord;
        position = obj.nb_edges_coords;
        return
      else % Find if edge_vertex has already been inserted.
        position = ismember(obj.edges_coords, edge_coord, 'rows');
        if sum(position) > 0
          % Vertexe has already been inserted
          position = find(position);
          return
        else
          % Add the vertexe:
          obj.nb_edges_coords = obj.nb_edges_coords + 1;
          obj.edges_coords(obj.nb_edges_coords,:) = edge_coord;
          position = obj.nb_edges_coords;
          return
        end
      end
    
  end
  
  end

end
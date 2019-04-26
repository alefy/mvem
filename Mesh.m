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

classdef Mesh < handle
  %
  % Class defining a mesh.
  %
  % The mesh class is mainly a cell container handling the elements.
  %
  
  %------------------------------------------------------------------------
  % PROPERTIES
  %------------------------------------------------------------------------
  
  properties
    
    dimension = 2;
    
    nb_cells = 0;
    nb_vertexes = 0;
    
    vertexes;
    connectivity; % CHECK: cell?
    cells % Will be a cell contaning all elements
    
    boundary
    
    hmax = 0;
    
  end
  
  %------------------------------------------------------------------------
  % METHODS
  %------------------------------------------------------------------------
  
  methods
    
    function obj = add_cell(obj, cell)
      % add_element(element)
      %    element is of what so ever type.
      %    TODO: Define an abstract class element such, type can be checked
      %
      % TODO BAD terminology: should a polygon (if I were to do a 3D version that would add_cell)
    
      if not(isa(cell,'Polygon'))
        
        obj.nb_cells = obj.nb_cells + 1;
        
        nb_vertices = cell.nb_vertices;
        
        poly_connectivity = zeros(nb_vertices,1);
        
        for k=1:nb_vertices
          position = obj.add_vertex(cell.vertices(k,:));
          poly_connectivity(k) = position;
        end
        % Add element in connectivity matrix
        obj.connectivity{obj.nb_cells} = poly_connectivity;
        
        % Add polygon
        obj.cells{obj.nb_cells} = cell;
        
      else
        error('cell not supported')
      end
      
    end
    
    function obj = add_polygon(obj, vertices, edges)
      %
      % edges has to be the "local" edges data structure, which we will
      % modify later.
      
      obj.nb_cells = obj.nb_cells + 1;
      
      nb_vertices = size(vertices,1);
      
      poly_connectivity = zeros(nb_vertices,1);
      
      for k=1:nb_vertices
        position = obj.add_vertex(vertices(k,:));
        poly_connectivity(k) = position;
      %  pos =  ismember(edges, k);
      %  edges(pos) = position;
      end
      
      % Add element in connectivity matrix
      obj.connectivity{obj.nb_cells} = poly_connectivity;
      
      % build polygon
      obj.cells{obj.nb_cells} = Polygon(obj.vertexes(poly_connectivity,:),edges);
    end

    function obj = remove_polygon(obj, vec)
    
      %polygon_to_remove = cells{k};
      
      % Remove polygon
      for k=1:length(vec)
        obj.cells{vec(k)} = 0;
%         obj.cells(vec(k),:) = 0;
        %obj.nb_cells = obj.nb_cells -1;
      end
      %obj.nb_cells = obj.nb_cells -1;
      %obj.nb_cells = obj.nb_cells -1;
      
      % Should remove dofs too
      
      % 
    end
    
    function obj = mesh_from_PolyMesher(obj, v, c)
      
      % Construct a mesh from voronin or PolyMesh
      % v is the vertex list, c is the connectivity table
      %
      
      for k=1:length(c)
        nb_vert = length(c{k});
        % Build edges
        edges = zeros(nb_vert,2);
        edges = [ 1:(nb_vert-1) ; 2:(nb_vert)]';
        edges(end+1,:) = [ nb_vert, 1];
        obj.add_polygon(v(c{k},:), edges); % Issue: copy vertexes data !!!
      end
      
      % Get the boundary
      %obj = build_boundary_alphaShape(obj);
      
      obj = compute_hmax(obj);
    end
    
    function obj = mesh_from_FreeFem(obj, file)
 
      mesh_data = dlmread(file);

      nb_vertexes = uint32(mesh_data(1,1)); % number of vertices
      nb_cells = uint32(mesh_data(1,2)); % number of cell
      %ne = uint32(mesh_data(1,3)); % number of edges

      % Vertices treatment
      vertexes = zeros(nb_vertexes,2);
      
      %vertices_bl = zeros(nv,1); % Boundary treatment

      vertexes = double(mesh_data(2:nb_vertexes+1,1:2));
      
      boundary = find(double(mesh_data(2:nb_vertexes+1,3))); % Boundary treatment

      % Elements treatment
      connectivity = zeros(nb_cells, 3);
       
      %elements_bl = zeros(ns, 1); % Boundary treatment

      connectivity_tmp = double(mesh_data(nb_vertexes...
        +2:nb_vertexes+nb_cells+1,1:3));
      connectivity = mat2cell( connectivity_tmp,...
        ones(size(connectivity_tmp,1),1),3 );
      clear connectivity_tmp
      
      %elements_rl = double(mesh_data(nv+2:nv+ns+1,4)); % Boundary
      %treatment

      % Edges treatment
      %edges = zeros(ne, 3);
      %edges_bl = zeros(ne, 1);

      %edges = double(mesh_data(ns+nv+2:ns+nv+ne+1,1:2));
      %edges_bl = double(mesh_data(ns+nv+2:ns+nv+ne+1,3));

      %clear mesh_data
      
      % Build polygonal structure
      for k=1:length(connectivity)
        %nb_vert = length(connectivity{k});
        nb_vert = 3; % Triangles only
        % Build edges
        edges = zeros(nb_vert,2);
        edges = [ 1:(nb_vert-1) ; 2:(nb_vert)]';
        edges(end+1,:) = [ nb_vert, 1];
        obj.add_polygon(vertexes(connectivity{k},:), edges);
      end
      
      % Get the boundary
      %obj = build_boundary_alphaShape(obj);
      [~, lib] = ismember(vertexes(boundary,:), obj.vertexes,'rows');
      obj.boundary = nonzeros(lib);
      
      
      obj = compute_hmax(obj);
    end
    
    function obj = build_boundary_alphaShape(obj,S)
      %
      % Function to build boundary connectivity using Matlab's alphaShape.
      %
      % Issue: The boundary depends on the shrinking factor S. Problematic,
      % highly user dependant.
      %
      % S shrinking factor
      
      if nargin == 1
        obj.boundary = boundary(obj.vertexes(:,1),obj.vertexes(:,2));
      else
        obj.boundary = boundary(obj.vertexes(:,1),obj.vertexes(:,2),S);
      end
        
    end
    
    function new_bc = immersed_mesh(obj, immersed_mesh, type)
      %
      % Immersed_mesh(obj, immersed_mesh, type)
      %
      % Type : 'inside' 
      % Type : 'outside'
      %
      % Remesh initial mesh with regard to an immersed mesh
      % The immersed mesh is sought, in 2D as a polygon, that is as a
      % Matlab polyshape
      
      if (strcmp(type,'inside'))
        new_bc = immersed_mesh_inside(obj, immersed_mesh);
      elseif (strcmp(type,'outside'))
        new_bc = immersed_mesh_outside(obj, immersed_mesh);
      else
        disp('Unsupported')
      end
      
    end
    
    function obj = compute_hmax(obj)
      for k=1:obj.nb_cells
        hk = obj.cells{k}.diameter;
        if (hk > obj.hmax)
          obj.hmax = hk;
        end
      end
    end
    
    function obj = plot_mesh(obj, points)
      if nargin == 1
        if obj.dimension == 2
          for k=1:obj.nb_cells
            %% old connectivity
           patch(obj.vertexes(obj.connectivity{k},1),...
            obj.vertexes(obj.connectivity{k},2),1);
%              if not(isa(obj.cells{k},'Polygon'))
%                continue;
%              end
%             obj.cells{k}.plot_polyshape()

            axis equal; axis tight;
          end
        end
      
      elseif nargin == 2
        if obj.dimension == 2
          for k=1:obj.nb_cells
            h = patch(obj.vertexes(obj.connectivity{k},1),...
              obj.vertexes(obj.connectivity{k},2),1);
            axis equal; axis tight;
          end
          hold on
          scatter(points(:,1),points(:,2));
          uistack(h, 'bottom')
        end
      end
    end
    
    function obj = plot_solution(obj, x)
      
      if obj.dimension == 2
        for k=1:obj.nb_cells
          patch(obj.vertexes(obj.connectivity{k},1),...
            obj.vertexes(obj.connectivity{k},2),...
            x(obj.connectivity{k}),x(obj.connectivity{k}));
          axis equal; axis tight;
        end
      end
      
    end
    
    
  end % End public method
  
  methods (Access = private)
    
    function position = add_vertex(obj, vertex)
      
      if isempty(obj.vertexes) % Init
        obj.nb_vertexes = obj.nb_vertexes + 1;
        obj.vertexes(obj.nb_vertexes,:) = vertex;
        position = obj.nb_vertexes;
        return
      else % Find if vertex has already been inserted.
%        position = ismember(obj.vertexes, vertex, 'rows');
        position = ismembertol(obj.vertexes, vertex, 'ByRows',true);
        if sum(position) > 0
          % Vertexe has already been inserted
          position = find(position);
          return
        else
          % Add the vertexe:
          obj.nb_vertexes = obj.nb_vertexes + 1;
          obj.vertexes(obj.nb_vertexes,:) = vertex;
          position = obj.nb_vertexes;
          return
        end
      end
      
    end
    
    function new_bc_index = immersed_mesh_inside(obj, immersed_mesh)
      
      %polygons_to_remove;
      nb_polygons_to_remove = 0;
      %polygons_to_add;
      nb_polygons_to_add = 0;
      
      new_bc_vertexes = [];
      
      % Loop over polygon cells
      for k=1:obj.nb_cells
        % Check if intersect immersed_mesh
        intersect = obj.cells{k}.query_polygon_intersection(immersed_mesh);
        
        if intersect
          %% Remesh
          %disp("Intersect")
          [refined_polygon, idxA, idxB] = obj.cells{k}.compute_intersection(immersed_mesh);
          % Compute number of regions
          %disp("Number of regions:")
          %refined_polygon.NumRegions
          % if zero eliminate polygon
          if refined_polygon.NumRegions == 0
            nb_polygons_to_remove = nb_polygons_to_remove + 1;
            polygons_to_remove(nb_polygons_to_remove) = k;
          elseif refined_polygon.NumRegions == 1
              nb_polygons_to_remove = nb_polygons_to_remove + 1;
              polygons_to_remove(nb_polygons_to_remove) = k;
              
              edges = [ (1:(length(refined_polygon.Vertices)-1))', (2:(length(refined_polygon.Vertices)))'];
              edges = [edges ; length(refined_polygon.Vertices) 1];
              
              nb_polygons_to_add = nb_polygons_to_add + 1;
              polygons_to_add_vertex{nb_polygons_to_add} = refined_polygon.Vertices;
              polygons_to_add_edges{nb_polygons_to_add} = edges;
%              obj.add_polygon(refined_polygon.Vertices, edges);

          % Update BC to include interior
          if isempty(new_bc_vertexes)
            new_bc_vertexes = refined_polygon.Vertices(((idxA == 2) | (idxA == 0)),:);
          else
            new_bc_vertexes = [new_bc_vertexes; refined_polygon.Vertices(((idxA == 2) | (idxA == 0)),:)];
          end
          end
          % for 1, 2, ..., j add j polygon
          
        else 
          obj.cells{k}.plot_polyshape();
        end % End intersection test
        
      end
      
      % Eliminate removed polygons
      obj.remove_polygon(polygons_to_remove);
      
      % redo mesh
      cells_temp = obj.cells;
      obj.cells = {};
      nb_cells = obj.nb_cells;
      obj.nb_cells = 0;
      
      obj.nb_vertexes = 0;
      obj.vertexes = [];
      obj.connectivity = {};
      
      for k=1:nb_cells
        if not(isa(cells_temp{k},'Polygon'))
          continue;
        end
        obj.add_polygon(cells_temp{k}.vertices, cells_temp{k}.edges);
      end
      
      % Add new ones
      for i=1:length(polygons_to_add_vertex)
        if ispolycw(polygons_to_add_vertex{i}(:,1),polygons_to_add_vertex{i}(:,2)) 
          % This line should not be necessary, issue within 
          [polyx, polyy] = poly2ccw(polygons_to_add_vertex{i}(:,1),polygons_to_add_vertex{i}(:,2));
          obj.add_polygon([polyx, polyy],polygons_to_add_edges{i});
        else
          obj.add_polygon(polygons_to_add_vertex{i},polygons_to_add_edges{i});
        end
      end
      
      % Update BC
      
      new_bc_bool = ismembertol(obj.vertexes, new_bc_vertexes, 'ByRows',true);
      new_bc_index = find(new_bc_bool==1);
    end % End immersed_mesh_inside
    
    %% NOT WORKING YET
    function new_bc_index = immersed_mesh_outside(obj, immersed_mesh)
      
      %polygons_to_remove;
      nb_polygons_to_remove = 0;
      %polygons_to_add;
      nb_polygons_to_add = 0;
      
      new_bc_vertexes = [];
      
      % Loop over polygon cells
      for k=1:obj.nb_cells
        % Check if intersect immersed_mesh
        intersect = obj.cells{k}.query_polygon_intersection(immersed_mesh);
        
        if intersect
          %% Remesh
          %disp("Intersect")
          [refined_polygon, idxA] = obj.cells{k}.compute_intersection_outside(immersed_mesh);
          % Compute number of regions
          %disp("Number of regions:")
          %refined_polygon.NumRegions
          % if zero eliminate polygon
          %if refined_polygon.NumRegions == 0

          if not(isequal(obj.cells{k}.polygon,refined_polygon))
            % Check number of added Polygons
            if refined_polygon.NumRegions == 1
              nb_polygons_to_remove = nb_polygons_to_remove + 1;
              polygons_to_remove(nb_polygons_to_remove) = k;
              
              edges = [ (1:(length(refined_polygon.Vertices)-1))', (2:(length(refined_polygon.Vertices)))'];
              edges = [edges ; length(refined_polygon.Vertices) 1];
              
              nb_polygons_to_add = nb_polygons_to_add + 1;
              polygons_to_add_vertex{nb_polygons_to_add} = refined_polygon.Vertices;
              polygons_to_add_edges{nb_polygons_to_add} = edges;
              %obj.add_polygon(refined_polygon.Vertices, edges);

              % Update BC to include interior
              if isempty(new_bc_vertexes)
                new_bc_vertexes = refined_polygon.Vertices(((idxA == 2) | (idxA == 0)),:);
              else
                new_bc_vertexes = [new_bc_vertexes; refined_polygon.Vertices(((idxA == 2) | (idxA == 0)),:)];
              end
            elseif (refined_polygon.NumRegions > 1)
              disp("Number of Polygons to add:")
              disp(refined_polygon.NumRegions);
              % Split Polygon connectivity table
              idx = all(isnan(refined_polygon.Vertices),2);
              idy = 1+cumsum(idx);
              idz = 1:size(refined_polygon.Vertices,1);
              C = accumarray(idy(~idx),idz(~idx),[],@(r){refined_polygon.Vertices(r,:)});
              
              for l=1:refined_polygon.NumRegions
                nb_polygons_to_remove = nb_polygons_to_remove + 1;
                polygons_to_remove(nb_polygons_to_remove) = k;
                
                edges = [ (1:(length(C{l})-1))', (2:(length(C{l})))'];
                edges = [edges ; length(C{l}) 1];
                
                nb_polygons_to_add = nb_polygons_to_add + 1;
                polygons_to_add_vertex{nb_polygons_to_add} = C{l};
                polygons_to_add_edges{nb_polygons_to_add} = edges;
              end % Loop over number of added polygons
              if isempty(new_bc_vertexes)
                new_bc_vertexes = refined_polygon.Vertices(((idxA == 2) | (idxA == 0)),:);
              else
                new_bc_vertexes = [new_bc_vertexes; refined_polygon.Vertices(((idxA == 2) | (idxA == 0)),:)];
              end
            end % End check number of region of added polygon
          elseif refined_polygon.NumRegions == 0
              disp("keep")
          end
          % for 1, 2, ..., j add j polygon
         
        elseif not(intersect)
          nb_polygons_to_remove = nb_polygons_to_remove + 1;
          polygons_to_remove(nb_polygons_to_remove) = k;
          obj.cells{k}.plot_polyshape();
        end % End intersection test
      end
      % Eliminate removed polygons
      obj.remove_polygon(polygons_to_remove);
      
      % redo mesh
      cells_temp = obj.cells;
      obj.cells = {};
      nb_cells = obj.nb_cells;
      obj.nb_cells = 0;
      
      obj.nb_vertexes = 0;
      obj.vertexes = [];
      obj.connectivity = {};
      
      for k=1:nb_cells
        if not(isa(cells_temp{k},'Polygon'))
          continue;
        end
        obj.add_polygon(cells_temp{k}.vertices, cells_temp{k}.edges);
      end
      
      % Add new ones
      for i=1:length(polygons_to_add_vertex)
        if ispolycw(polygons_to_add_vertex{i}(:,1),polygons_to_add_vertex{i}(:,2)) 
          % This line should not be necessary, issue within 
          [polyx, polyy] = poly2ccw(polygons_to_add_vertex{i}(:,1),polygons_to_add_vertex{i}(:,2));
          obj.add_polygon([polyx, polyy],polygons_to_add_edges{i});
        else
          obj.add_polygon(polygons_to_add_vertex{i},polygons_to_add_edges{i});
        end
      end
      
      % Find BC indexes
      new_bc_bool = ismembertol(obj.vertexes, new_bc_vertexes, 'ByRows',true);
      new_bc_index = find(new_bc_bool==1);
    end
    
  end % End private method
  
end % End class mesh
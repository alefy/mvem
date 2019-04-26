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

classdef Polygon
  %
  % Class defining a polygon. The polygon shall not selfintersect.
  %
  
  %------------------------------------------------------------------------
  % PROPERTIES
  %------------------------------------------------------------------------
  
  properties (SetAccess = public)
    dimension = 2;
    
    % Polygon definition
    vertices
    edges
    
    %Some data
    nb_vertices
    nb_edges
    
    % 
    area
    centroid
    diameter
    
    % number of holes
    nb_holes = 0;
    holes_nb_vertices;
    
    % Parameters regarding polygon tesselation
    tesselated = false;
    tesselation % Object of delaunayTriangulation class
    tesselation_areas
    
    plotting = 1
    
    polygon % private access to Matlab polyshape structure
  end
  
  properties (Access = private)
    index_loop
    
  end
   
  %------------------------------------------------------------------------
  % METHODS
  %------------------------------------------------------------------------
  
  methods (Access = public)
    
    function obj = Polygon(vertices, edges)
      %
      % Polygon constructor.
      % Input:
      %   - vertices ([x_1,y_1 ; x_2, y_2 ; ... ; x_n,y_n])
      %   - connectivity table (Convention: counter clockwise)
      %
      
      obj.nb_vertices = size(vertices,1);
      obj.nb_edges = size(edges,1);
      
      obj.vertices = zeros(size(vertices));
      obj.edges = zeros(size(edges));
      
      obj.vertices = vertices;
      obj.edges = edges;
      
      obj.polygon = polyshape(vertices(edges(:,1), 1),vertices(edges(:,1), 2),'SolidBoundaryOrientation','ccw');
      
      obj.index_loop = [1:obj.nb_vertices,1];
      
      obj = compute_Area(obj); 
      % Check orientation
      if obj.area < 0
        warning( [ 'The orientation is not counter clockwise: !' ...
                   'Perform reorientation' ]);
        %obj.vertices = obj.vertices(end:-1:1,:);
        %obj.area = -obj.area;
      end
      obj = compute_Centroid(obj);
      obj = compute_Diameter(obj);
    end
    
    function n = get_Normal(obj,k)
      % This function return the normal of facet (edge in 2D ) k
      if obj.dimension == 2
        x1 = obj.vertices( obj.edges(k,1), 1 );
        x2 = obj.vertices( obj.edges(k,2), 1 );
        y1 = obj.vertices( obj.edges(k,1), 2 );
        y2 = obj.vertices( obj.edges(k,2), 2 );
        n = [ ( y2-y1 ) ; -( x2 - x1 ) ] / sqrt( (x2-x1)^2 + (y2-y1)^2 );
      else
        error('Normal not defined')
      end
    end
    
    function plot_polygon(obj, points, var)
      
      if nargin == 1
        if obj.nb_holes < 1
          patch(obj.vertices(:,1),obj.vertices(:,2), ...
            zeros(size(obj.vertices(:,2))));
        else
          for p = 1:(obj.nb_holes+1)
            patch(obj.vertices(obj.holes_nb_vertices{p},1),...
                  obj.vertices(obj.holes_nb_vertices{p},2), ...
                  0);
            axis equal; axis tight;
            hold on
          end
        end
      axis equal; axis tight;
      elseif nargin == 2
        h = patch(obj.vertices(:,1),obj.vertices(:,2), ...
          zeros(size(obj.vertices(:,2))));
        axis equal; axis tight;
        hold on
        scatter(points(:,1),points(:,2));
        uistack(h, 'bottom')
      else
        h = patch(obj.vertices(:,1),obj.vertices(:,2), ...
          zeros(size(obj.vertices(:,2))));
        axis equal; axis tight;
        hold on
        scatter(points(:,1),points(:,2), var{:});
        uistack(h, 'bottom')
      end
    end
    
    function plot_polyshape(obj)
      %figure
      plot(obj.polygon)
      hold on
      axis equal
    end
    
    function obj = add_hole(obj, vertices, edges)
      
      obj.nb_holes = obj.nb_holes + 1;
      
      % Build temporary polygon
      hole = Polygon(vertices, edges);
      % edges(end:-1:1) To avoid automatic reorientation 
      % (we assume clockwise order for hole)
      
      
      if obj.nb_holes == 1
        obj.holes_nb_vertices{1} = 1:size(obj.vertices,1);
      end
      obj.holes_nb_vertices{obj.nb_holes+1} = ...
        (length(obj.holes_nb_vertices{obj.nb_holes})+1):...
        (length(obj.holes_nb_vertices{obj.nb_holes})+size(vertices,1));
      
      obj.edges = [ obj.edges ; size(obj.vertices,1) + edges];
      obj.vertices = [ obj.vertices ; hole.vertices ];
      
      obj.nb_vertices = obj.nb_vertices + size(vertices,1);
      obj.nb_edges = obj.nb_edges + size(edges,1);
      clear vertices edges
      
      obj.centroid = (obj.centroid*obj.area-hole.centroid*hole.area) ...
        / (hole.area + obj.area);
      obj.area = obj.area + hole.area;
      
    end
    
    function obj = tesselate(obj)
      
      obj = obj.compute_Delaunay();
      obj.tesselated = true;
      
    end
    
    function polygons_list = refine(obj,n)
      
      polygons_list = refine_by_elements(obj,n);
      
    end
    
    function test_query = query_polygon_intersection(obj, polygon)
      % Function to query if this polygon intersect another polygon
      test_query = overlaps(obj.polygon, polygon);
    end
    
    function [refined_polygon, idxA, idxB] = compute_intersection(obj, polygon)
      % Function to query if this polygon intersect another polygon
      [refined_polygon, idxA, idxB] = subtract(obj.polygon, polygon);
      %disp("Number of region: ")
      %refined_polygon.NumRegions
      %plot(refined_polygon)
    end
    
    function [refined_polygon, idxA, idxB] = compute_intersection_outside(obj, polygon)
      % Function to query if this polygon intersect another polygon
      [refined_polygon, idxA, idxB] = intersect(obj.polygon, polygon);
      %disp("Number of region: ")
      %refined_polygon.NumRegions
      %plot(refined_polygon)
    end
    
  end
  
  methods (Access = private)
    
    function obj = compute_Area(obj)
      %
      % Polyarea: matlab built in computationa geometry function.
      %
      
      % Issue: always positive
      %obj.area = polyarea(obj.vertices(:,1),obj.vertices(:,2));
      
      % Use this construction to check correct orientation
       obj.area = 0.5 * sum( obj.vertices(obj.index_loop(1:end-1),1) .* ...
                             obj.vertices(obj.index_loop(2:end),2) - ...
                             obj.vertices(obj.index_loop(2:end),1) .* ...
                             obj.vertices(obj.index_loop(1:end-1),2) );
%        obj.area = area(obj.polygon); % New from Matlab polyshape
    end
    
    function obj = compute_Centroid(obj)
      %
      % c_x = frac{1}{6\cdot area}\sum{i=1}^{N}(x_i+x_{i+1})(x_iy_{x_i} 
      %       - x_{i+1}y_i)
      % c_y = frac{1}{6\cdot area}\sum{i=1}^{N}(y_i+y_{i+1})(x_iy_{x_i} 
      %       - x_{i+1}y_i)
      %
      % Cf. Paul Bourke
      
      C = 1/(6*obj.area);
      
      c_x = C*sum( ( obj.vertices(obj.index_loop(1:end-1),1) + ...
                     obj.vertices(obj.index_loop(2:end),1) ) .* ...
                   ( obj.vertices(obj.index_loop(1:end-1),1) .* ...
                     obj.vertices(obj.index_loop(2:end),2) - ...
                     obj.vertices(obj.index_loop(2:end),1) .* ...
                     obj.vertices(obj.index_loop(1:end-1),2) ) );
           
      c_y = C*sum( ( obj.vertices(obj.index_loop(1:end-1),2) + ...
                     obj.vertices(obj.index_loop(2:end),2) ) .* ...
                   ( obj.vertices(obj.index_loop(1:end-1),1) .* ...
                     obj.vertices(obj.index_loop(2:end),2) - ...
                     obj.vertices(obj.index_loop(2:end),1) .* ...
                     obj.vertices(obj.index_loop(1:end-1),2) ) );
      
      obj.centroid = [c_x,c_y];
      
      % New version
     % obj.centroid = centroid(obj.polygon);

    end
    
    function obj = compute_Diameter(obj)
      %
      % This function computes the diameter of the polygon.
      % Simple home made algorithm, figure out how efficient !
      %
      
      % Algorithm
      % First we take the convexhull of the vertices of the polygon
      K = convhull(obj.vertices(:,1),obj.vertices(:,2),'simplify',true);
      
      points = obj.vertices(K(1:end-1),:); % the last point is the first
                                           % thus eliminate it
      
      % Loop over the convexhull points
      obj.diameter = 0;
      nb_convexhull_points = size(points,1);
      for i=1:nb_convexhull_points
        % compute max distance with other points
        dmax = max( sqrt( (points(2:end,1) - points(1,1)).^2 + ...
                          (points(2:end,2) - points(1,2)).^2) );
         if dmax > obj.diameter
           obj.diameter = dmax;
         end
         % Eliminate first point
         points = points(2:end,:);
                          
      end
      
    end
  
    function obj = compute_Delaunay(obj)
      
      delaunay = ...
        delaunayTriangulation(obj.vertices,obj.edges);
      
      inout = isInterior(delaunay);
      
      % Add support for holes:
      obj.tesselation = triangulation(delaunay.ConnectivityList(inout,:), ...
        obj.vertices);
      
      obj.tesselation_areas = ...
        obj.compute_area(obj.tesselation.Points, ...
                         obj.tesselation.ConnectivityList);
                       
    end
    
    function area = compute_area(obj, vertexes, connectivity)
      %
      % Compute the area of the triangle k of the local triangulation
      %
      % Could be added as a static method in the class triangle

      b2 = vertexes(connectivity(:,3),2)-vertexes(connectivity(:,1),2 );
      b3 = vertexes(connectivity(:,1),2)-vertexes(connectivity(:,2),2);

      % c1 = x3-x2
      c2 = vertexes(connectivity(:,1),1)-vertexes(connectivity(:,3),1);
      c3 = vertexes(connectivity(:,2),1)-vertexes(connectivity(:,1),1);

      area = -0.5*(-c3.*b2+b3.*c2);

    end
    
    function poly_list = refine_by_internal_edges(obj,n)
      
      % Get number of elements in tesselation
      %nb_elements_tesselation = size(obj.tesselation.ConnectivityList,1);
      
      % Get internal edges
      edges_tes = obj.tesselation.edges();
%       edges_tes = unique(...
%               [obj.tesselation.ConnectivityList(:,[2,3]); ...
%                obj.tesselation.ConnectivityList(:,[3,1]); ...
%                obj.tesselation.ConnectivityList(:,[1,2])],'rows','legacy');
      toexclude = obj.tesselation.freeBoundary();
      internal_edges = setdiff( edges_tes,...
                                [ toexclude(:,1:2) ; toexclude(:,2:-1:1)], ...
                                'rows','legacy');
%      internal_edges = sort(internal_edges,2);            
%       
      nb_edges_tesselation = size(internal_edges,1);
                              
      nb_sub_edges = floor(nb_edges_tesselation/n);
%       if ~isinteger(nb_sub_edges)
%         n = n - 1;
%         nb_sub_edges = nb_edges_tesselation/n; 
%       end
      

            % Init build patchs
      num_patch = 1;
      num_elem_in_patch = 1;
      num_picked_elements = 1;
                              
      already_picked_elements = zeros(nb_edges_tesselation,1);
      
      idx_1 = 1:nb_sub_edges:nb_edges_tesselation;
      
      elem_local_list = [];
      
      for i=1:(size(idx_1,2)-1)
        
        idx = idx_1(i):(idx_1(i)+nb_sub_edges-1);
        
        neighbors_tes = obj.tesselation.edgeAttachments(...
                                        internal_edges(idx,1),...
                                        internal_edges(idx,2));
        % Init
        if i==1
          elem_local_list{num_patch} = neighbors_tes{1};
          already_picked_elements = neighbors_tes{1};
          if size(neighbors_tes,1) > 1
            for j=2:size(neighbors_tes,1)
              toadd = setdiff(neighbors_tes{j},already_picked_elements);
              elem_local_list{num_patch} = union(elem_local_list{num_patch},toadd);
              already_picked_elements = union(already_picked_elements,neighbors_tes{j});
            end
          end
        % Loop
        else
          elem_local_list{num_patch} = setdiff(neighbors_tes{1},already_picked_elements);
          already_picked_elements = union(already_picked_elements,elem_local_list{num_patch});
          if size(neighbors_tes,1) > 2
            for j=2:size(neighbors_tes,1)
              toadd = setdiff(neighbors_tes{j},already_picked_elements);
              elem_local_list{num_patch} = union(elem_local_list{num_patch},toadd);
              already_picked_elements = union(already_picked_elements,neighbors_tes{j});
            end
          end
        end
        
        
        num_patch = num_patch + 1;
        
   
      end
      % Final patch
      idx = (idx_1(i)+nb_sub_edges):(nb_edges_tesselation);
      if ~isempty(idx) 
        neighbors_tes = obj.tesselation.edgeAttachments(...
        internal_edges(idx,1),...
        internal_edges(idx,2));
      
      elem_local_list{num_patch} = setdiff(neighbors_tes{1},already_picked_elements);
      already_picked_elements = union(already_picked_elements,elem_local_list{num_patch});
      if size(neighbors_tes,1)>1
        for j=2:size(neighbors_tes,1)
          toadd = setdiff(neighbors_tes{j},already_picked_elements);
          elem_local_list{num_patch} = union(elem_local_list{num_patch},toadd);
        end
      end
      end
      % Build polygon edges

      
      %% Export
      for k=1:size(elem_local_list,2)

        tri_patch = triangulation(...
          obj.tesselation.ConnectivityList(elem_local_list{k},:), ...
          obj.tesselation.Points(:,1),...
          obj.tesselation.Points(:,2));
        
        [fe, pts ] = freeBoundary(tri_patch);
        
        poly_list{k} = Polygon(pts,fe);
        
        if obj.plotting
         %triplot(obj.tesselation)
         %hold on
         patch(pts(fe',1), pts(fe',2), k);
         
        end
        
      end
      
    end
    
    function poly_list = refine_by_elements(obj,n)
      %
      % The refine function returns a list of polygons
      % (including its tesseletion, which we get for free at this stage)
      %
      % WARNING: Works only for enclosed boundaries
      %
      
      % Get number of elements in tesselation
      nb_elements_tesselation = size(obj.tesselation.ConnectivityList,1);
      
      % Compute pattern: number of sub elements by new polygons
      nb_sub_elements = ceil(nb_elements_tesselation/n);
      
  %    if ~isinteger(nb_elements_tesselation/n)
  %      division_not_integer = true;
  %    end
      
      idx_1 = 1:nb_sub_elements:nb_elements_tesselation;
      
      already_picked_elements = zeros(nb_elements_tesselation,1);
      
      % Init build patchs
      num_patch = 1;
      num_elem_in_patch = 1;
      num_picked_elements = 1;
      neighbors_tes = obj.tesselation.neighbors(1);
      
      % Build patch: by looping over elements
      for i = 1:(nb_elements_tesselation+1); % Why +1?
      %while already_picked_elements(end) ~= 0
      % Pick an element
        
        for j = 1:length(neighbors_tes)
          if isnan(neighbors_tes(j))
            continue
          elseif sum(ismember(already_picked_elements,neighbors_tes(j))) == 0
            picked_element = neighbors_tes(j);
            elem_local_list{num_patch}(num_elem_in_patch) = picked_element;
            
            already_picked_elements(num_picked_elements) = picked_element;
            num_elem_in_patch = num_elem_in_patch + 1;
            num_picked_elements = num_picked_elements + 1;
            break
          end
        end
        

        % First try for enclosed elements:
        if sum(isnan(neighbors_tes)) == 2
          for j = 1:length(neighbors_tes)
            if isnan(neighbors_tes(j))
              continue
            else
              picked_element = neighbors_tes(j);
              break
            end
          end
        end
        
        neighbors_tes = obj.tesselation.neighbors(picked_element);
         
        if (num_elem_in_patch > nb_sub_elements) && ...
            (num_picked_elements <= (nb_sub_elements*n) )
          num_patch = num_patch + 1;
          num_elem_in_patch = 1;
        end
       

        %i = i + 1;
      end % End loop
      
      % Check if some elements were missed:
      nb_missed_elements = size(...
        setdiff( 1:nb_elements_tesselation,...
                 already_picked_elements,'stable'),2 );
      while nb_missed_elements > 0
        idx_missed = setdiff( 1:nb_elements_tesselation,...
                 already_picked_elements,'stable');
        %idx_missed = find(bool_missed);
        for i = 1:nb_missed_elements
          
          missed_elem = idx_missed(i);
          % Find its neighbors
          neighbors_tes = obj.tesselation.neighbors(missed_elem);
          % Find in which patch this element is:          
          missed_patch = max( find(find(ismember(already_picked_elements,...
                              neighbors_tes(~isnan(neighbors_tes))) ) >= idx_1));
                            
          % Add element to patch
          elem_local_list{missed_patch}(end+1) = missed_elem;
          already_picked_elements(num_picked_elements) = missed_elem;
          num_picked_elements = num_picked_elements + 1;
        end
        nb_missed_elements = size(...
            setdiff( 1:nb_elements_tesselation,...
                 already_picked_elements,'stable'),2 );
      end
      
      for k=1:size(elem_local_list,2)
        
        %idx = idx_1(k):(idx_1(k)+nb_sub_elements-1);
        
        tri_patch = triangulation(...
          obj.tesselation.ConnectivityList(elem_local_list{k},:), ...
          obj.tesselation.Points(:,1),...
          obj.tesselation.Points(:,2));
  
        [fe, pts ] = freeBoundary(tri_patch);
  
        poly_list{k} = Polygon(pts,fe);
        
        
        if obj.plotting
          patch(pts(:,1), pts(:,2), k);
          axis tight; axis equal;
        end
      end
    end
  end
    
end
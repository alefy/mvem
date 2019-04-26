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

classdef Vem
  %
  % Class Vem defines local matrices.
  %
  
  %------------------------------------------------------------------------
  % PROPERTIES
  %------------------------------------------------------------------------
  
  properties
    dimension = 2;
    
    % Dofs definitions: nd : number of dofs
    nd_V % Vertices dofs (equal to the number of vertices)
    nd_E % Edge dofs (equal to k-1 internal points of the k+1 
         % Gauss-Lobatto quadrature on edge)
    nd_M % Virtual dofs 
    
    nd_total % dim(FESpace) i.e. nd_V + nd_E + nd_M
    
    degree = 1;
    
    % dim( P_k(E) ) = dim( M_k(E) ) = (k+1)(k+2)/2 (in 2D !!! )
    nk % Dimension of the polygonal space
    
    polygon % Associated local element
    
    basis_type = 'Monomial';
    
    % Matrices for Pi^grad
    G % nk times nk
    D % nk times nd_total
    B % nk times nd_total 
    H % nk times nk
 
  end
  
  properties (Access = public)
    monomial_index;
    laplacian_monomial_index;    
    
    quadrature_segment
    quadrature_triangle
    
    local_dofs
    edge_points
    edge_dofs
    
    local_triangulation
    areas_local_triangulation
  end
  
  %------------------------------------------------------------------------
  % METHODS
  %------------------------------------------------------------------------
  
  methods
    
    function obj = Vem(polygon, k)
      %
      % Build first order VEM 
      
      if nargin > 1
        obj.degree = k;
      end
      
      obj.polygon = polygon;
      
      % Build local number of dofs
      obj.nd_V = obj.polygon.nb_vertices;
      obj.nd_E = obj.polygon.nb_vertices*(obj.degree-1);
      obj.nd_M = obj.dimension_polynomial_space(obj.degree-2);
      obj.nd_total = obj.nd_V + obj.nd_E + obj.nd_M;
      
      % Or: nd_V*obj.degree + ((obj.degree - 1 ) obj.degree) / 2 in 2D
      obj.nk = obj.dimension_polynomial_space(obj.degree); 
      
      obj = build_monomial_index(obj);
      
      if obj.degree > 1
        
        obj = build_laplacian_monomial_index(obj);
        obj = build_high_order_dofs(obj);
        % Build local dofs for high order
        
      end
      
      obj = build_D(obj);
      obj = build_G(obj);
      obj = build_B(obj);
      obj = build_H(obj); % Extend to k=1
    end
    
  end
  
  methods (Access = public) % Public in debug mode, private otherwise
    
%     function y = basis(obj,X,k)
%       %
%       % For the VEM the monomial one;
%       y = obj.monomial_basis(X,k);
%     end
    
    function obj = build_monomial_index(obj)
      [X,Y] = meshgrid(0:obj.degree,0:obj.degree);
      tmp = [X(:) Y(:)];
      total_degree = sum(tmp,2);
      [~,idx] = sort(total_degree);
      tmp = tmp(idx,:);
      total_degree = sum(tmp,2);
      obj.monomial_index = tmp(total_degree < (obj.degree+1),:);
    end
    
    function obj = build_laplacian_monomial_index(obj)
      [X,Y] = meshgrid(0:obj.degree-2,0:obj.degree-2);
      tmp = [X(:) Y(:)];
      total_degree = sum(tmp,2);
      [~,idx] = sort(total_degree);
      tmp = tmp(idx,:);
      total_degree = sum(tmp,2);
      obj.laplacian_monomial_index = tmp(total_degree < (obj.degree-1),:);
    end
    
    function y = monomial_basis(obj,X,k)
      %
      % Return the monomial basis value at point X for base k
      y = prod( ( (X-obj.polygon.centroid)/obj.polygon.diameter ) .^ ...
        obj.monomial_index(k,:) );
    end
    
    function Dy = grad_monomial_basis(obj,X,k)
       %
       % Return the gradient of monomial basis value at point X for base k
       Dy = zeros(obj.dimension,1);
       
       for i=1:obj.dimension
         if ( obj.monomial_index(k,i) == 0 )
           Dy(i,1) = 0;
         else
            alpha = obj.monomial_index(k,:);
            alpha(i) = alpha(i) - 1; 
            Dy(i,1) = obj.monomial_index(k,i)/obj.polygon.diameter*...
              prod( ( (X-obj.polygon.centroid)/obj.polygon.diameter ) .^ alpha );   
         end
       end
       
    end
    
    function y = laplacian_monomial_basis(obj,X,k)
      %
      % Return the monomial basis value at point X for base k
      y = prod( ( (X-obj.polygon.centroid)/obj.polygon.diameter ) .^ ...
        obj.laplacian_monomial_index(k,:) );
    end
    
    function P0m = P0m(obj,i)
       %
       % Note that the first column is always equal to one.
       %
       % Fix constant for the monomial basis.
       P0m = 0;
       if obj.degree == 1
         for j = 1:obj.nd_V
           P0m = P0m + obj.monomial_basis( obj.polygon.vertices(j,:), i );
         end
         P0m = P0m/obj.polygon.nb_vertices;
       elseif ( obj.degree > 1 )% && ( obj.degree < 3)
         % On the limitation of the degree: the highest degree integration
         % precision I have for triangles is 6 !!! 
         % 
         % Actually for the case k=2 integration over the whole element
         % could be avoided quite easily since after integration by parts 
         % the laplacian term is a constant, thus we need to find a
         % primitive of a polynome.
         %
         % Loop over triangles
         for k=1:obj.polygon.tesselation.size(1) % Could be vectorized
           for g=1:obj.quadrature_triangle.nb_points
           
             pt = barycentricToCartesian( obj.polygon.tesselation, k, ...
               obj.quadrature_triangle.points(g,:) );
           
             P0m = P0m + obj.polygon.tesselation_areas(k).*...
               obj.quadrature_triangle.weights(g).* ...
               obj.monomial_basis( pt , i );
           
           end
         end
         P0m = P0m/obj.polygon.area;
       else
         error('Construction of P0m: issue.')
       end
       
    end
     
    function obj = build_G(obj)
       %
       % Build the local G matrix
       obj.G = zeros( obj.nk , obj.nk );
       
       % Gradient terms
       if obj.degree == 1
         % First line
         for j=1:obj.nk
           obj.G(1,j) = obj.P0m(j);
         end
         % No need to integrate over the domain in case of k=1
         for i=2:obj.nk
           for j=2:obj.nk
             for k=1:obj.polygon.nb_edges
               v1 = obj.polygon.vertices(obj.polygon.edges(k,1),:); % Vertex 1
               v2 = obj.polygon.vertices(obj.polygon.edges(k,2),:); % Vertex 2
               normal = obj.polygon.get_Normal(k); % Normal of the edge k
               
               distance = sqrt( (v1(1)-v2(1))^2 + (v1(2)-v2(2))^2 );
               
               % Use trapezoium rule
               grad_m_i_v1 = obj.grad_monomial_basis(v1,i);
               m_j_v1 = obj.monomial_basis(v1,j);
               
               grad_m_i_v2 = obj.grad_monomial_basis(v2,i);
               m_j_v2 = obj.monomial_basis(v2,j);
               
               obj.G(i,j) = obj.G(i,j) + distance/2 * ( ...
                 sum( grad_m_i_v1.*normal )*m_j_v1 + ...
                 sum( grad_m_i_v2.*normal )*m_j_v2 );
             end
           end
         end
       elseif ( obj.degree > 1 )
         % First line
         for j=1:obj.nk
           obj.G(1,j) = obj.P0m(j);
         end
         % Note that I could use a lower quadrature rule for this case.
         for i=2:obj.nk
           for j=2:obj.nk
             for k=1:obj.polygon.tesselation.size(1) % Could be vectorized
               for g=1:obj.quadrature_triangle.nb_points
           
                 pt = barycentricToCartesian( obj.polygon.tesselation, k, ...
                   obj.quadrature_triangle.points(g,:) );
           
                 obj.G(i,j) = obj.G(i,j) + obj.polygon.tesselation_areas(k).*...
                   obj.quadrature_triangle.weights(g).* ...
                   ( dot( obj.grad_monomial_basis(pt,i), ...
                          obj.grad_monomial_basis(pt,j) ) );
                 
               end
             end
           end
         end
       else
         error('Build G issue')
       end
       
    end
     
    function obj = build_D(obj)
       % Build the local D matrix
       obj.D = zeros( obj.nd_total, obj.nk );
       
       if (obj.degree == 1)
         for i=1:(obj.nd_V)
           for j=1:obj.nk
             obj.D(i,j) = obj.monomial_basis( obj.polygon.vertices(i,:), j);
           end
         end
       elseif (obj.degree >= 1)
         % Build for the Vextexes dofs
         for i=1:(obj.nd_V)
           for j=1:obj.nk
             obj.D(i,j) = obj.monomial_basis( obj.polygon.vertices(i,:), j);
           end
         end
         % Build for the Edges dofs
         for i=(obj.nd_V+1):(obj.nd_V+obj.nd_E)
           for j=1:obj.nk
             obj.D(i,j) = obj.monomial_basis( ...
               obj.edge_points(i,:), j);
           end
         end
         % Build for the Internal dofs
         for i=(obj.nd_V+obj.nd_E+1):(obj.nd_V+obj.nd_E+obj.nd_M)
           for j=1:obj.nk
             for k=1:obj.polygon.tesselation.size(1) % Could be vectorized
               for g=1:obj.quadrature_triangle.nb_points
           
                 pt = barycentricToCartesian( obj.polygon.tesselation, k, ...
                   obj.quadrature_triangle.points(g,:) );
           
                 obj.D(i,j) = obj.D(i,j) + obj.polygon.tesselation_areas(k).*...
                   obj.quadrature_triangle.weights(g).* ...
                   ( obj.monomial_basis( pt , j ) .* ...
                     obj.laplacian_monomial_basis( pt , ...
                     i - (obj.nd_V+obj.nd_E) ) ) / obj.polygon.area;
               end
             end
           end
         end
         
       end
    end
    
    function obj = build_B(obj)
       % Build the local B matrix
       obj.B = zeros( obj.nk, obj.nd_total);
       
       % First line
       if (obj.degree == 1)
         for j=1:obj.nd_total
           obj.B(1,j) = 1/obj.nd_V;
         end
       
         for i=2:obj.nk
           for j=1:obj.nd_total
             for k=1:obj.polygon.nb_edges
               v1 = obj.polygon.vertices(obj.polygon.edges(k,1),:); % Vertex 1
               v2 = obj.polygon.vertices(obj.polygon.edges(k,2),:); % Vertex 2
               normal = obj.polygon.get_Normal(k); % Normal of the edge k
               
               distance = sqrt( (v1(1)-v2(1))^2 + (v1(2)-v2(2))^2 );
               
               % Use trapezoium rule
               grad_m_i_v1 = obj.grad_monomial_basis(v1,i);
               %phi_j_v1 = obj.phi_i(j,v1);
               phi_j_v1 = delta(j,obj.polygon.edges(k,1));
               
               grad_m_i_v2 = obj.grad_monomial_basis(v2,i);
               %phi_j_v2 = obj.phi_i(j,v2);
               phi_j_v2 = delta(j,obj.polygon.edges(k,2));
     
               obj.B(i,j) = obj.B(i,j) + distance/2 * ( ...
                 sum( grad_m_i_v1.*normal ) * phi_j_v1 + ...
                 sum( grad_m_i_v2.*normal ) * phi_j_v2 );
              
             end
           end
         end
         
       else
         % First line
         for j=(obj.nd_V*obj.degree+1):obj.nd_total
           obj.B(1,j) = obj.P0m( j-(obj.nd_V*obj.degree) );
         end
         
         for i=2:obj.nk
           d_alpha_beta = obj.compute_d_alpha_beta(i);
           for j=1:obj.nd_total
             % First compute the components owning the Laplacian
             if j > (obj.nd_V*obj.degree)
             obj.B(i,j) = obj.B(i,j) - ...
               obj.polygon.area * ( d_alpha_beta(j - (obj.nd_V*obj.degree) ) );
             end
             % The boundary part
             % Loop over the edges
             for k = 1:obj.polygon.nb_edges
               v1 = obj.polygon.vertices(obj.polygon.edges(k,1),:); % Vertex 1
               v2 = obj.polygon.vertices(obj.polygon.edges(k,2),:); % Vertex 2
               normal = obj.polygon.get_Normal(k); % Normal of the edge k
             
               distance = sqrt( (v1(1)-v2(1))^2 + (v1(2)-v2(2))^2 );
             
               % Loop over Gauss-Lobatto points
               for g=1:obj.quadrature_segment.nb_points
                 grad_m_i = obj.grad_monomial_basis(...
                   obj.edge_points(obj.edge_dofs(k,g),:),i);
                 phi_j = delta(j,obj.edge_dofs(k,g));
                 
                 % Case monomial (internal) dofs
                 if j > (obj.nd_V*obj.degree+1)
                   phi_j = 0;
                 end
                 
                 obj.B(i,j) = obj.B(i,j) + ...
                   distance * ( obj.quadrature_segment.weights(g)* sum( grad_m_i.*normal ) * phi_j ); % NOT FINISH
               end
            end
           end
         end
       end
    end
    
    function obj = build_H(obj)
      obj.H = zeros( obj.nk , obj.nk );
      
      if obj.degree == 1
        
        obj.polygon = obj.polygon.tesselate(); % TODO: Check copying obj?
        
        % Build quadrature rule for triangles
        obj.quadrature_triangle = ...
          Quadrature('Gauss-Legendre', 'triangle', 3 );
        

      end
        
      for i=1:obj.nk
        for j=1:obj.nk
          for k=1:obj.polygon.tesselation.size(1) % Could be vectorized
            for g=1:obj.quadrature_triangle.nb_points
              
              pt = barycentricToCartesian( obj.polygon.tesselation, k, ...
                obj.quadrature_triangle.points(g,:) );
              
              obj.H(i,j) = obj.H(i,j) + obj.polygon.tesselation_areas(k).*...
                obj.quadrature_triangle.weights(g).* ...
                ( obj.monomial_basis(pt,i)* obj.monomial_basis(pt,j) );
              
            end
          end
        end
      end
      
    end
    
    function r = phi_tmp_i(obj,i,X)
      %
      % Temporary function
      %
      % Degrees of freedom of Phi_i\in V_k(E):
      % i.e.,: dof_i(X_j) = delta_{ij}
      if obj.polygon.vertices(i,:) == X
        r = 1;
        return
      else
        r = 0;
        return
      end
    end
    
    function r = phi_i(obj,i,X)
       %
       % Degrees of freedom of Phi_i\in V_k(E):
       % i.e.,: dof_i(X_j) = delta_{ij}
       if obj.polygon.vertices(i,:) == X
         r = 1;
         return
       else
         r = 0;
         return 
       end
    end
     
    function obj = build_high_order_dofs(obj)
      %
      % This function constructs the edges dofs
      %
      %
      
      obj.polygon = obj.polygon.tesselate();
      
      % Build Gauss points
      obj.quadrature_segment = ...
        Quadrature('Gauss-Lobatto', 'segment', obj.degree+1);
      
      % First build edge nodes
      
      % Attention: edge_points is a tensor:
      %  ( Edges x Number of points x spatial_dimension)
      %obj.edge_points = zeros( obj.polygon.nb_edges, ...
      %  obj.dimension ); 
      
      % For simplicity put together vertices and edge points
      obj.edge_points = zeros( obj.nd_V*obj.degree , obj.dimension );
      obj.edge_points(1:obj.nd_V,:) = obj.polygon.vertices; % WARN: Copy data
      
      pt_idx = reshape((obj.nd_V+1):(obj.nd_E+obj.nd_V),...
        [obj.degree-1,obj.polygon.nb_edges])';
      
      for k=1:obj.polygon.nb_edges
        
        obj.edge_points(pt_idx(k,:),1:obj.dimension) = Mappings.unity_to_segment( ...
          obj.polygon.vertices(obj.polygon.edges(k,1),:), ...
          obj.polygon.vertices(obj.polygon.edges(k,2),:), ...
          obj.quadrature_segment.points(3:end));
        
      end
     
      
      
      % Using alpha shapes: (for use with Pacman but issue tayloring Alpha)
      %al = alphaShape(obj.polygon.vertices); al.Alpha = 0.55;
       
      % Build quadrature rule for triangles
      obj.quadrature_triangle = ...
        Quadrature('Gauss-Legendre', 'triangle', 7 ); 
      % ATTENTION: Test for P2 only: otherwise change the number of Gauss
      % points
      
      % Build edge to dofs matrix.
      obj.edge_dofs = zeros( obj.polygon.nb_edges, obj.degree+1 );
      
      obj.edge_dofs(:,1:2) = obj.polygon.edges;
      obj.edge_dofs(:,3:(3+obj.degree-2)) = reshape((obj.nd_V+1):(obj.nd_E+obj.nd_V),...
        [obj.degree-1,obj.polygon.nb_edges])';
      
    end
    
    function plot_dofs(obj)
      % This function plot the polygon dofs:
      
      if obj.degree == 1
        %var{1} = 'filled';
%        obj.polygon.plot_polygon(obj.polygon.vertices, var)
         obj.polygon.plot_polyshape();
      else
        %h = patch(obj.polygon.vertices(:,1),obj.polygon.vertices(:,2), ...
        %  zeros(size(obj.polygon.vertices(:,2))));
        obj.polygon.plot_polyshape();
        axis equal; axis tight;
         hold on
         var1 = {'filled', 'r'};
         scatter(obj.polygon.vertices(:,1),obj.polygon.vertices(:,2),...
                 var1{:});
         %edges_dofs_toplot = reshape( obj.edge_points, ...
         %  [ size(obj.edge_points,1).*size(obj.edge_points,2),2] );
         var2 = {'d', 'filled' ,'b'};
         %scatter( edges_dofs_toplot(:,1), edges_dofs_toplot(:,2),var2{:});
         scatter( obj.edge_points((obj.nd_V+1):end,1),...
           obj.edge_points((obj.nd_V+1):end,2),var2{:});
         uistack(h, 'bottom')
      
      end
      
    end
    
    function d = compute_d_alpha_beta(obj,k)
      
      nkm2 = obj.dimension_polynomial_space(obj.degree-2); 
      % Dimension of M_{k-2}: n_{k-2}
      
      d = zeros(nkm2,1);
      
      h = obj.polygon.diameter;
      
      % Build dx^2 monomial
      dxx_monomial_index = ...
        [ obj.monomial_index(k,1)-2, obj.monomial_index(k,2) ];
      if (dxx_monomial_index(1) >= 0)
        d1 = obj.monomial_index(k,1)*(obj.monomial_index(k,1)-1)/h^2;
        % Get position in laplacian_monomial_index
        dxx_base_position = ...
          ismember(obj.laplacian_monomial_index,dxx_monomial_index,'rows');
        d(dxx_base_position) = d1;
      end
      % Build dy^2 monomial
      dyy_monomial_index = ...
        [ obj.monomial_index(k,1), obj.monomial_index(k,2)-2 ];
      if (dyy_monomial_index(2) >= 0)
        d2 = obj.monomial_index(k,2)*(obj.monomial_index(k,2)-1)/h^2;
        % Get position in laplacian_monomial_index
        dyy_base_position = ...
          ismember(obj.laplacian_monomial_index,dyy_monomial_index,'rows');
        d(dyy_base_position) = d(dyy_base_position) + d2;
      end
 
    end
      
    function nk = dimension_polynomial_space(obj,k)
      if obj.dimension == 2
        nk = (k+1)*(k+2)/2; % In 2D only
      else
        error('Dimension not support!')
      end
    end
    
    function bool = check_Pi_grad(obj)
      
      if isequal(ones(obj.nk,obj.nk), abs(obj.G-obj.B*obj.D) < sqrt(eps));
        bool = true;
      else
        bool = false;
      end
      
    end
    
    

  end
  
end
 


function d = delta(i,j)
if i == j
  d = 1;
  return
else
  d = 0;
end
end
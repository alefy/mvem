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

classdef Operators < handle
  %Operators The operators class builds left and right hand sides operators
  % and performs the assembly of their matricial representation.
  %   
  
  %------------------------------------------------------------------------
  % PROPERTIES
  %------------------------------------------------------------------------
  
  properties
    
    dofs % Holds the dofs class object 
    
    K % Stiffness matrix
    M % Mass matrix
    
    % Boolean checks if the operators have been built.
    built_K = false;
    built_M = false;
    
  end
  
  %------------------------------------------------------------------------
  % METHODS
  %------------------------------------------------------------------------
  
  methods
    
    % Constructor
    function obj = Operators(dofs)
      obj.dofs = dofs;
    end
    
    % Public methods
    function obj = build_K(obj)
      
      % Init
      obj.K = sparse(obj.dofs.nb_total_dofs,obj.dofs.nb_total_dofs);
      
      for k = 1:obj.dofs.mesh.nb_cells;
        
        % Build local stiffness matrix
        if isa(obj.dofs.mesh.cells{k},'Polygon') 
          % VEM
          
          Pi_nabla_star = obj.dofs.elements{k}.G\obj.dofs.elements{k}.B;
          Pi_nabla = obj.dofs.elements{k}.D*Pi_nabla_star;
          G_tilde = [ zeros(1, size(obj.dofs.elements{k}.G,2)) ;...
                      obj.dofs.elements{k}.G(2:end,:) ];
          
          alpha = 1;          
                    
          K_local = Pi_nabla_star'*G_tilde*Pi_nabla_star + ...
                    alpha*(eye(size(Pi_nabla))-Pi_nabla)'* ...
                    (eye(size(Pi_nabla))-Pi_nabla);
          
          % Assembly
          for i=1:obj.dofs.elements{k}.nd_total
            for j=1:obj.dofs.elements{k}.nd_total
              
              obj.K = obj.K + sparse( obj.dofs.get_dofs(k, i),...
                                      obj.dofs.get_dofs(k, j),...
                                      K_local(i,j), ...
                                      obj.dofs.nb_total_dofs,...
                                      obj.dofs.nb_total_dofs );
                                    
            end
          end
        end
        
      end
      obj.built_K = true;
    end
    
    function obj = build_M(obj)
      
      beta = 1;
      
      % Init
      obj.M = sparse(obj.dofs.nb_total_dofs,obj.dofs.nb_total_dofs);

      for k = 1:obj.dofs.mesh.nb_cells;
        if obj.dofs.elements{k}.degree  < 3
          % In this case: Pi^Nabla = Pi^0
          Pi_nabla_star = obj.dofs.elements{k}.G\obj.dofs.elements{k}.B;
          Pi_nabla = obj.dofs.elements{k}.D*Pi_nabla_star;
          
          M_local = Pi_nabla_star'*obj.dofs.elements{k}.H*Pi_nabla_star;
        else
          error('The L^2 operator is not yet supported for k>2.')
        end
        % Assembly
        for i=1:obj.dofs.elements{k}.nd_total
          for j=1:obj.dofs.elements{k}.nd_total
            
            obj.M = obj.M + sparse( obj.dofs.get_dofs(k, i),...
              obj.dofs.get_dofs(k, j),...
              M_local(i,j), ...
              obj.dofs.nb_total_dofs,...
              obj.dofs.nb_total_dofs );
            
          end
        end
                
      end

      obj.built_M = true;
      
    end
    

    
  end
  
end


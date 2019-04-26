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

classdef Mappings
  %
  % The class Mappings defines projections from reference geometries to
  % current ones. The main function is to map Gauss points.
  %
  % At first only affine mappings are supported. 
 
  %------------------------------------------------------------------------
  % METHODS
  %------------------------------------------------------------------------
  
  methods(Static)
    
    function Y = unity_to_segment(x1,x2,xi)
      % This function projects the point xi in [0,1] 
      % to current segment a point X
      %
      % It supports xi as a set of points
      
      Y = zeros([length(xi),2]);
      for i=1:length(xi) % How can I get ride of that loop? 
                         % (should hold for in 3D two!)
        Y(i,:) = x1*(1-xi(i)) + x2*xi(i);
      end
    end
    
  end
  
  
end
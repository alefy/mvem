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

classdef Load < handle
  % The load class will deal with the construction of the load (lhs) term
  %   
  
  %------------------------------------------------------------------------
  % PROPERTIES
  %------------------------------------------------------------------------
  
  properties
    b % Right Hand Side vector
  end
  
  %------------------------------------------------------------------------
  % METHODS
  %------------------------------------------------------------------------
  
  methods
    % Only a homogeneous constant load is supported.
    function obj = build_constant_load(obj, operators, coef)
      % Requires the operator class for the L^2 projection
      one = coef.*ones(length(operators.M),1);
      obj.b = operators.M*one;
    end
    
  end
  
  
end


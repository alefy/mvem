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


% Create a polygonal mesh of a circle: to be used with PolyMesher 
% (version 1.1)
%

function x = PolyCircle(Demand, Arg)
BdBox = [-1 1 -1 1];

switch(Demand)
  case('Dist');  x = DistFnc(Arg, BdBox);
  case('BC');    x = BndryCnds(Arg{:}, BdBox);
  case('BdBox'); x = BdBox;
  case('PFix');  x = FixedPoints(BdBox);
end
end

function Dist = DistFnc(P, BdBox)
Dist = dCircle(P,0,0,1);
end

function x = BndryCnds(Node, Element, BdBox)

  eps = 0.5*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  CircleNodes = find(abs(sqrt(Node(:,1).^2+Node(:,2).^2)-1)<eps);
  
  Supp = ones(size(CircleNodes,1),3);
  Supp(:,1) = CircleNodes;

  Load = [1, 0, 0];
  x = {Supp,Load};
  
end

function [PFix] = FixedPoints(BdBox)
  PFix = [];
end
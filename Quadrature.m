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

classdef Quadrature
  %
  % The class Quadrature manages numerical integration.
  %
  %   Supported geometry:
  %     > Segments (0,1)
  %     Supported integration methods:
  %       - Gauss-Legendre (Arbitrary precision)
  %       - Gauss-Lobatto (up to precision 5)
  %     > Triangles
  %       - up to precision 6
  %     > Quadrilaterals (to be implemented)
  %
  
  properties
    geometry % segment (0,1) triangle quadrilateral
    
    type; % First implemented: Gauss-Legendre
    
    nb_points
    integration_order % integration order
    
    points
    weights 
  end
  
  methods
    
    % Constructor
    function obj = Quadrature( type, geometry, nb_points )
        
      obj.geometry = geometry;
      obj.type = type;
      obj.nb_points = nb_points;
      
      
      if strcmp(obj.geometry,'segment')
          switch obj.type
            case 'Gauss-Legendre'
              obj.integration_order = 2*obj.nb_points - 1 ;
              if obj.nb_points < 5
                [obj.points, obj.weights] = ...
                  Quadrature.Gauss_Legendre_exact( 0,1,obj.nb_points );
              else
                [obj.points, obj.weights] = ...
                  Quadrature.Gauss_Legendre_generate( ...
                  0, 1, obj.nb_points, sqrt(eps) );
              end
            case 'Gauss-Lobatto'
              if obj.nb_points < 2
                error('Gauss-Lobatto formulae not enought points: (<3-1)')
              elseif obj.nb_points > 6
                error( [ 'Gauss-Lobatto formulae not supported ' ...
                         'for this number of points.'])
              end
              obj.integration_order = 2*obj.nb_points - 3 ;
              [obj.points, obj.weights] = ...
                Quadrature.Gauss_Lobatto_exact( 0,1,obj.nb_points );
            otherwise
              error('Integration formulae not supported.')
          end
      elseif strcmp(obj.geometry, 'triangle')
        switch obj.type
          case 'Gauss-Legendre'
            [obj.points, obj.weights, obj.integration_order] = ...
              Quadrature.Gauss_triangle_barycentric(nb_points);
          otherwise
            error('Numberical integration for triangle not supported.')
        end
      else
        error('Numerical integration on such a geometry is not supported.')
      end
    end
    
  end
  
  methods(Static)
    % I leave access from outside to quadrature rules.
    function [p, w] = Gauss_Legendre_exact( a, b, nb_points)
      % Gaussian rules (cf Ern & Guermond 2004)
      m = 0.5*(a+b);
      d = b-a;
      
      switch nb_points

        case 1
          p = m;
          w = d;
        case 2
          p(1) = m-d/2*sqrt(3)/3;
          p(2) = m+d/2*sqrt(3)/3;
          w(1) = 0.5*d;
          w(2) = 0.5*d;
        case 3
          p(1) = m - d/2*sqrt(3/5);
          p(2) = m;
          p(3) = m + d/2*sqrt(3/5);
          w(1) = 5/18*d;
          w(2) = 8/18*d;
          w(3) = 5/18*d;
        case 4
          p(1) = m - d/2*sqrt( ( 15 + 2*sqrt(30) ) / 35 );
          p(2) = m + d/2*sqrt( ( 15 + 2*sqrt(30) ) / 35 );
          p(3) = m - d/2*sqrt( ( 15 - 2*sqrt(30) ) / 35 );
          p(4) = m + d/2*sqrt( ( 15 - 2*sqrt(30) ) / 35 );
          w(1) = ( 1/4 - sqrt(5/6)/12 )*d;
          w(2) = ( 1/4 - sqrt(5/6)/12 )*d;
          w(3) = ( 1/4 + sqrt(5/6)/12 )*d;
          w(4) = ( 1/4 + sqrt(5/6)/12 )*d;
        otherwise
          error('Number of Gauss points no supported for exact rule')
      end
    end
    
    function [p, w] = Gauss_Legendre_generate( a, b, n_pts, tol)
      % Generate Gauss points using the Newton-Method
      % Check for right number of arguments

      p = zeros(n_pts,1);
      w = zeros(n_pts,1);
      
      if(nargin==3)
        tol = eps;
      end
  
      % Initalize variables
  
      m = floor((n_pts+1)/2);
      xm = (b+a)/2;
      xl = (b-a)/2;
   
      for i = 1:m
    
        % Initial guess of root (starting value)
        
        z = cos(pi*(i-1/4)/(n_pts+1/2));
    
        delta = tol+1;
        while(delta > tol)
          
          p1 = 0;
          p2 = 1;
          
          for k = 0:(n_pts-1)
        
            % Computing value of n-th Legendre polynomial at point z using the
            % recursion:
            %
            %   (j+1)*P_(j+1)(z) = (2*j+1)*z*P_(j)(z)-j*P_(j-1)(z)
          
            p3 = ((2*k+1)*z*p2-k*p1)/(k+1);
        
            % Computing value of first derivative of n-th Legendre polynomial
            % at point z using the recursion:
            %
            %   (1-z^2)*P'_(j)(z) = j*[z*P_(j)(z)-P_(j-1)(z)]
            
            dp = n_pts*(z*p3-p2)/(z^2-1);
            p1 = p2;
            p2 = p3;
            
          end    
      
          % Performing Newton update
      
          z_old = z;
          z = z_old-p3/dp;
      
          delta = abs(z-z_old);
      
        end
    
        % Computing wieghts and abscissae
    
        p(i) = xm-xl*z;
        p(n_pts+1-i) = xm+xl*z;
        w(i) = 2*xl/((1-z^2)*dp^2);
        w(n_pts+1-i) = w(i);
        
      end
    
    end
    
    function [p, w] = Gauss_Lobatto_exact( a, b, nb_points)
      m = 0.5*(a+b);
      d = b-a;
      switch nb_points
        case 2
          % Trapezium rule
          p(1) = a;
          p(2) = b;
          w(1) = 1/2;
          w(2) = 1/2;
        case 3
          p = zeros(1,nb_points);
          w = zeros(1,nb_points);
          % First compute Gauss points with one Gauss points with exact
          % rule
          p(1) = a;
          p(2) = b;
          p(3) = m;
          
          w(1) = d/6;
          w(2) = d/6;
          w(3) = 2*d/3;
        case 4
          p(1) = a;
          p(2) = b;
          p(3) = m - d/10*sqrt(5);
          p(4) = m + d/10*sqrt(5);
          
          w(1) = d*1/12;
          w(2) = d*1/12;
          w(3) = d*5/12;
          w(4) = d*5/12;
        case 5
          p(1) = a; 
          p(2) = b;
          p(3) = m - d/14*sqrt(21);
          p(4) = m;
          p(5) = m + d/14*sqrt(21);
          
          w(1) = d/20;
          w(2) = d/20;
          w(3) = d*49/180;
          w(4) = d*32/90;
          w(5) = d*49/180;
          
        case 6
          p(1) = a;
          p(2) = b;
          p(3) = m - d/2*sqrt( (7+2*sqrt(7))/21 );
          p(4) = m - d/2*sqrt( (7-2*sqrt(7))/21 );
          p(5) = m + d/2*sqrt( (7-2*sqrt(7))/21 );
          p(6) = m + d/2*sqrt( (7+2*sqrt(7))/21 );
          
          w(1) = d/30;
          w(2) = d/30;
          w(3) = d*(14-sqrt(7))/60;
          w(4) = d*(14+sqrt(7))/60;
          w(5) = d*(14+sqrt(7))/60;
          w(6) = d*(14-sqrt(7))/60;
          
        otherwise
          error('Gauss-Lobatto integration error.')
      end
    end
    
    function [p, w, precision] = Gauss_triangle_barycentric(nb_points)
      % p: matrix (nb_points times 3) (p=[r ; s ; t]')
      % Such a choice of the organisation of the matrix p result from the
      % construction of the class delaunayTriangulation cf. 
      % r, s and t: barycentric Gauss points coordinates
      % w: vector (length nb_points)
      switch nb_points
        case 1 % degree of precision 0
          precision = 0;
          
          r = 1/3;
          s = 1/3;
          t = 1/3;
          
          w = 1.0;
          
          p = [ r ; s ; t]';
        case 2
          error('No integration rules for triangles for 2 Gauss points')
        case 3 % degree of precision 2
          precision = 2;
          % r = [0  1/2  1/2];
          % s = [1/2 0   1/2];
          % t = [ 1/2  1/2 0];
          
          r = [2/3 1/6 1/6];
          s = [1/6 2/3 1/6];
          t = [1/6 1/6 2/3];
          %
          % Vertex rule precision 1
          % r = [1 0 0];
          % s = [0 1 0];
          % t = [0 0 1];
        
          w = [1/3 1/3 1/3];

          p = [ r ; s ; t]';
          
        case 4 % degree of precision 3
          precision = 3;
          
          r = [1/3 0.6 0.2 0.2];
          s = [1/3 0.2 0.6 0.2];
          t = [1/3 0.2 0.2 0.6];
          
          w = [-27/48 25/48 25/48 25/48];
          
          p = [ r ; s ; t]';
        case 5
          error('No integration rules for triangles for 5 Gauss points')
        case 6 % degree of precision: 4
          
          precision = 4;
          
          r = [0.816847572980459 0.091576213509771 0.091576213509771 ...
               0.108103018168070 0.445948490915965 0.445948490915965];
          s = [0.091576213509771 0.816847572980459 0.091576213509771 ...
               0.4459484909115965 0.108103018170459 0.445948490915965];
          t = [0.091576213509771 0.091576213509771 0.816847572980459 ...
               0.4459484909115965 0.445948490915965 0.108103018168070];

          w = [0.109951743655322 0.109951743655322 0.109951743655322 ...
               0.223381589678011 0.223381589678011 0.223381589678011];
          
          p = [ r ; s ; t]';
        case 7 % Precision Degree: 5
          
          precision = 5;
          
          r = [1/3 ...
               0.059715871789770 0.470142064105115 0.470142064105115 ...
               0.797426985353087 0.101286507323456 0.101286507323456 ];
          s = [1/3 ...
               0.470142064105115 0.059715871789770 0.470142064105115 ...
               0.101286507323456 0.797426985353087 0.101286507323456 ];
          t = [1/3 ...
               0.470142064105115 0.470142064105115 0.059715871789770 ...
               0.101286507323456 0.101286507323456 0.797426985353087 ];
          
          w = [0.225000000000000 ...
               0.132394152788506 0.132394152788506 0.132394152788506 ...
               0.125939180544827 0.125939180544827 0.125939180544827 ];
          
          p = [ r ; s ; t]';
       
        case 9
          
          precision = 6;
          
          r = [ 0.124949503233232 0.437525248383384 0.437525248383384...
                0.797112651860071 0.037477420750088 0.165409927389841...
                0.797112651860071 0.037477420750088 0.165409927389841 ];
          
          s = [ 0.437525248383384 0.124949503233232 0.437525248383384...
                0.165409927389841 0.797112651860071 0.037477420750088...
                0.165409927389841 0.797112651860071 0.037477420750088 ];
          
          t = [ 0.437525248383384 0.437525248383384 0.124949503233232...
                0.037477420750088 0.165409927389841 0.797112651860071...
                0.037477420750088 0.165409927389841 0.797112651860071 ];
          
          w = [ 0.205950504760887 0.205950504760887 0.205950504760887 ...
                0.063691414286223 0.063691414286223 0.063691414286223 ...
                0.063691414286223 0.063691414286223 0.063691414286223 ];
          
          p = [ r ; s ; t]';
          
        case 12 % Precision Degree: 6
          
          precision = 6;
          
          r = [ 0.873821971016996 0.063089014491502 0.063089014491502 ...
                0.501426509658179 0.249286745170911 0.249286745170910 ...
                0.636502499121399 0.053145049844816 0.310352451033785 ...
                0.636502499121399 0.053145049844816 0.310352451033785 ];
            
          s = [ 0.063089014491502 0.873821971016996 0.063089014491502 ...
                0.249286745170910 0.501426509658179 0.249286745170911 ...
                0.310352451033785 0.636502499121399 0.053145049844816 ...
                0.310352451033785 0.636502499121399 0.053145049844816 ];
          
          t = [ 0.063089014491502 0.063089014491502 0.873821971016996 ...
                0.249286745170911 0.249286745170910 0.501426509658179 ...
                0.053145049844816 0.310352451033785 0.636502499121399 ...
                0.053145049844816 0.310352451033785 0.636502499121399 ];
              
          w = [ 0.050844906370207 0.05084490637020 0.050844906370207 ...
                0.116786275726379 0.116786275726379 0.116786275726379 ...
                0.082851075618374 0.082851075618374 0.082851075618374 ...
                0.082851075618374 0.082851075618374 0.082851075618374 ];
          
          p = [ r ; s ; t]';
          
        case 13
          precision = 7;
          
          r = [ 1/3 ...
                0.479308067841923 0.260345966079038 0.260345966079038 ...
                0.869739794195568 0.065130102902216 0.065130102902216 ...
                0.638444188569809 0.086903154253160 0.312865496004875 ...
                0.638444188569809 0.086903154253160 0.312865496004875 ];
          
          s = [ 1/3 ...
                0.260345966079038 0.479308067841923 0.260345966079038 ...
                0.065130102902216 0.869739794195568 0.065130102902216 ...
                0.312865496004875 0.638444188569809 0.086903154253160 ...
                0.312865496004875 0.638444188569809 0.086903154253160 ];
          
          t = [ 1/3 ...
                0.260345966079038 0.260345966079038 0.479308067841923 ...
                0.065130102902216 0.065130102902216 0.869739794195568 ...
                0.086903154253160 0.312865496004875 0.638444188569809 ...
                0.086903154253160 0.312865496004875 0.638444188569809 ];
          
          w = [ -0.149570044467670 ...
                0.175615257433204 0.175615257433204 0.175615257433204 ...
                0.053347235608839 0.053347235608839 0.053347235608839 ...
                0.077113760890257 0.077113760890257 0.077113760890257 ...
                0.077113760890257 0.077113760890257 0.077113760890257 ];
               
          p = [ r ; s ; t]';
        otherwise
          error('Error in the choice of number of Gauss points')
      end
      
      
      
    end
    
  end

end
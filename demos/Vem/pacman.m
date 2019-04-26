function [vertices, edges] = pacman(nb_points)

close all
%nb_points = 30;

radius = 0.5;

aperture = pi/6;
alpha = (2*pi-2*aperture)/(nb_points-3);

vertices = zeros(nb_points-1,2);
edges = zeros( nb_points, 2);

vertices(1,:) = [ 0 , 0];
vertices(2,:) = [ radius*cos(aperture) , radius*sin(aperture)];
for k = 3:(nb_points-1)
  vertices(k,1) = radius*cos((k-2)*alpha+aperture);
  vertices(k,2) = radius*sin((k-2)*alpha+aperture);
end

scatter(vertices(:,1),vertices(:,2));
axis([-1,1,-1,1])
axis equal;

% Build vertices
edges = [ 1:(nb_points-2) ; 2:(nb_points-1)]';
edges(end+1,:) = [ nb_points-1, 1];

end
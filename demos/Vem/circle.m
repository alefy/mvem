function [vertices, edges] = circle(radius, nb_points)

%close all
%nb_points = 30;

aperture = pi/6;
alpha = (2*pi)/(nb_points);

vertices = zeros(nb_points,2);
edges = zeros( nb_points+1, 2);

vertices(1,:) = [ 0 , 0];
vertices(2,:) = [ cos(aperture) , sin(aperture)];
for k = 1:(nb_points)
  vertices(k,1) = radius*cos(k*alpha);
  vertices(k,2) = radius*sin(k*alpha);
end

scatter(vertices(:,1),vertices(:,2));
axis([-1,1,-1,1])
axis equal;

% Build vertices
edges = [ 1:(nb_points-2) ; 2:(nb_points-1)]';
edges(end+1,:) = [ nb_points-1, 1];

end
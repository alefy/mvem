% Pentagon

% Cf. The hitchhiker's guide to the virtual element method

%clear
close all

vert_penta = [0,0;3,0;3,2;3/2,4;0,4];
edge_penta = [ 1 2 ; 2 3; 3 4; 4 5 ; 5 1 ];

penta_poly = Polygon(vert_penta, edge_penta);
vem_pentagon = Vem(penta_poly,3); % Quadratic VEM

vem_pentagon.plot_dofs()
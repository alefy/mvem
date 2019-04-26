% Unit square 
%
% Cf. The hitchhiker's guide to the virtual element method

clear
close all

vert_square = [0,0 ; 1,0 ; 1,1 ; 0,1];
edge_square = [ 1 2 ; 2 3 ; 3 4 ; 4 1 ];

square = Polygon(vert_square, edge_square);
square.plot_polyshape();

vem_square = Vem(square,2); % Quadratic VEM

vem_square.plot_dofs
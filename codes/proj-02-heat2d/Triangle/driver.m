clear all; clc;

kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

% quadrature rule for triangles
% Using 1-point quadrature for linear triangles
xi = [1/3, 1/3];
weight = 1/2;

% mesh generation
n_en   = 3;               % number of nodes in an element (triangle)
n_el_x = 60;               % number of elements in x-dir (originally quadrilaterals)
n_el_y = 60;               % number of elements in y-dir
n_el   = n_el_x * n_el_y * 2; % total number of triangular elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array for triangular elements
IEN = zeros(n_el, n_en);
count = 1;
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    n1 = (ey-1)*n_np_x + ex;
    n2 = n1 + 1;
    n3 = n1 + n_np_x;
    n4 = n3 + 1;
    
    % First triangle
    IEN(count, :) = [n1, n2, n4];
    count = count + 1;
    
    % Second triangle
    IEN(count, :) = [n1, n4, n3];
    count = count + 1;
  end
end

% ID array
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter;

LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over elements to assemble the matrix and vector
for ee = 1 : n_el
  nodes = IEN(ee, 1:n_en);
  x_ele = x_coor(nodes);
  y_ele = y_coor(nodes);
  
  k_ele = zeros(n_en, n_en); % element stiffness matrix
  f_ele = zeros(n_en, 1);    % element load vector
  
  % Compute area of the triangle
  detJ = (x_ele(2)-x_ele(1))*(y_ele(3)-y_ele(1)) - (x_ele(3)-x_ele(1))*(y_ele(2)-y_ele(1));
  area = abs(detJ) / 2;
  
  % Loop over quadrature points
  % For linear triangles with 1-point quadrature
  L1 = xi(1);
  L2 = xi(2);
  L3 = 1 - L1 - L2;
  
  % Compute shape functions at quadrature point
  N = [Triangle(1, L1, L2, L3), Triangle(2, L1, L2, L3), Triangle(3, L1, L2, L3)];
  
  % Compute physical coordinates of quadrature point
  x_l = N(1)*x_ele(1) + N(2)*x_ele(2) + N(3)*x_ele(3);
  y_l = N(1)*y_ele(1) + N(2)*y_ele(2) + N(3)*y_ele(3);
  
  % Compute gradients of shape functions
  % For linear triangles, gradients are constant
  % Using the formula: grad N = [dN/dx, dN/dy]
  % dN/dx = (y2 - y3) / (2*area)
  % dN/dy = (x3 - x2) / (2*area)
  
  dN_dx = [(y_ele(2) - y_ele(3)), (y_ele(3) - y_ele(1)), (y_ele(1) - y_ele(2))] / detJ;
  dN_dy = [(x_ele(3) - x_ele(2)), (x_ele(1) - x_ele(3)), (x_ele(2) - x_ele(1))] / detJ;
  
  % Compute f_ele
  f_val = f(x_l, y_l);
  f_ele = f_val * N' * weight * area * 2; % weight=1/2 for area
  
  % Compute k_ele
  for aa = 1 : n_en
    for bb = 1 : n_en
      k_ele(aa, bb) = k_ele(aa, bb) + kappa * (dN_dx(aa)*dN_dx(bb) + dN_dy(aa)*dN_dy(bb)) * area * 2; % area * 2 = |detJ|
    end
  end
  
  % Assemble into global K and F
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

% save the solution vector and number of elements to disp with name
% HEAT.mat
save("HEAT", "disp", "n_el_x", "n_el_y");

% EOF

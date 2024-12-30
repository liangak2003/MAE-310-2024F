clear all; clc;

% Define exact solution and its derivatives
exact = @(x,y) sin(pi*x).*sin(pi*y); % Manufactured solution
exact_x = @(x,y) pi*cos(pi*x).*sin(pi*y);
exact_y = @(x,y) pi*sin(pi*x).*cos(pi*y);

kappa = 1.0; % conductivity

% Source term based on the manufactured solution
f = @(x,y) 2.0*kappa*pi^2*sin(pi*x).*sin(pi*y); % source term derived from exact solution

% Quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% Define different mesh sizes for convergence study
mesh_sizes = [10, 20, 40, 80]; % Number of elements in each direction
num_mesh = length(mesh_sizes);

% Preallocate error arrays
L2_errors = zeros(num_mesh,1);
H1_errors = zeros(num_mesh,1);
h_values = zeros(num_mesh,1);

for m = 1:num_mesh
    n_el_x = mesh_sizes(m);
    n_el_y = mesh_sizes(m);
    h = 1.0 / n_el_x;
    h_values(m) = h;
    
    n_el   = n_el_x * n_el_y; % total number of elements
    
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
    
    % IEN array
    IEN = zeros(n_el, 4);
    for ex = 1 : n_el_x
      for ey = 1 : n_el_y
        ee = (ey-1) * n_el_x + ex; % element index
        IEN(ee, 1) = (ey-1) * n_np_x + ex;
        IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
        IEN(ee, 3) =  ey    * n_np_x + ex + 1;
        IEN(ee, 4) =  ey    * n_np_x + ex;
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
      x_ele = x_coor( IEN(ee, 1:4) );
      y_ele = y_coor( IEN(ee, 1:4) );
      
      k_ele = zeros(4, 4); % element stiffness matrix
      f_ele = zeros(4, 1);    % element load vector
      
      for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : 4
          x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
          y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
          [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
          dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
          dx_deta = dx_deta + x_ele(aa) * Na_eta;
          dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
          dy_deta = dy_deta  + y_ele(aa) * Na_eta;
        end
        
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
        
        for aa = 1 : 4
          Na = Quad(aa, xi(ll), eta(ll));
          [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
          Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
          Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
          
          f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
          
          for bb = 1 : 4
            Nb = Quad(bb, xi(ll), eta(ll));
            [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
            Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
            Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
            
            k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
          end % end of bb loop
        end % end of aa loop
      end % end of quadrature loop
     
      for aa = 1 : 4
        PP = LM(ee, aa);
        if PP > 0
          F(PP) = F(PP) + f_ele(aa);
          
          for bb = 1 : 4
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
    disp_num = zeros(n_np, 1);
    
    for ii = 1 : n_np
      index = ID(ii);
      if index > 0
        disp_num(ii) = dn(index);
      else
        % modify disp with the g data. Here it does nothing because g is zero
      end
    end
    
    % Calculate errors
    [L2_error, H1_error] = error_calc(x_coor, y_coor, disp_num, exact, exact_x, exact_y, n_el_x, n_el_y);
    
    L2_errors(m) = L2_error;
    H1_errors(m) = H1_error;
    
    fprintf('Mesh size: %d x %d, h = %.5f, L2 error = %.5e, H1 error = %.5e\n', ...
            n_el_x, n_el_y, h, L2_error, H1_error);
    
end

% Plot convergence
figure;
loglog(h_values, L2_errors, '-o', 'LineWidth',2);
hold on;
loglog(h_values, H1_errors, '-s', 'LineWidth',2);
xlabel('Mesh size h');
ylabel('Error');
legend('L2 Error', 'H1 Error');
title('Convergence of FEM Solution');
grid on;
hold off;

% Save the solution for the finest mesh (optional)
% save("HEAT_convergence.mat", "disp_num", "n_el_x", "n_el_y");

% EOF

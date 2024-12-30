clear all; clc;

% Define Exact Solution and Its Derivatives
exact = @(x,y) x.*(1-x).*y.*(1-y); % Manufactured Exact Solution
exact_x = @(x,y) (1-2*x).*y.*(1-y);
exact_y = @(x,y) x.*(1-x).*(1-2*y);

kappa = 1.0; % Thermal Conductivity

% Source Term Based on Manufactured Solution
f = @(x,y) 2.0*kappa*(x.*(1 - x) + y.*(1 - y)); % Derived from Exact Solution

% Define Different Mesh Sizes for Convergence Study
mesh_sizes = [10, 20, 40, 80]; % Number of Elements in Each Direction
num_mesh = length(mesh_sizes);

% Preallocate Error Arrays
L2_errors = zeros(num_mesh,1);
H1_errors = zeros(num_mesh,1);
h_values = zeros(num_mesh,1);

for m = 1:num_mesh
    n_el_x = mesh_sizes(m);
    n_el_y = mesh_sizes(m);
    h = 1.0 / n_el_x;
    h_values(m) = h;
    
    n_el = n_el_x * n_el_y * 2; % Total Number of Triangle Elements
    
    n_np_x = n_el_x + 1;      % Number of Nodes in x-direction
    n_np_y = n_el_y + 1;      % Number of Nodes in y-direction
    n_np   = n_np_x * n_np_y; % Total Number of Nodes
    
    x_coor = zeros(n_np, 1);
    y_coor = zeros(n_np, 1);
    
    hx = 1.0 / n_el_x;        % Mesh Size in x-direction
    hy = 1.0 / n_el_y;        % Mesh Size in y-direction
    
    % Generate Node Coordinates
    for ny = 1 : n_np_y
        for nx = 1 : n_np_x
            index = (ny-1)*n_np_x + nx; % Node Index
            x_coor(index) = (nx-1) * hx;
            y_coor(index) = (ny-1) * hy;
        end
    end
    
    % IEN Array (Element-to-Node Mapping) Using Triangle Elements
    IEN = zeros(n_el, 3);
    count = 1;
    for ey = 1 : n_el_y
        for ex = 1 : n_el_x
            n1 = (ey-1)*n_np_x + ex;
            n2 = n1 + 1;
            n3 = n1 + n_np_x;
            n4 = n3 + 1;
            
            % First Triangle
            IEN(count, :) = [n1, n2, n4];
            count = count + 1;
            
            % Second Triangle
            IEN(count, :) = [n1, n4, n3];
            count = count + 1;
        end
    end
    
    % ID Array (Identification of Interior Nodes)
    ID = zeros(n_np,1);
    counter = 0;
    for ny = 2 : n_np_y - 1
        for nx = 2 : n_np_x - 1
            index = (ny-1)*n_np_x + nx;
            counter = counter + 1;
            ID(index) = counter;  
        end
    end
    
    n_eq = counter; % Number of Equations
    
    LM = ID(IEN); % Linking Matrix
    
    % Allocate Stiffness Matrix and Load Vector
    K = spalloc(n_eq, n_eq, 9 * n_eq);
    F = zeros(n_eq, 1);
    
    % Assemble Stiffness Matrix and Load Vector
    for ee = 1 : n_el
        nodes = IEN(ee, :);
        x_ele = x_coor(nodes);
        y_ele = y_coor(nodes);
        
        k_ele = zeros(3, 3); % Element Stiffness Matrix
        f_ele = zeros(3, 1); % Element Load Vector
        
        % Compute Jacobian Determinant and Element Area
        detJ = (x_ele(2)-x_ele(1))*(y_ele(3)-y_ele(1)) - (x_ele(3)-x_ele(1))*(y_ele(2)-y_ele(1));
        area = abs(detJ) / 2;
        
        % Compute Shape Function Gradients
        dN_dx = [(y_ele(2) - y_ele(3)), (y_ele(3) - y_ele(1)), (y_ele(1) - y_ele(2))] / detJ;
        dN_dy = [(x_ele(3) - x_ele(2)), (x_ele(1) - x_ele(3)), (x_ele(2) - x_ele(1))] / detJ;
        
        % Use Single Point Gaussian Integration (Centroid)
        L1 = 1/3;
        L2 = 1/3;
        L3 = 1 - L1 - L2;
        
        % Shape Functions at Centroid
        N = [Triangle(1, L1, L2, L3), Triangle(2, L1, L2, L3), Triangle(3, L1, L2, L3)];
        
        % Physical Coordinates at Centroid
        x_l = N(1)*x_ele(1) + N(2)*x_ele(2) + N(3)*x_ele(3);
        y_l = N(1)*y_ele(1) + N(2)*y_ele(2) + N(3)*y_ele(3);
        
        % Compute Load Vector
        f_val = f(x_l, y_l);
        f_ele = f_val * N' * area; % Correctly Apply Area Factor
        
        % Compute Stiffness Matrix
        for aa = 1 : 3
            for bb = 1 : 3
                k_ele(aa, bb) = k_ele(aa, bb) + kappa * (dN_dx(aa)*dN_dx(bb) + dN_dy(aa)*dN_dy(bb)) * area; % Correctly Apply Area Factor
            end
        end
        
        % Assemble into Global K and F
        for aa = 1 : 3
            PP = LM(ee, aa);
            if PP > 0
                F(PP) = F(PP) + f_ele(aa);
                
                for bb = 1 : 3
                    QQ = LM(ee, bb);
                    if QQ > 0
                        K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
                    else
                        % Boundary conditions are zero or homogeneous; no action needed
                    end
                end  
            end
        end
    end
    
    % Solve Stiffness Matrix
    dn = K \ F;
    
    % Assign Solution Back to All Nodes
    disp_num = zeros(n_np, 1);
    
    for ii = 1 : n_np
        index = ID(ii);
        if index > 0
            disp_num(ii) = dn(index);
        else
            % Boundary conditions are zero or homogeneous; no action needed
        end
    end
    
    % Compute Errors
    [L2_error, H1_error] = error_calc(x_coor, y_coor, disp_num, exact, exact_x, exact_y, n_el_x, n_el_y, IEN);
    
    L2_errors(m) = L2_error;
    H1_errors(m) = H1_error;
    
    fprintf('Mesh size: %d x %d, h = %.5f, L2 error = %.5e, H1 error = %.5e\n', ...
            n_el_x, n_el_y, h, L2_error, H1_error);
end

% Plot Convergence Graph
figure;
loglog(h_values, L2_errors, '-o', 'LineWidth',2);
hold on;
loglog(h_values, H1_errors, '-s', 'LineWidth',2);

% Add Theoretical Convergence Rate Lines
ref_h = h_values;
loglog(ref_h, (h_values(1)/h_values).^2 * L2_errors(1), '--k', 'LineWidth',1); % O(h^2)
loglog(ref_h, (h_values(1)/h_values) * H1_errors(1), '--k', 'LineWidth',1);   % O(h)

xlabel('Mesh Size h');
ylabel('Error');
legend('L2 Error', 'H1 Error');
title('Convergence of Finite Element Solution Using Triangular Elements');
grid on;
hold off;

% EOF

% e.m
clear all; clc; clf; % Clear memory, command window, and figure window

% Problem definition
f = @(x) -20*x.^3; % Source function f(x)
g = 1.0;           % Dirichlet boundary condition u(1) = g
h = 0.0;           % Neumann boundary condition -u_x(0) = h

% Using cubic elements
pp = 3; % Cubic elements
n_el = 4; % Number of elements (adjust as needed)
n_en = pp + 1; % Number of nodes per element
n_np = n_el * pp + 1; % Total number of nodes
n_eq = n_np - 1; % Degrees of freedom (excluding Dirichlet node)

hh = 1.0 / (n_np - 1); % Mesh size
x_coor = linspace(0, 1, n_np); % Node coordinates

% Construct the IEN array
IEN = zeros(n_el, n_en);
for ee = 1 : n_el
    for aa = 1 : n_en
        IEN(ee, aa) = (ee - 1) * pp + aa;
    end
end

% Setup ID array
ID = 1 : n_np;
ID(end) = 0; % Dirichlet boundary condition applied at the last node

% Test different numbers of quadrature points
quadrature_points = [1, 2, 3, 4, 5, 6];
errors_L2 = zeros(length(quadrature_points), 1);

for idx = 1:length(quadrature_points)
    n_int = quadrature_points(idx);
    [xi, weight] = Gauss(n_int, -1, 1);

    % Allocate stiffness matrix and load vector
    K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
    F = zeros(n_eq, 1);

    % Assembly process
    for ee = 1 : n_el
        k_ele = zeros(n_en, n_en); % Element stiffness matrix
        f_ele = zeros(n_en, 1);    % Element load vector

        x_ele = x_coor(IEN(ee,:)); % Element node coordinates

        % Quadrature loop
        for qua = 1 : n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            for aa = 1 : n_en
                N = PolyShape(pp, aa, xi(qua), 0);
                dN_dxi = PolyShape(pp, aa, xi(qua), 1);
                x_l    = x_l    + x_ele(aa) * N;
                dx_dxi = dx_dxi + x_ele(aa) * dN_dxi;
            end
            dxi_dx = 1.0 / dx_dxi;

            for aa = 1 : n_en
                N = PolyShape(pp, aa, xi(qua), 0);
                dN_dxi = PolyShape(pp, aa, xi(qua), 1);
                f_ele(aa) = f_ele(aa) + weight(qua) * N * f(x_l) * dx_dxi;
                for bb = 1 : n_en
                    dN_dxi_bb = PolyShape(pp, bb, xi(qua), 1);
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * dN_dxi * dN_dxi_bb * dxi_dx;
                end
            end
        end

        % Assemble into global matrix
        for aa = 1 : n_en
            P = ID(IEN(ee,aa));
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : n_en
                    Q = ID(IEN(ee,bb));
                    if(Q > 0)
                        K(P, Q) = K(P, Q) + k_ele(aa, bb);
                    else
                        F(P) = F(P) - k_ele(aa, bb) * g; % Handle Dirichlet BC
                    end
                end
            end
        end
    end

    % Apply Neumann boundary condition at x = 0
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

    % Solve the linear system
    d_temp = K \ F;
    disp = [d_temp; g]; % Append Dirichlet BC

    % Compute errors
    [e_L2, ~] = compute_errors(n_el, n_en, pp, IEN, disp, x_coor);
    errors_L2(idx) = e_L2;
    fprintf('Quadrature points: %d, L2 Error: %e\n', n_int, e_L2);
end

% Plot error graph
figure;
plot(quadrature_points, errors_L2, '-o', 'LineWidth', 2);
xlabel('Number of Quadrature Points');
ylabel('Relative L^2 Error');
title('Effect of Quadrature Points on Error (Cubic Elements)');
grid on;

% Embed the compute_errors function
function [e_L2, e_H1] = compute_errors(n_el, n_en, pp, IEN, disp, x_coor)
    syms x_sym;
    u_exact = x_sym^5;
    u_exact_dx = diff(u_exact, x_sym);
    L2_num = 0;
    L2_den = 0;
    H1_num = 0;
    H1_den = 0;
    n_int = pp + 2; % Ensure sufficient quadrature points
    [xi, weight] = Gauss(n_int, -1, 1);
    for ee = 1:n_el
        x_ele = x_coor(IEN(ee,:));
        u_ele = disp(IEN(ee,:));
        for qua = 1:n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            u_h = 0.0;
            u_h_xi = 0.0;
            for aa = 1:n_en
                N = PolyShape(pp, aa, xi(qua), 0);
                dN_dxi = PolyShape(pp, aa, xi(qua), 1);
                x_l = x_l + x_ele(aa) * N;
                dx_dxi = dx_dxi + x_ele(aa) * dN_dxi;
                u_h = u_h + u_ele(aa) * N;
                u_h_xi = u_h_xi + u_ele(aa) * dN_dxi;
            end
            dxi_dx = 1.0 / dx_dxi;
            u_h_x = u_h_xi * dxi_dx;
            % Compute exact solution and derivatives
            u_ex = double(subs(u_exact, x_sym, x_l));
            u_ex_x = double(subs(u_exact_dx, x_sym, x_l));
            % Compute errors
            L2_num = L2_num + weight(qua) * (u_h - u_ex)^2 * dx_dxi;
            L2_den = L2_den + weight(qua) * (u_ex)^2 * dx_dxi;
            H1_num = H1_num + weight(qua) * (u_h_x - u_ex_x)^2 * dx_dxi;
            H1_den = H1_den + weight(qua) * (u_ex_x)^2 * dx_dxi;
        end
    end
    e_L2 = sqrt(L2_num / L2_den);
    e_H1 = sqrt(H1_num / H1_den);
end

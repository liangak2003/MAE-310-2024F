function [L2_error, H1_error] = error_calc(x_coor, y_coor, disp, exact, exact_x, exact_y, n_el_x, n_el_y, IEN)
    L2_error_sq = 0;
    H1_error_sq = 0;
    
    n_el = size(IEN, 1); % Number of Elements
    
    % Loop Through Each Element
    for ee = 1:n_el
        nodes = IEN(ee, :);
        x_ele = x_coor(nodes);
        y_ele = y_coor(nodes);
        u_num_ele = disp(nodes);
        
        % Compute Jacobian Determinant and Element Area
        detJ = (x_ele(2)-x_ele(1))*(y_ele(3)-y_ele(1)) - (x_ele(3)-x_ele(1))*(y_ele(2)-y_ele(1));
        area = abs(detJ) / 2;
        
        % Barycentric Coordinates at Centroid
        L1 = 1/3;
        L2 = 1/3;
        L3 = 1 - L1 - L2;
        
        % Shape Functions at Centroid
        N = [Triangle(1, L1, L2, L3), Triangle(2, L1, L2, L3), Triangle(3, L1, L2, L3)];
        
        % Physical Coordinates at Centroid
        x_l = N(1)*x_ele(1) + N(2)*x_ele(2) + N(3)*x_ele(3);
        y_l = N(1)*y_ele(1) + N(2)*y_ele(2) + N(3)*y_ele(3);
        
        % Numerical and Exact Solutions
        u_num = N * u_num_ele;
        u_ex = exact(x_l, y_l);
        
        % Accumulate L2 Error
        L2_error_sq = L2_error_sq + (u_ex - u_num)^2 * area;
        
        % Compute Numerical Solution Gradient
        dN_dx = [(y_ele(2) - y_ele(3)), (y_ele(3) - y_ele(1)), (y_ele(1) - y_ele(2))] / detJ;
        dN_dy = [(x_ele(3) - x_ele(2)), (x_ele(1) - x_ele(3)), (x_ele(2) - x_ele(1))] / detJ;
        grad_u_num_x = dN_dx * u_num_ele;
        grad_u_num_y = dN_dy * u_num_ele;
        
        % Exact Gradient
        grad_u_ex_x = exact_x(x_l, y_l);
        grad_u_ex_y = exact_y(x_l, y_l);
        
        % Accumulate H1 Error
        H1_error_sq = H1_error_sq + ((grad_u_ex_x - grad_u_num_x)^2 + (grad_u_ex_y - grad_u_num_y)^2) * area;
    end
    
    % Compute Final Errors
    L2_error = sqrt(L2_error_sq);
    H1_error = sqrt(H1_error_sq);
end

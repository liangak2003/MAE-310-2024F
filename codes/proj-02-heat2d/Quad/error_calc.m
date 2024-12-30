function [L2_error, H1_error] = error_calc(x_coor, y_coor, disp, exact, exact_x, exact_y, n_el_x, n_el_y)

    n_np_x = n_el_x + 1;
    n_np_y = n_el_y + 1;
    n_np = n_np_x * n_np_y;
    
    % Compute exact solution at nodes
    u_exact = exact(x_coor, y_coor);
    
    % Compute numerical solution
    u_num = disp;
    
    % Compute L2 error
    L2_error = sqrt(sum((u_exact - u_num).^2) * (1/n_np));
    
    % Compute gradients of exact solution
    u_exact_x = exact_x(x_coor, y_coor);
    u_exact_y = exact_y(x_coor, y_coor);
    
    % Approximate gradients of numerical solution using finite differences
    u_num_x = zeros(n_np,1);
    u_num_y = zeros(n_np,1);
    
    hx = 1.0 / n_el_x;
    hy = 1.0 / n_el_y;
    
    for j = 1:n_np_y
        for i = 1:n_np_x
            index = (j-1)*n_np_x + i;
            % Compute u_num_x
            if i == 1
                u_num_x(index) = (u_num(index+1) - u_num(index)) / hx;
            elseif i == n_np_x
                u_num_x(index) = (u_num(index) - u_num(index-1)) / hx;
            else
                u_num_x(index) = (u_num(index+1) - u_num(index-1)) / (2*hx);
            end
            
            % Compute u_num_y
            if j == 1
                u_num_y(index) = (u_num(index+n_np_x) - u_num(index)) / hy;
            elseif j == n_np_y
                u_num_y(index) = (u_num(index) - u_num(index-n_np_x)) / hy;
            else
                u_num_y(index) = (u_num(index+n_np_x) - u_num(index-n_np_x)) / (2*hy);
            end
        end
    end
    
    % Compute H1 error
    H1_error = sqrt(sum((u_exact_x - u_num_x).^2 + (u_exact_y - u_num_y).^2) * (1/n_np));

end

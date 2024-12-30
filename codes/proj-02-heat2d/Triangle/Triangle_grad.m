function [dN_dL1, dN_dL2, dN_dL3] = Triangle_grad(aa)
    % aa: node number (1, 2, or 3)
    % For linear triangles, gradients are constants

    if aa == 1
        dN_dL1 = 1;
        dN_dL2 = 0;
        dN_dL3 = 0;
    elseif aa == 2
        dN_dL1 = 0;
        dN_dL2 = 1;
        dN_dL3 = 0;
    elseif aa == 3
        dN_dL1 = 0;
        dN_dL2 = 0;
        dN_dL3 = 1;
    else
        error('Error: value of aa should be 1, 2, or 3.');
    end
end

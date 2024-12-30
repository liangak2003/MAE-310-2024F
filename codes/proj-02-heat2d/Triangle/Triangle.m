function N = Triangle(aa, L1, L2, L3)
    % aa: node number (1, 2, or 3)
    % L1, L2, L3: barycentric coordinates

    if aa == 1
        N = L1;
    elseif aa == 2
        N = L2;
    elseif aa == 3
        N = L3;
    else
        error('Error: value of aa should be 1, 2, or 3.');
    end
end

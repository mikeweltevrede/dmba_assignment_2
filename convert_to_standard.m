function [f, Aeq, beq, lb, ub] = convert_to_standard(f, A, b, lb, ub, sign)
    
    % Make sure that the input parameter sign has value geq, leq, or eq.
    assert(any(strcmp(sign, ["leq", "geq", "eq"])), ...
        sprintf(['"%s" is not a valid input. Choose sign as '...
        '"leq", "geq", or "eq".'], sign))

    if sign == "leq"
        Aeq = cat(2, A, eye(size(A,1)));
    elseif sign == "geq"
        Aeq = cat(2, A, -eye(size(A,1)));
    elseif sign == "eq"
        Aeq = A;
    else
        disp("")
    end

    % Construct all other variables
    f = cat(2, f, zeros(1, size(A, 1)));
    beq = b;
    lb = zeros(1, size(Aeq, 2));
    ub = [];
end
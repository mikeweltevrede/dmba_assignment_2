function [f, Aeq, beq, lb, ub] = ...
    convert_to_standard(f, Aleq, bleq, Ageq, bgeq, Aeq, beq, lb, ub)
    % This function assumes that there is no upper bound to the variables
    % and that the lower bound for the slack variables is 0 (as is true by
    % definition).

    % Assert that f, lb, ub, and the b's are row vectors
    if size(f, 1) > 1
        f = f';
    end

    if size(lb, 1) > 1
        lb = lb';
    end
    
    if size(ub, 1) > 1
        ub = ub';
    end
    
    if size(bleq, 1) > 1
        bleq = bleq';
    end
    
    if size(bgeq, 1) > 1
        bgeq = bgeq';
    end
    
    if size(beq, 1) > 1
        beq = beq';
    end
    
    % Assert concurring shapes
    if ~isempty(Aleq) && ~isempty(Ageq)
        assert(size(Aleq, 2) == size(Ageq, 2), ...
            "Aleq and Ageq have different number of columns")
    end
    
    if ~isempty(Aleq) && ~isempty(Aeq)
        assert(size(Aleq, 2) == size(Aeq, 2), ...
            "Aleq and Aeq have different number of columns")
    end
    
    if ~isempty(Ageq) && ~isempty(Aeq)
        assert(size(Ageq, 2) == size(Aeq, 2), ...
            "Ageq and Aeq have different number of columns")
    end
    
    if ~isempty(Ageq)
        assert(size(Ageq, 2) == size(f, 2), ...
            "f does not have the same number of entries as coefficient matrices have columns")
    end
    
    if ~isempty(Ageq) && ~isempty(lb)
        % Since the assertions above have passed, we take one of the
        % matrices to represent all matrices
        assert(size(Ageq, 2) == size(lb, 2), ...
            "lb does not have the same number of entries as coefficient matrices have columns")
    end
    
    if ~isempty(Aleq) && ~isempty(bleq)
        assert(size(Aleq, 1) == size(bleq, 2), ...
            "Aleq and bleq have different number of rows vs. entries")
    end
    
    if ~isempty(Ageq) && ~isempty(bgeq)
        assert(size(Ageq, 1) == size(bgeq, 2), ...
            "Ageq and bgeq have different number of rows vs. entries")
    end
    
    if ~isempty(Aeq) && ~isempty(beq)
        assert(size(Aeq, 1) == size(beq, 2), ...
            "Aeq and beq have different number of rows vs. entries")
    end
    
    % Create coefficients for slack variables
    leq = eye(size(Aleq, 1)); 
    geq = -eye(size(Ageq, 1));
    slacks = blkdiag(leq, geq);
    slacks = cat(1, slacks, zeros(size(Aeq, 1), size(slacks, 2)));
    
    % Append all arrays
    Aeq = cat(2, cat(1, Aleq, Ageq, Aeq), slacks);
    beq = cat(2, bleq, bgeq, beq);
    
    % Modify f and lb with zeros
    f = cat(2, f, zeros(1, size(slacks, 2)));
    lb = cat(2, lb, zeros(1, size(slacks, 2)));
    disp(Aeq)
end
%% Homework 2 DMBA
% We specify an example for the ReDog problem
f = [2 1];
A = [3 2; 3 1; 1 2];
b = [80 50 60];

%% Exercise 1A
% We first add slack variables. We need to add zeroes to the cost vector
% and the identity matrix on the right side of our A matrix. Moreover, we
% need to specify bounds on these slack variables. By definition, slack
% variables are bounded from below by 0 and have no upper bound.
% For this, we use the function constructed to convert any problem into
% standard form for exercise 2A
[f, Aeq, beq, lb, ub] = convert_to_standard(f, A, b, lb, ub, "leq");

% Solve this problem for an optimal basis
[x,z] = linprog(-f, [], [], Aeq, beq, lb, ub);

bv = find(x); % Nonzero entries of x
nbv = find(~x); % Zero entries of x

B = Aeq(:, bv);
N = Aeq(:, nbv);

disp('Optimal basis matrix:')
disp(B)

%% Exercise 1B
B_inv = inv(B);

% Get rid of floating point errors around 0
% TODO: We now round this to take care of zeroes being stored as
% 0.0000000001 etcetera. How can we fix this?
close_to_zero = ismembertol(B_inv, 0, 10e-5);
B_inv(close_to_zero) = round(B_inv(close_to_zero));

I = eye(size(A, 1));
x_bv = x(bv);

% Initialise arrays to store bounds
lb_b = [];
ub_b = [];

for i = 1:size(A, 1)

    % Select the column for that epsilon
    B_inv_col = B_inv * I(:, i);
    
    eps = -x_bv ./ B_inv_col;

    lower =[];
    upper = [];
    for j = 1:size(x_bv, 1)  %whether it is upper or lower only depends on 
                                %B_inv_col and not on x_bv
        if B_inv_col(j) >= 0
            lower(:, j) = eps(j, :);
        else
            upper(:, j) = eps(j, :);
        end
    end

    % To remove the zeroes, we use find() 
    lower_bound = lower(:, find(lower));
    upper_bound = upper(:, find(upper));
    
    % If there is no upper bound, set it to Inf. Idem for lower bound and
    % -Inf
    if isempty(upper_bound)
        upper_bound = Inf;
    end
    if isempty(lower_bound)
        lower_bound = -Inf;
    end
    
    % Calculate lower and upper bounds of b and store these
    lb_b(i,:) = b(:,i) + max(lower_bound);
    ub_b(i,:) = b(:,i) + min(upper_bound);

end

bounds_b = [lb_b b' ub_b];
col_names = {'lower_bound', 'current_value', 'upper_bound'};
bounds_of_b = array2table(bounds_b, 'VariableNames', col_names);

disp(bounds_of_b)

%% Exercise 3

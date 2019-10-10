%% Homework 2 DMBA
% Note: this script runs clear at the top of the file
clear

% We specify an example for the ReDog problem
f = [2 1];
A = [3 2; 3 1; 1 2];
b = [80 50 60];
lb = zeros(1, size(A,2));
ub = [];
size_cvector = [1:size(f,2)];

%% Exercise 1A
% We first add slack variables. We need to add zeroes to the cost vector
% and the identity matrix on the right side of our A matrix. Moreover, we
% need to specify bounds on these slack variables. By definition, slack
% variables are bounded from below by 0 and have no upper bound.
% For this, we use the function constructed to convert any problem into
% standard form for exercise 2A: convert_to_standard. Since our current
% problem only has 'less than or equal to' constraints, we call the
% function as specified below.

[f, Aeq, beq, lb, ub] = ...
    convert_to_standard(f, A, b, [], [], [], [], lb, ub);

% Solve this problem for an optimal basis
[x,z] = linprog(-f, [], [], Aeq, beq, lb, ub,...
    optimoptions('linprog','Display','none')); % Don't print "Optimal solution found."

bv = find(x); % Nonzero entries of x
nbv = find(~x); % Zero entries of x

B = Aeq(:, bv);
N = Aeq(:, nbv);

disp("============== Exercise 1A) ==============")
disp("------ Identifying the optimal basis ------")
fprintf(2,'\n')

disp('Optimal basis matrix:')
disp(B)

%% Exercise 1B
B_inv = inv(B);
x_bv = x(bv);

% Get rid of floating point errors around 0
close_to_zero = ismembertol(B_inv, 0, 10e-5);
B_inv(close_to_zero) = round(B_inv(close_to_zero));

% Construct identity matrix to later retrieve unit vectors from
I = eye(size(A, 1));

% Initialise arrays to store bounds
lb_b = zeros(size(A, 1), 1);
ub_b = zeros(size(A, 1), 1);

for i = 1:size(A, 1)

    % Select the column for that epsilon
    B_inv_col = B_inv * I(:, i);
    
    eps = -x_bv ./ B_inv_col;

    % Preallocate vectors for upper and lower bounds
    lower_bound = [];
    upper_bound = [];
    
    for j = 1:size(x_bv, 1)  %whether it is upper or lower only depends on 
                                %B_inv_col and not on x_bv
        if B_inv_col(j) >= 0
            % If the scalar is nonnegative, then it provides a lower bound
            lower_bound = cat(2, lower_bound, eps(j, :));
        else
            % If negative, it gives us an upper bound.
            upper_bound = cat(2, upper_bound, eps(j, :));
        end
    end
    
    % If there is no lower bound, set it to -Inf. Idem for upper bound and
    % Inf.
    if isempty(lower_bound)
        lower_bound = -Inf;
    end
    
    if isempty(upper_bound)
        upper_bound = Inf;
    end
    
    % Calculate lower and upper bounds of b and store these
    lb_b(i, :) = b(:, i) + max(lower_bound);
    ub_b(i, :) = b(:, i) + min(upper_bound);

end

bounds_b = [lb_b b' ub_b];
col_names = {'lower_bound', 'current_value', 'upper_bound'};
bounds_of_b = array2table(bounds_b, 'VariableNames', col_names);

disp("=============== Exercise 1B) ================")
disp("---- Computing the lower and upper bounds ----")
disp("---- for the entries of the b vector such ----")
disp("--- that the optimal basis is not changed. ---")
fprintf(2, '\n')

disp(bounds_of_b)

%% Exercise 1C
total_B = [N B_inv]; %misschien is dit onnodig maar ik wist even niet hoe 
                     %het anders moest haha
B_inv_nbv = total_B(:, nbv);
coef_nbv_z = (f(bv) * B_inv_nbv - f(nbv))';

% Initialise arrays to store bounds
lb_c = zeros(size(f(bv), 2), 1);
ub_c = zeros(size(f(bv), 2), 1);

for i = 1:size(f(bv), 2)
    
    % Select the column for that epsilon
    B_inv_col = B_inv_nbv(i, :)';

    eps = -coef_nbv_z ./ B_inv_col;

    % Preallocate vectors for upper and lower bounds
    lower_bound = [];
    upper_bound = [];
    
    for j = 1:size(B_inv_nbv, 2)
        if B_inv_col(j) >= 0
            % If the scalar is nonnegative, then it provides a lower bound
            lower_bound = cat(2, lower_bound, eps(j, :));
        else
            % If negative, it gives us an upper bound.
            upper_bound = cat(2, upper_bound, eps(j, :));
        end
    end

    % If there is no upper bound, set it to Inf. Idem for lower bound and
    % -Inf.
    if isempty(upper_bound)
        upper_bound = Inf;
    end

    if isempty(lower_bound)
        lower_bound = -Inf;
    end

    % Calculate lower and upper bounds of b and store these.
    lb_c(i, :) = f(:, i) + max(lower_bound);
    ub_c(i, :) = f(:, i) + min(upper_bound);

end
    
bounds_c = [lb_c(size_cvector) f(size_cvector)' ub_c(size_cvector)]
col_names = {'lower_bound', 'current_value', 'upper_bound'};
bounds_of_c = array2table(bounds_c, 'VariableNames', col_names);

disp("=============== Exercise 1C) ================")
disp("---- Computing the lower and upper bounds ----")
disp("--- for the entries of the cost vector such ---")
disp("--- that the optimal basis is not changed. ---")
fprintf(2, '\n')

disp(bounds_of_c)

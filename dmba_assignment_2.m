%% Homework 2 DMBA
% Note: this script runs clear and clc at the top of the file
clear
clc
 
% We specify an example for the problem from page 255 from the book
% (Winston, W. L. (2004). Operations research applications and algorithms.
% Belmont (Calif.): Duxbury Press.), figure 17 (LINDO Output for HAL).
f = [600 1000 800 1300 -20];
Aleq = [2 3 2 4 -1; 1 0 1 0 0; 0 1 0 1 0; 1 1 0 0 0; 0 0 1 1 0; 0 0 0 0 1];
bleq = [0 800 1000 900 900 4000];
Ageq = [];
bgeq = [];
Aeq = [];
beq = [];
 
lb = zeros(1, size(f, 2));
ub = [];
 
original_number_of_variables = size(f, 2);
 
%% Exercise 1A
% We first add slack variables. We need to add zeroes to the cost vector
% and the identity matrix on the right side of our A matrix. Moreover, we
% need to specify bounds on these slack variables. By definition, slack
% variables are bounded from below by 0 and have no upper bound.
% For this, we use the function constructed to convert any problem into
% standard form for exercise 2A: convert_to_standard. Since our current
% problem only has 'less than or equal to' constraints, we call the
% function with empty arrays for Aqeg, bgeq, Aeq, and beq.
 
[f, Aeq, beq, lb, ub] = ...
    convert_to_standard(f, Aleq, bleq, Ageq, bgeq, Aeq, beq, lb, ub);
 
% Solve this problem for an optimal basis
[x,z] = linprog(-f, [], [], Aeq, beq, lb, ub,...
    optimoptions('linprog','Display','none'));
 
% Find the basic (nonzero in x) and nonbasic variables, along with their
% current values and costs.
bv = find(x);
nbv = find(~x);
x_bv = x(bv); % x(nbv) is not needed, so we do not create it
f_bv = f(bv);
f_nbv = f(nbv);
 
B = Aeq(:, bv);
N = Aeq(:, nbv);
 
disp("============== Exercise 1A) ===============")
disp("------ Identifying the optimal basis ------")
fprintf(2,'\n')
 
disp('Basic variables in optimal basis:')
disp(bv')

disp('Nonbasic variables in optimal basis:')
disp(nbv')
 
disp('Optimal basis matrix:')
disp(B)
 
%% Exercise 1B
B_inv = inv(B);
 
% Get rid of floating point errors around 0
close_to_zero = ismembertol(B_inv, 0, 10e-5);
B_inv(close_to_zero) = round(B_inv(close_to_zero));
 
% Construct identity matrix to later retrieve unit vectors from
I = eye(size(Aeq, 1));
 
% Initialise arrays to store bounds
lb_b = zeros(size(Aeq, 1), 1);
ub_b = zeros(size(Aeq, 1), 1);
 
for i = 1:size(Aeq, 1)
 
    % Select the column for that epsilon
    B_inv_col = B_inv * I(:, i);
    eps = -x_bv ./ B_inv_col;
 
    % Preallocate vectors for upper and lower bounds
    lower_bound = [];
    upper_bound = [];
    
    for j = 1:size(x_bv, 1)
        if B_inv_col(j) >= 0
            % If the scalar is nonnegative, then it provides a lower bound
            lower_bound = cat(2, lower_bound, eps(j, :));
        else
            % If negative, it gives us an upper bound
            upper_bound = cat(2, upper_bound, eps(j, :));
        end
    end
    
    % If there is no lower bound, set it to -Inf. Idem for upper bound and
    % Inf
    if isempty(lower_bound)
        lower_bound = -Inf;
    end
    
    if isempty(upper_bound)
        upper_bound = Inf;
    end
    
    % Calculate lower and upper bounds of b and store these
    lb_b(i, :) = beq(:, i) + max(lower_bound);
    ub_b(i, :) = beq(:, i) + min(upper_bound);
 
end
 
bounds_b = [lb_b beq' ub_b];
col_names = {'LowerBound', 'CurrentValue', 'UpperBound'};
bounds_of_b = array2table(bounds_b, 'VariableNames', col_names);
 
disp("=============== Exercise 1B) =================")
disp("---- Computing the lower and upper bounds ----")
disp("---- for the entries of the b vector such ----")
disp("--- that the optimal basis is not changed. ---")
fprintf(2, '\n')
 
disp(bounds_of_b)
 
%% Exercise 1C 
% Number of (non)basic variables
num_bv = size(bv, 1);
num_nbv = size(nbv, 1);
 
% Initialise arrays to store bounds
lb_c_bv = zeros(num_bv, 1);
ub_c_bv = zeros(num_bv, 1);
 
% The lower bounds for nonbasic variables are always -Inf
lb_c_nbv = -Inf * ones(num_nbv, 1);
ub_c_nbv = zeros(num_nbv, 1);
 
B_inv_aj = B_inv * N;
coef_nbv = f_bv*B_inv_aj - f_nbv;
 
% Find the bounds for the basic variables
for i = 1:num_bv
    
    % Select the column for that epsilon
    B_inv_aj_row = B_inv_aj(i,:);
    eps = -coef_nbv ./ B_inv_aj_row;

    % Preallocate vectors for upper and lower bounds
    lower_bound = [];
    upper_bound = [];
    
    for j = 1:size(B_inv_aj, 2)
        if B_inv_aj_row(j) >= 0
            % If the scalar is nonnegative, then it provides a lower bound
            lower_bound = cat(2, lower_bound, eps(:, j));
        else
            % If negative, it gives us an upper bound.
            upper_bound = cat(2, upper_bound, eps(:, j));
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
 
    % Calculate lower and upper bounds of b and store these.
    lb_c_bv(i, :) = f_bv(:, i) + max(lower_bound);
    ub_c_bv(i, :) = f_bv(:, i) + min(upper_bound);
 
end
 
bounds_c_bv = [lb_c_bv f_bv' ub_c_bv];
 
% Find the (upper) bounds for the nonbasic variables (recall that lower
% bounds for nonbasic variables are always -Inf)
for i = 1:num_nbv
    coef_nbv = f_bv*B_inv_aj(:,i) - f_nbv(:,i);
    ub_c_nbv(i, :) = f_nbv(:, i) + coef_nbv ;
end
 
bounds_c_nbv = [lb_c_nbv f_nbv' ub_c_nbv];
 
% Create a table with the bounds and the variable "names" (bv/nbv)
table = cat(2, cat(1, bv, nbv), cat(1, bounds_c_bv, bounds_c_nbv));
table = sortrows(table);
col_names = {'Variable', 'LowerBound', 'CurrentValue', 'UpperBound'};
 
% We only show the bounds for the nonslack variables, as these have cost 0
% by definition, so a changed cost for slack variables does not make sense
bounds_of_c = array2table(table(1:original_number_of_variables, :),...
    'VariableNames', col_names);
 
disp("================ Exercise 1C) =================")
disp("---- Computing the lower and upper bounds -----")
disp("--- for the entries of the cost vector such ---")
disp("--- that the optimal basis is not changed. ----")
fprintf(2, '\n')
 
disp(bounds_of_c)
 

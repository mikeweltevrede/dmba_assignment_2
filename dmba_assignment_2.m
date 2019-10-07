%% Homework 2 DMBA
%% Exercise 1

f = [2 1];
A = [3 2; 3 1; 1 2];
b = [80 50 60];

f = cat(2, f, zeros(1, size(A,1)))
Aeq = cat(2, A, eye(size(A,1)))
beq = b;
lb = zeros(1, size(Aeq, 2));
ub = [];

[x,z] = linprog(-f,[],[],Aeq,beq,lb,ub)

bv = find(x);
nbv = find(~x);

B = Aeq(:, bv)
N = Aeq(:, nbv)

%% Exercise 2
B_inv = inv(B);
I = eye(size(A,1));

lb_b = [];
ub_b = [];
for i = 1:size(A,1)     %%VRAAG::: vanaf i=2 loopt hij goed, maar voor i=1 
                        % twijfel ik, wat denk jij?

    B_inv_col = B_inv * I(:,i); %select the column for that epsilon
    x_bv = x(bv);
    eps = -x_bv./B_inv_col      %%VERVOLG OP VRAAG::: Want voor i=1 zegt hij 
                                % hier namelijk [-inf, 9.3675, -0.000], terwijl 
                                % het [inf, inf, 4] zou moeten zijn toch? Of zit 
                                % ik nou helemaal verkeerd?

    lower =[];
    upper = [];
    for j = 1:size(x_bv,1)  %whether it is upper or lower only depends on 
                            %B_inv_col and not on x_bv
        if B_inv_col(j) >= 0
            lower(:,j) = eps(j,:);
        else
            upper(:,j) = eps(j,:);
        end
    end

    lower1 = find(lower); % To remove the zeros 
    upper1 = find(upper);
    lower_bound = lower(:, lower1);
    upper_bound = upper(:, upper1);
  
    lb_b(i,:) = b(:,i) + max(lower_bound); %calculate lower and upper bounds of b
    ub_b(i,:) = b(:,i) + min(upper_bound);

end

bounds_b = [lb_b b' ub_b];
col_names = {'lower_bound', 'current_value', 'upper_bound'};
bounds_of_b = array2table(bounds_b, 'VariableNames', col_names)

%% Exercise 3

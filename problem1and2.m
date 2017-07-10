% UCLA EE 239
% Project 3 Collaborative Filtering

clc
clear

% Part 1
% Load Dataset
R_raw = importdata ('u.data');
user = R_raw(:,1);
item = R_raw(:,2);
rating = R_raw(:,3);
R_row = max(user);
R_column = max(item);
R = NaN(R_row, R_column);
W = zeros(R_row, R_column);
for i = 1:size(user)
    R(user(i),item(i)) = rating(i);
    W(user(i),item(i)) = 1;
end
Rpart1 = R;
% Least Square Factorization
k1 = 10;
[U1,V1] = wnmfrule(R,k1);
Rpred1 = U1*V1;
k2 = 50;
[U2,V2] = wnmfrule(R,k2);
Rpred2 = U2*V2;
k3 = 100;
[U3,V3] = wnmfrule(R,k3);
Rpred3 = U3*V3;

for i = 1:1:numel(R)
    if Rpred1(i) > 5
        Rpred1(i) = 5;
    end
    if Rpred2(i) > 5
        Rpred2(i) = 5;
    end
    if Rpred3(i) > 5
        Rpred3(i) = 5;
    end
    if isnan(Rpart1(i)) == 1
        Rpart1(i) = 0;
        Rpred1(i) = 0;
        Rpred2(i) = 0;
        Rpred3(i) = 0;
    end        
end
E1 = W.*((Rpart1-Rpred1).^2);
error1 = sum(E1(:));
E2 = W.*((Rpart1-Rpred2).^2);
error2 = sum(E2(:));
E3 = W.*((Rpart1-Rpred3).^2);
error3 = sum(E3(:));

% Part 2
% 10-fold Cross Validation
known_indices = find(isnan(R)==0);
N = 100000;
prm = randperm(N);
K = 10;
% K = 50;
% K = 100;
abs_error = zeros(1,10);
for i=1:1:10
    test_indices = known_indices(prm((i-1)*N/10+1:N/10*i));
    Rtrain = R;
    Rtrain(test_indices) = NaN;
    [U,V] = wnmfrule(Rtrain,K);
    Rprediction = U*V;
    Rtest = Rprediction(test_indices);
    for j = 1:1:numel(Rtest)
        if Rtest(j)> 5
            Rtest(j) = 5;
        end
    end
    abs_error(i) = sum(abs(Rtest-R(test_indices)))/(N/10);
end





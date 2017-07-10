clc
clear
k1 = 50;

R_raw = importdata ('u.data');
R_row = max(R_raw(:,1));
R_column = max(R_raw(:,2));
R = NaN(R_row,R_column);
R_2 = NaN(R_row,R_column,10);

predW = zeros(R_row,R_column);
predR = zeros(R_row,R_column);
for i = 1:100000;
        R(R_raw(i,1),R_raw(i,2)) = R_raw(i,3);
end
entryAvail = find(isnan(R)==0);
N = length(entryAvail);
Mat1 = randperm(N);
for i=1:10
     for j = 1:N
        R_2(R_raw(j,1),R_raw(j,2),i) = R_raw(j,3);
     end
     for k = 1:10
        R_3 = R_2;
        testset = entryAvail(Mat1(N/10*(k-1)+1:N/10*k));
        R_3(testset*i) = 0;
        R_3(R_3==0)=NaN;
     end
     [W(:,:,i),U(:,:,i),V(:,:,i),numIter,tElapsed,finalResidual] = wnmfrulep5(R_3(:,:,i),k1);    
      TempR(:,:,i) = U(:,:,i)*V(:,:,i);
     for k = 1:10
        testset = entryAvail(Mat1(N/10*(k-1)+1:N/10*k));
        predW(testset) = TempR(testset);
     end
end 
predR = predW.* R;
L1=5;
rec = zeros(L1,R_row);
maxr = zeros(L1,R_row);
count = 0;
for x = 1:L1
    [maxr(x,:), rec(x,:)] = max(predR,[],2);
    for i = 1:R_row
        if R(i,rec(x,i)) >= 4
           count = count + 1;
        end
    end
end
precision = count / ( L1*R_row);

L2 = 50;
counth = zeros(L2);
countf = zeros(L2);
rec = zeros(L1,R_row);
maxr = zeros(L1,R_row);
for L = 1:L2
    TempR = predR;
    for y = 1:L
        [maxr(y,:), rec(y,:)] = max(TempR,[],2);
        for i = 1:R_row
            if R(i,rec(y,i)) >= 4
                counth(L) = counth(L) + 1;
            elseif R(i,rec(y,i)) <=3
                countf(L) = countf(L) + 1;
            end
        end
    end
end
rateh = counth/100000;
ratef = countf/100000;

plot (ratef,rateh);
ylabel('hit rate');
xlabel('falsealarm rate');
title('hitrate vs falsealarm rate');


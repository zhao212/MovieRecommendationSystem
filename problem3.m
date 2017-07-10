%problem 3

clc
clear
close all


R_raw = importdata ('u.data');
user = R_raw(:,1);
mov = R_raw(:,2);
rate = R_raw(:,3);
R_row = max(user);
R_column = max(mov);
R = NaN(R_row,R_column);

for i = 1:size(user)
    R(user(i),mov(i)) = rate(i);
end
R_2 = R;
k1 = 10;
entryAvail = find(isnan(R)==0);
N = length(entryAvail);
prm = randperm(N);
testset = entryAvail(prm(1:floor(N/10)));
R_2(testset) = 0;
R_2(R_2==0)=NaN;
[U1,V1,numIter1,tElapsed1,finalResidual1] = wnmfrule(R_2,k1);
Rtemp1 = U1*V1;
x = 1:0.01:4;
j=1;
for i = 1:0.01:4
    precision1(j) = length(find((R(testset)>3)&Rtemp1(testset)>i))/length(find(Rtemp1(testset)>i));
    recall1(j) = length(find((R(testset)>3)&Rtemp1(testset)>i))/length(find(R(testset)>3));
    j = j+1;
end
figure
plot(x,precision1)
trapz(x,precision1)
xlabel('threshold');
ylabel('precision');
title('precision vs threshold k = 10');
figure
plot(x,recall1)
trapz(x,recall1)
xlabel('threshold');
ylabel('recall');
title('recall vs threshold k = 10');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2 = 50;
[U1,V1,numIter1,tElapsed1,finalResidual1] = wnmfrule(R_2,k2);
Rtemp1 = U1*V1;
j=1;
for i = 1:0.01:4
    precision1(j) = length(find((R(testset)>3)&Rtemp1(testset)>i))/length(find(Rtemp1(testset)>i));
    recall1(j) = length(find((R(testset)>3)&Rtemp1(testset)>i))/length(find(R(testset)>3));
    j = j+1;
end
figure
plot(x,precision1)
trapz(x,precision1)
xlabel('threshold');
ylabel('precision');
title('precision vs threshold k = 50');
figure
plot(x,recall1)
trapz(x,recall1)
xlabel('threshold');
ylabel('recall');
title('recall vs threshold k = 50');

k3 = 100;
[U1,V1,numIter1,tElapsed1,finalResidual1] = wnmfrule(R_2,k3);
Rtemp1 = U1*V1;
j=1;
for i = 1:0.01:4
    precision1(j) = length(find((R(testset)>3)&Rtemp1(testset)>i))/length(find(Rtemp1(testset)>i));
    recall1(j) = length(find((R(testset)>3)&Rtemp1(testset)>i))/length(find(R(testset)>3));
    j = j+1;
end
figure
plot(x,precision1)
trapz(x,precision1)
xlabel('threshold');
ylabel('precision');
title('precision vs threshold k = 100');
figure
plot(x,recall1)
trapz(x,recall1)
xlabel('threshold');
ylabel('recall');
title('recall vs threshold k = 100');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

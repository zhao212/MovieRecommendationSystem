%clc
clear
close all

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
end
W((R>0)) = 1;

Rpart1 = R;
% Least Square Factorization
k1 = 10;
[U1,V1] = wnmfrule40(R,k1);
Rpred1 = U1*V1;
k2 = 50;
[U2,V2] = wnmfrule40(R,k2);
Rpred2 = U2*V2;
k3 = 100;
[U3,V3] = wnmfrule40(R,k3);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Comparing with part 1, take k = 100
% %lamda = 0.01
% [U41,V41] = wnmfrulep41 (0.01,R,100);
% E41 = W.*((R-U41*V41).^2);
% error41 = sum(E41(:));
% 
% %lamda = 0.1
% [U42,V42] = wnmfrulep41 (0.1,R,100);
% E42 = W.*((R-U42*V42).^2);
% error42 = sum(E42(:));
% 
% %lamda = 1
% [U43,V43] = wnmfrulep41 (1,R,100);
% E43 = W.*((R-U43*V43).^2);
% error43 = sum(E43(:));
% 
% %Compared with previous part, take k =10;
% %lamda = 0.01
% [U44,V44] = wnmfrulep42 (0.01,R,10);
% E44 = R.*((W-U44*V44).^2);
% error44 = sum(E44(:));
% 
% %lamda = 0.1
% [U45,V45] = wnmfrulep42 (0.1,R,10);
% E45 = R.*((W-U45*V45).^2);
% error45 = sum(E45(:));
% 
% %lamda = 1
% [U46,V46] = wnmfrulep42 (1,R,10);
% E46 = R.*((W-U46*V46).^2);
% error46 = sum(E46(:));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R_2 = R;
% k1 = 10;
% known_indices = find(isnan(R)==0);
% N = length(known_indices);
% prm = randperm(N);
% test_indices = known_indices(prm(1:floor(N/10)));
% R_2(test_indices) = 0;
% R_2(R_2==0)=NaN;
% [U2,V2] = wnmfrulep43(R_2,k1);
% 
% %Rtemp1 = U1*V1;
% Rtemp2 = U2*V2;
% Rtemp3 = Rtemp2 .* R;
% j=1;
% for i = 1:0.01:4
%     precision2(j) = length(find((R(test_indices)>i)&Rtemp3(test_indices)>i))/length(find(Rtemp3(test_indices)>i));
%     recall2(j) = length(find((R(test_indices)>i)&Rtemp3(test_indices)>i))/length(find(R(test_indices)>i));
%     j = j+1;
% end
% x = 1:0.01:4;
% figure
% plot(x,precision2)
% trapz(x,precision2)
% xlabel('threshold');
% ylabel('precision');
% title('precision vs threshold k = 10');
% figure
% plot(x,recall2)
% trapz(x,recall2)
% xlabel('threshold');
% ylabel('recall');
% title('recall vs threshold k = 10');
% 
% k1 = 50;
% [U2,V2] = wnmfrulep43(R_2,k1);
% 
% %Rtemp1 = U1*V1;
% Rtemp2 = U2*V2;
% Rtemp3 = Rtemp2 .* R;
% j=1;
% for i = 1:0.01:4
%     precision2(j) = length(find((R(test_indices)>i)&Rtemp3(test_indices)>i))/length(find(Rtemp3(test_indices)>i));
%     recall2(j) = length(find((R(test_indices)>i)&Rtemp3(test_indices)>i))/length(find(R(test_indices)>i));
%     j = j+1;
% end
% x = 1:0.01:4;
% figure
% plot(x,precision2)
% trapz(x,precision2)
% xlabel('threshold');
% ylabel('precision');
% title('precision vs threshold k = 50');
% figure
% plot(x,recall2)
% trapz(x,recall2)
% xlabel('threshold');
% ylabel('recall');
% title('recall vs threshold k = 50');
% 
% k1 = 50;
% [U2,V2] = wnmfrulep43(R_2,k1);
% 
% %Rtemp1 = U1*V1;
% Rtemp2 = U2*V2;
% Rtemp3 = Rtemp2 .* R;
% j=1;
% for i = 1:0.01:4
%     precision2(j) = length(find((R(test_indices)>i)&Rtemp3(test_indices)>i))/length(find(Rtemp3(test_indices)>i));
%     recall2(j) = length(find((R(test_indices)>i)&Rtemp3(test_indices)>i))/length(find(R(test_indices)>i));
%     j = j+1;
% end
% x = 1:0.01:4;
% figure
% plot(x,precision2)
% trapz(x,precision2)
% xlabel('threshold');
% ylabel('precision');
% title('precision vs threshold k = 100');
% figure
% plot(x,recall2)
% trapz(x,recall2)
% xlabel('threshold');
% ylabel('recall');
% title('recall vs threshold k = 100');
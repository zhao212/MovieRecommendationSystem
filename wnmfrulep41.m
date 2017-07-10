function [A,Y,numIter,tElapsed,finalResidual]=wnmfrulep41(lamda,X,k,option)
% Sparse NMF, X=AY, using multiplicative update rules. Both basis matrix and coefficient matrix are constrained by
% l_2 and l_1 norm.
% Definition:
%     [A,Y,numIter,tElapsed,finalResidual]=sparsenmf2rule(X,k)
%     [A,Y,numIter,tElapsed,finalResidual]=sparsenmf2rule(X,k,option)
% X: matrix, dataset to factorize, each column is a sample, and each row is a feature.
% k: number of clusters.
% option: struct:
% option.alpha2: non-negative scalar, control the smoothness and scale of
% the basis vectors. The default is 0.
% option.alpha1: non-negative scalar, control the sparsity of the basis
% vectors. The default is 0.
% option.lambda2: non-negative scalar, control the smoothness and scale of
% the coefficient vectors. The default is 0.
% option.lambda1: non-negative scalar, control the sparsity of the
% coefficient vectors. The default is 0.
% option.iter: max number of interations. The default is 1000.
% option.dis: boolen scalar, It could be 
%     false: not display information,
%     true: display (default).
% option.residual: the threshold of the fitting residual to terminate. 
%    If the ||X-XfitThis||<=option.residual, then halt. The default is 1e-4.
% option.tof: if ||XfitPrevious-XfitThis||<=option.tof, then halt. The default is 1e-4.
% A: matrix, the basis matrix.
% Y: matrix, the coefficient matrix.
% numIter: scalar, the number of iterations.
% tElapsed: scalar, the computing time used.
% finalResidual: scalar, the fitting residual.
% References:
%  []\
%%%%
% Copyright (C) <2012>  <Yifeng Li>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% Feb. 26, 2013
%%%%

tStart=tic;
optionDefault.alpha2=1;
optionDefault.alpha1=0.01;
optionDefault.lambda2=0;
optionDefault.lambda1=0.01;
optionDefault.iter=1000;
optionDefault.dis=true;
optionDefault.residual=1e-4;
optionDefault.tof=1e-4;
if nargin<4
   option=optionDefault;
else
    option=mergeOption(option,optionDefault);
end

[numFeat,numS]=size(X);

%compute the weight matrix
W=zeros(length(X(:,1)),length(X(1,:)));
W((X>0))=1;

W=sqrt(W);


% iter: number of iterations
[r,c]=size(X); % c is # of samples, r is # of features
Y=rand(k,c);
% Y(Y<eps)=0;
Y=max(Y,eps);
A=X/Y;
% A(A<eps)=0;
A=max(A,eps);
XfitPrevious=Inf;
% E1=option.alpha1.*ones(numFeat,k);
% E2=option.lambda1.*ones(k,numS);
ep=0;
for i=1:option.iter
    %very important
    A=A.*(((W.*X)*Y')./(W.*(A*Y)*Y'+lamda*A));
    %             A(A<eps)=0;
    A=max(A,ep);
    Y=Y.*((A'*(W.*X))./(A'*(W.*(A*Y))+lamda*Y));
    %             Y(Y<eps)=0;
    Y=max(Y,ep);
    if mod(i,100)==0 || i==option.iter
        if option.dis
            disp(['Iterating >>>>>> ', num2str(i),'th']);
        end
        % very important
        XfitThis=A*Y;
        fitRes1=matrixNorm(W.*(XfitPrevious-XfitThis));
        fitRes = sqrt((fitRes1.^2)+lamda*(sum(A(:).^2)+sum(Y(:).^2)));
        XfitPrevious=XfitThis;
        temp = norm(W.*(X-XfitThis),'fro');
        temp2= temp^2 + lamda*(sum(A(:).^2)+sum(Y(:).^2));
        curRes = sqrt(temp2);
        if option.tof>=fitRes || option.residual>=curRes || i==option.iter
            s=sprintf('Mutiple update rules based Sparse NMF on Both Factors successes! \n # of iterations is %0.0d. \n The final residual is %0.4d.',i,curRes);
            disp(s);
            numIter=i;
            finalResidual=curRes;
            break;
        end
    end
end
tElapsed=toc(tStart);
end

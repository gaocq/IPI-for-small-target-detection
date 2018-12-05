function [A_hat E_hat] = winRPCA_median(I, opt)
% This code is the core part of the small target detecton algorithm in our 
% published paper: Chenqiang Gao, Deyu Meng, Yi Yang, et al., "Infrared Patch-Image Model for 
% Small Target Detection in a Single Image, " Image Processing, IEEE Transactions on, vol. 22, no. 12,
% pp. 4996-5009, 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code in your publications, please cite:
% @article{Gao2013Infrared,
%    author = {Gao, Chenqiang and Meng, Deyu and Yang, Yi and Wang, Yongtao and Zhou, Xiaofang and Hauptmann, Alex},
%    title = {Infrared Patch-Image Model for Small Target Detection in a Single Image},
%    journal = {Image Processing, IEEE Transactions on},
%    volume = {22},
%    number = {12},
%    pages = {4996-5009},
%    year = {2013}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paremeter Explanation
%
% I - m x n matrix of an infrared image
% opt - a struct with at least four following members:
%     dw, dh: the width and height of the patch
%     x_step, y_step: the sliding steps of patch along x-axis and y-axis            
%
% [A_hat, E_hat] - estimates for background image and target image, respectively
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Chenqiang Gao 
% Email: gaochenqiang@gmail.com or gaocq@cqupt.edu.cn
% Copyright:  Chongqing University of Posts and Telecommunications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(I, 3) == 3;
    I = rgb2gray(I);
end
[m n] = size(I);
A_hat = zeros(m,n);
E_hat = zeros(m ,n);
C = zeros(m, n);
dw = opt.dw;
dh = opt.dh;
x_step = opt.x_step;
y_step = opt.y_step;

I1 = I ; % imfilter(I, fspecial('gaussian', 5));
D = [];
for i = 1:y_step:m-dh+1
    for j = 1:x_step:n-dw+1
        temp = I1(i:i+dh-1, j:j+dw-1);
        D = [D, temp(:)];
    end
end
D = mat2gray(D);
% D = D - mean2(D);
[m1 n1] = size(D);
lambda = 1 / sqrt(max(m, n));

[A1, E1] = APG_IR(D, lambda);

AA = zeros(m, n, 100);
EE = zeros(m, n, 100);
temp = zeros(dh, dw);
temp1 = zeros(dh, dw);
index = 0;
for i = 1:y_step:m-dh+1
    for j = 1:x_step:n-dw+1
        index = 1+index;
        temp(:) = A1(:, index);
        C(i:i+dh-1, j:j+dw-1) = C(i:i+dh-1, j:j+dw-1)+1;
        temp1(:) = E1(:, index);
        for ii = i:i+dh-1
            for jj = j:j+dw-1
                AA(ii,jj, C(ii,jj)) = temp(ii-i+1, jj-j+1);
                EE(ii,jj, C(ii,jj)) = temp1(ii-i+1, jj-j+1);
            end
        end
    end
end
%     C(find(C==0)) = 1000;
for i=1:m
    for j=1:n
        if C(i,j) > 0
            A_hat(i,j) = median(AA(i,j,1:C(i,j)));
            E_hat(i,j) = median(EE(i,j,1:C(i,j)));
        end
    end
end


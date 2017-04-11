close all;
clear all;
clc

J = imread('ski_image.jpg');
I = imresize(J, 0.5);
[row , col , depth] = size(I);
x1 = reshape(I(:,:,1), 1, row*col);
x2 = reshape(I(:,:,2), 1, row*col);
x3 = reshape(I(:,:,3), 1, row*col);
X = [x1; x2; x3];
X = double(X);
X = X/255;
mu1 = [0.47; 0.47; 0.47];
mu2 = [0.05; 0.05; 0.05];
mu3 = [0.7; 0.7; 0.7];
mu = [mu1 mu2 mu3];
Kmu = mu;
C1 = eye(3);
C2 = eye(3);
C3 = eye(3);
C = [C1 C2 C3];
w1 = 1/3;
w2 = 1/3;
w3 = 1/3;
w = [w1 w2 w3];
resp = zeros(row*col,3);
clus = zeros(row*col,3);
for m = 1:70
for i = 1:row*col
    for k = 1:3
        resp(i,k) = w(k)*gaussian(X(:,i),mu(:,k),C(:,3*(k-1)+1:3*k));
        dis(k) = sum((X(:,i)-Kmu(:,k)).^2);
    end
    resp(i,:) = resp(i,:)/sum(resp(i,:));
    [val, ind] = min(dis);
    clus(i,ind) = 1;
end
clus_im = reshape(clus, [row , col , 3]);
img = imresize(clus_im, 2);
figure(1)
imshow(img)
for k =1:3
    Kmu(:,k) = X*clus(:,k)/length(find(clus(:,k)>0));
end
n = sum(resp);
N = [n; n; n];
mu = (X*resp)./N; 
c = zeros(3,9);
for k =1:3
    for i = 1:row*col
        c(:,3*(k-1)+1:3*k) = c(:,3*(k-1)+1:3*k) + (resp(i,k)*(X(:,i)-mu(:,k))*(X(:,i)-mu(:,k))')./n(k);
    end
end
C = c;
w = n/(row*col);
l=0;
for i = 1:row*col   
        l = l + log(w*[gaussian(X(:,i),mu(:,1),C(:,1:3));gaussian(X(:,i),mu(:,2),C(:,4:6));gaussian(X(:,i),mu(:,3),C(:,7:9))]);
end
figure(2)
plot(m,l,'rx')
hold on
end

fprintf('For first cluster, the final values are\n');
fprintf('Mean\n');
fprintf('%f\n',mu(:,1));
fprintf('Covariance Matrix\n');
fprintf('%f\t',C(1,1:3));
fprintf('\n');
fprintf('%f\t',C(2,1:3));
fprintf('\n');
fprintf('%f\t',C(3,1:3));
fprintf('\n');
fprintf('Prior Weight\n%f\n',w(1));

fprintf('For second cluster, the final values are\n');
fprintf('Mean\n');
fprintf('%f\n',mu(:,2));
fprintf('Covariance Matrix\n');
fprintf('%f\t',C(1,4:6));
fprintf('\n');
fprintf('%f\t',C(2,4:6));
fprintf('\n');
fprintf('%f\t',C(3,4:6));
fprintf('\n');
fprintf('Prior Weight\n%f\n',w(2));

fprintf('For third cluster, the final values are\n');
fprintf('Mean\n');
fprintf('%f\n',mu(:,3));
fprintf('Covariance Matrix\n');
fprintf('%f\t',C(1,7:9));
fprintf('\n');
fprintf('%f\t',C(2,7:9));
fprintf('\n');
fprintf('%f\t',C(3,7:9));
fprintf('\n');
fprintf('Prior Weight\n%f\n',w(3));

close all;
clear all;
clc;

X_full = zeros(200,10304);
X_test = zeros(200,10304);
image = zeros(28*10, 23*20);
W = zeros(10304,199);

dir_list = dir('./gallery');
dir_list_test = dir('./Probe');
for i = 3:42
    img_list = dir(strcat('./gallery/' , dir_list(i).name));
    img_list_test = dir(strcat('./Probe/' , dir_list_test(i).name));
    for j = 3:7
        A = imread(strcat('./gallery/' , dir_list(i).name , '/' , img_list(j).name));
        A_test = imread(strcat('./Probe/' , dir_list_test(i).name , '/' , img_list_test(j).name));
        X_full((i-3)*5+(j-2),:) = reshape(A,1,112*92);
        X_test((i-3)*5+(j-2),:) = reshape(A_test,1,112*92);
    end
end

Mean = mean(X_full);
X = X_full - (ones(size(X_full))*diag(Mean));
[evec, eval] = eigs(X*X'/200 , 199);


for i = 1:199
    W(:,i) = X'*evec(:,i);    
    W(:,i) = W(:,i)/norm(W(:,i));
end

z = X*W(:,1:5);

i=1;
for r = 1:28:28*10
    for c = 1:23:23*20
        if i<201
   b = z(i,:)*W(:,1:5)';
   a = reshape(b, [112 92]);
   p = imresize(a,0.25);
   image(r:r+27,c:c+22) = p;
   i = i+1;
        end
    end
end
figure(1)
imshow(image)

v = zeros(1,10304);
eval1 = eigs(X'*X/200, 10304);
eval1_sum = sum(eval1);
for i = 1:10304
    v(i) = (sum(eval1(1:i))/eval1_sum)*100;
end
figure(2)
plot(1:10304,v,'rx')

im = imread('face_input_1.pgm');
im = double(im);
x1 = reshape(im, [1 10304]);
x2 = x1 - Mean;
z1 = x2*W(:,1:15);
b = z1*W(:,1:15)';
a = reshape(b, [112 92]);
figure(3)
imshow(a)
imwrite(a, 'recon_face_input_1.pgm');

X_test = double(X_test);
Mean_test = mean(X_test);
X1 = X_test - (ones(size(X_test))*diag(Mean_test));

z2 = X*W(:,1:25);
z_test = X1*W(:,1:25);

dis = zeros(200,200);
for i = 1:200
    for j = 1:200
        dis(i,j) = (sum((z_test(i,:) - z2(j,:)).^2))^0.5;
    end
end

least_dis = ones(200,3,2).*10000;
for all=1:200
    for i=1:200
        j = 1;
        flag = 0;
        while(j<4 && flag==0)
            if( dis(all,i)<least_dis(all,j,2) )
                val = dis(all,i);
                val_ind = i;
                k = j;
                while(k<4)
                    temp = least_dis(all,k,2);
                    temp_ind = least_dis(all,k,1);
                    least_dis(all,k,2) = val;
                    least_dis(all,k,1) = val_ind;
                    val = temp;
                    val_ind = temp_ind;
                    k = k+1;
                end
                flag = 1;
            end
            j = j+1;
        end
    end
end
least_dis(:,:,1) = floor((least_dis(:,:,1)-1)/5) +1;

accuracy = 0;
dis_sort = sort(least_dis,2);
class = zeros(200,2);
for i = 1:200
    class(i,1)=dis_sort(i,1,1);
    if class(i,1) ~= dis_sort(i,2,1)
        class(i,1) = dis_sort(i,2,1);
        if class(i,1) ~= dis_sort(i,3,1)
            class(i,1) = least_dis(i,1,1);
        end
    end
    class(i,2) = floor((i-1)/5) + 1;
    if class(i,1) == class(i,2)
        accuracy = accuracy+1;
    end
end          
per_accuracy = accuracy/2;
fprintf('the percenteage accuracy obtained by PCA is %f%%\n',per_accuracy);
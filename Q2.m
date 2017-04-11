close all;
clc;

Mean_class = zeros(40,199);
SW = zeros(199,199);
Z = zeros(5,199);

z_fis = X*W;
z_fis_test = X1*W;

Mean_fis = mean(z_fis);

for i = 1:40
    z_sum = zeros(1,199);
    for j = 1:5
        z_sum = z_sum + z_fis((i-1)*5 + j,:);
    end
    Mean_class(i,:) = z_sum / 5;
end

M = Mean_class - (ones(size(Mean_class))*diag(Mean_fis));

SB = 5*(M')*M;

for i = 1:40
    for j = 1:5
        Z(j,:) = z_fis((i-1)*5 + j,:) - Mean_class(i,:);
    end
    SW = SW + Z' * Z;
end

[fis_vec , fis_val] = eigs(pinv(SW)*SB , 39);

W_fis = zeros(199,39);
for i = 1:39
    W_fis(:,i) = fis_vec(:,i);    
    W_fis(:,i) = W_fis(:,i)/norm(W_fis(:,i));
end
z_ftrain = z_fis - ones(size(z_fis))*diag(Mean_fis);
z_ftest = z_fis_test - ones(size(z_fis_test))*diag(mean(z_fis_test));

traind = z_ftrain*W_fis;
testd = z_ftest*W_fis;

fdis = zeros(200,200);
for i = 1:200
    for j = 1:200
        fdis(i,j) = (sum((testd(i,:) - traind(j,:)).^2))^0.5;
    end
end

least_fdis = ones(200,3,2).*10000;
for all=1:200
    for i=1:200
        j = 1;
        flag = 0;
        while(j<4 && flag==0)
            if( fdis(all,i)<least_fdis(all,j,2) )
                val = fdis(all,i);
                val_ind = i;
                k = j;
                while(k<4)
                    temp = least_fdis(all,k,2);
                    temp_ind = least_fdis(all,k,1);
                    least_fdis(all,k,2) = val;
                    least_fdis(all,k,1) = val_ind;
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
least_fdis(:,:,1) = floor((least_fdis(:,:,1)-1)/5) +1;

faccuracy = 0;
fdis_sort = sort(least_fdis,2);
fclass = zeros(200,2);
for i = 1:200
    fclass(i,1)=fdis_sort(i,1,1);
    if fclass(i,1) ~= fdis_sort(i,2,1)
        fclass(i,1) = fdis_sort(i,2,1);
        if fclass(i,1) ~= fdis_sort(i,3,1)
            fclass(i,1) = least_fdis(i,1,1);
        end
    end
    fclass(i,2) = floor((i-1)/5) + 1;
    if fclass(i,1) == fclass(i,2)
        faccuracy = faccuracy+1;
    end
end          
per_faccuracy = faccuracy/2;
fprintf('the percenteage accuracy obtained by FDA is %f%%\n',per_faccuracy);
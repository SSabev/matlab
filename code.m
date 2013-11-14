clear all;
clc()

load data_90.mat

[B, mean1, mean2] = init_values(data_90);

data_size=size(data_90(:,1));

means = [mean1; mean2];
[clusters2, centres2, sum_squared_func2] = mykmeans(data_90, means);

% figure
% scatter3(data_90(:,1),data_90(:,2),data_90(:,3), 100, clusters2, 'filled')
% xlabel SL, ylabel SW, zlabel PL

[mean3] = NM(data_90,data_size, means);
means = [means; mean3];
[clusters3, centres3, sum_squared_func3] = mykmeans(data_90, means);


load true_90.mat;

true = 0;
false = 0;

for i=1:data_size(1)
    if clusters3(i) == 2
        clusters3(i) = 3;
    elseif clusters3(i) == 3
        clusters3(i) = 2;
    end
end

for i=1:data_size(1)
    if clusters3(i) ~= true_90(i)
        false = false+1;
    else
        true = true+1;
    end
end

true
false


mycov1 = zeros(3,3);
mycov2 = zeros(3,3);
mycov3 = zeros(3,3);

sum1 = sum(clusters3 == 1)
sum2 = sum(clusters3 == 2)
sum3 = sum(clusters3 == 3)

for i=1:data_size
    if clusters3(i) == 1
        mycov1 = mycov1 + (data_90(i,:) - centres3(1,:))'*(data_90(i,:) - centres3(1,:));
    elseif clusters3(i) == 2
        mycov2 = mycov2 + (data_90(i,:) - centres3(3,:))'*(data_90(i,:) - centres3(3,:));      
    elseif clusters3(i) == 3
        mycov3 = mycov3 + (data_90(i,:) - centres3(2,:))'*(data_90(i,:) - centres3(2,:));
    end
end


mycov1 = mycov1./sum1
mycov2 = mycov2./sum3
mycov3 = mycov3./sum2

% figure
% scatter3(data_90(:,1),data_90(:,2),data_90(:,3), 100, true_90, 'filled')
% xlabel SL, ylabel SW, zlabel PL
% 
% figure
% scatter3(data_90(:,1),data_90(:,2),data_90(:,3), 100, clusters3, 'filled')
% xlabel SL, ylabel SW, zlabel PL

[mean4] = NM(data_90, data_size, means);
means = [means; mean4];
[clusters4, centres4, sum_squared_func4] = mykmeans(data_90, means);

[mean5] = NM(data_90, data_size, means);
means = [means; mean5];
[clusters5, centres5, sum_squared_func5] = mykmeans(data_90, means);


% figure
% subplot(2,2,1); plot(sum_squared_func2)
% title('Sum-squared error function for 2 clusters')
% subplot(2,2,2); plot(sum_squared_func3)
% title('Sum-squared error function for 3 clusters')
% subplot(2,2,3); plot(sum_squared_func4)
% title('Sum-squared error function for 4 clusters')
% subplot(2,2,4); plot(sum_squared_func5)
% title('Sum-squared error function for 5 clusters')

load data_900.mat;

[Distances, mean1, mean2] = init_values(data_900);

data_size=size(data_900(:,1));

gclusters = zeros(data_size(1),1);

for i=1:data_size(1)
    prob1 = (1/(2*pi)^1.5)*(det(mycov1)^1.5)*exp((-0.5)*(data_900(i,:)-centres3(1,:))*mycov1^(-1)*(data_900(i,:)-centres3(1,:))');
    prob2 = (1/(2*pi)^1.5)*(det(mycov3)^1.5)*exp((-0.5)*(data_900(i,:)-centres3(3,:))*mycov3^(-1)*(data_900(i,:)-centres3(3,:))');
    prob3 = (1/(2*pi)^1.5)*(det(mycov2)^1.5)*exp((-0.5)*(data_900(i,:)-centres3(2,:))*mycov2^(-1)*(data_900(i,:)-centres3(2,:))');
    
    if prob1 == max([prob1, prob2, prob3])
        gclusters(i) = 1;
    elseif prob2 == max([prob1, prob2, prob3])
        gclusters(i) = 2;
    elseif prob3 == max([prob1, prob2, prob3])
        gclusters(i) = 3;
    end
end


confusion_matrix_gaussian = zeros(3,3);
load true_900.mat;

gclusters

for i=1:data_size(1)
    if true_900(i) == gclusters(i) && true_900(i)==1
        confusion_matrix_gaussian(1,1) = confusion_matrix_gaussian(1,1) + 1;
    elseif true_900(i) == gclusters(i) && true_900(i)==2;
        confusion_matrix_gaussian(2,2) = confusion_matrix_gaussian(2,2) + 1;
    elseif true_900(i) == gclusters(i) && true_900(i)==3;
        confusion_matrix_gaussian(3,3) = confusion_matrix_gaussian(3,3) + 1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 1 && gclusters(i) == 2
        confusion_matrix_gaussian(2,1) = confusion_matrix_gaussian(2,1) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 1 && gclusters(i) == 3
        confusion_matrix_gaussian(3,1) = confusion_matrix_gaussian(3,1) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 2 && gclusters(i) == 1
        confusion_matrix_gaussian(1,2) = confusion_matrix_gaussian(1,2) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 2 && gclusters(i) == 3
        confusion_matrix_gaussian(3,2) = confusion_matrix_gaussian(3,2) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 3 && gclusters(i) == 1
        confusion_matrix_gaussian(1,3) = confusion_matrix_gaussian(3,1) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 3 && gclusters(i) == 2
        confusion_matrix_gaussian(2,3) = confusion_matrix_gaussian(2,3) +1;
    end
end

confusion_matrix_gaussian

means = [mean1; mean2];
[mean3] = NM(data_900,data_size, means);
means = [means; mean3];

[clusters3, centres3, sum_squared_func3_2] = mykmeans(data_900, means);


for i=1:data_size(1)
    if true_900(i) == gclusters(i) && true_900(i)==1
        confusion_matrix_gaussian(1,1) = confusion_matrix_gaussian(1,1) + 1;
    elseif true_900(i) == gclusters(i) && true_900(i)==2;
        confusion_matrix_gaussian(2,2) = confusion_matrix_gaussian(2,2) + 1;
    elseif true_900(i) == gclusters(i) && true_900(i)==3;
        confusion_matrix_gaussian(3,3) = confusion_matrix_gaussian(3,3) + 1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 1 && gclusters(i) == 2
        confusion_matrix_gaussian(2,1) = confusion_matrix_gaussian(2,1) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 1 && gclusters(i) == 3
        confusion_matrix_gaussian(3,1) = confusion_matrix_gaussian(3,1) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 2 && gclusters(i) == 1
        confusion_matrix_gaussian(1,2) = confusion_matrix_gaussian(1,2) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 2 && gclusters(i) == 3
        confusion_matrix_gaussian(3,2) = confusion_matrix_gaussian(3,2) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 3 && gclusters(i) == 1
        confusion_matrix_gaussian(1,3) = confusion_matrix_gaussian(1,3) +1;
    elseif true_900(i)~=gclusters(i) && true_900(i) == 3 && gclusters(i) == 2
        confusion_matrix_gaussian(2,3) = confusion_matrix_gaussian(2,3) +1;
    end
end



for i=1:data_size(1)
    if clusters3(i) == 2
        clusters3(i) = 3;
    elseif clusters3(i) == 3
        clusters3(i) = 2;
    end
end

confusion_matrix_kmeans = zeros(3,3);
for i=1:data_size(1)
    if true_900(i) == clusters3(i) && true_900(i)==1
        confusion_matrix_kmeans(1,1) = confusion_matrix_kmeans(1,1) + 1;
    elseif true_900(i) == clusters3(i) && true_900(i)==2;
        confusion_matrix_kmeans(2,2) = confusion_matrix_kmeans(2,2) + 1;
    elseif true_900(i) == clusters3(i) && true_900(i)==3;
        confusion_matrix_kmeans(3,3) = confusion_matrix_kmeans(3,3) + 1;
    elseif true_900(i)~=clusters3(i) && true_900(i) == 1 && clusters3(i) == 2
        confusion_matrix_kmeans(1,2) = confusion_matrix_kmeans(1,2) +1;
    elseif true_900(i)~=clusters3(i) && true_900(i) == 1 && clusters3(i) == 3
        confusion_matrix_kmeans(1,3) = confusion_matrix_kmeans(1,3) +1;
    elseif true_900(i)~=clusters3(i) && true_900(i) == 2 && clusters3(i) == 1
        confusion_matrix_kmeans(2,1) = confusion_matrix_kmeans(2,1) +1;
    elseif true_900(i)~=clusters3(i) && true_900(i) == 2 && clusters3(i) == 3
        confusion_matrix_kmeans(2,3) = confusion_matrix_kmeans(2,3) +1;
    elseif true_900(i)~=clusters3(i) && true_900(i) == 3 && clusters3(i) == 1
        confusion_matrix_kmeans(3,1) = confusion_matrix_kmeans(3,1) +1;
    elseif true_900(i)~=clusters3(i) && true_900(i) == 3 && clusters3(i) == 2
        confusion_matrix_kmeans(3,2) = confusion_matrix_kmeans(3,2) +1;
    end
end

confusion_matrix_kmeans

% figure
% scatter3(data_900(:,1),data_900(:,2),data_900(:,3), 100, true_900, 'filled')
% xlabel SL, ylabel SW, zlabel PL
% 
% figure
% scatter3(data_900(:,1),data_900(:,2),data_900(:,3), 100, clusters3, 'filled')
% xlabel SL, ylabel SW, zlabel PL


Cov = [1 0.5 0.3
       0.5 2 0
       0.3 0 3];
mu = [1 2 3]';

[U,L] = eig(Cov);
% L: eigenvalue diagonal matrix
% U: eigen vector matrix, each column is an eigenvector

% For N standard deviations spread of data, the radii of the eliipsoid will
% be given by N*SQRT(eigenvalues).

N = 1; % choose your own N
radii = N*sqrt(diag(L));

% generate data for "unrotated" ellipsoid
[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));

% rotate data with orientation matrix U and center mu
a = kron(U(:,1),xc);
b = kron(U(:,2),yc);
c = kron(U(:,3),zc);

data = a+b+c; n = size(data,2);

x = data(1:n,:)+mu(1);
y = data(n+1:2*n,:)+mu(2);
z = data(2*n+1:end,:)+mu(3);

% now plot the rotated ellipse
% sc = surf(x,y,z); shading interp; colormap copper
h = surfl(x, y, z); colormap copper
title('actual ellipsoid represented by mu and Cov')
axis equal
alpha(0.7)
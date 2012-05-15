function [ B, mean1, mean2 ] = init_values( A )
%UNTITLED Summary of this function goes here
%S.Sabev 2012

row_size = size(A(1,:));
column_size = size(A(:,1));

B = zeros(column_size(1), column_size(1));

for i=1:column_size
    for j=1:column_size
        B(i,j) = sqrt(sum((A(i,:)-A(j,:)).^2));
    end
end

[~, s_i] =max(B(:));
[mean1 mean2] = ind2sub(size(B), s_i);

mean1 = A(mean1,:);
mean2 = A(mean2,:);
end


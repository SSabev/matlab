function [new_mean] = NM(A,column_size, means)
%Given a dataset of points, this function finds an addition mean
%S.Sabev 2012
temp = 1;
NewMeans=zeros(column_size);

for i=1:column_size
    for j=1:size(means(:,1))
        temp = temp*sqrt(sum((A(i,:)-means(j,:)).^2));
    end
    NewMeans(i) = temp;
    temp = 1;
end

[~, s_i]=max(NewMeans(:));
[new_mean] = ind2sub(size(NewMeans), s_i);
new_mean = A(new_mean,:);

end
function [clusters, centres, sum_squared_func] = mykmeans(A, means)
%My own version of the well-known kmeans algorightm
%S.Sabev 2012

data_size = size(A(:,1));
means_size = size(means(:,1));
clusters = zeros(data_size(1),1);
temp_clusters = zeros(data_size(1),1);
centres = means;
iterations = 0;
data_row_size = size(means(:,1));
temp_distances = zeros(data_row_size(1), 1);
finished = 0;

sum_squared_func = 0;


while finished~=1
    iterations = iterations +1;
    for i=1:data_size(1)
        for j=1:means_size
           temp_distances(j) = sqrt(sum((A(i,:)-centres(j,:)).^2)); 
        end
        [~, s_i]=min(temp_distances);
        [min_value] = ind2sub(size(temp_distances), s_i);
        clusters(i) = min_value;
    end
    centres = zeros(means_size(1), 3);
    temp_cnt = zeros(means_size(1), 1);
    for i=1:data_size(1)
        centres(clusters(i), :) = centres(clusters(i),:) + A(i,:);
        temp_cnt(clusters(i)) = temp_cnt(clusters(i)) + 1;
    end
    for i=1:means_size(1)
        centres(i, :)=centres(i, :)./temp_cnt(i);
    end
    if temp_clusters==clusters
        finished = 1;
    else
        temp_clusters = clusters;
    end
end

for i = 1:data_size(1)
    sum_squared_func = sum_squared_func + (sum(sum(A(i, : ) - centres(clusters(i), : ))))^2;
end

sum_squared_func = sum_squared_func./data_size(1);

end
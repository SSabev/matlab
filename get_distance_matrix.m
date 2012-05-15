function [ distance ] = get_distance_matrix( input_matrix )
%get_distance_matrix Gets the distance matrix from a normal matrix
distance = zeros(90,90);
for i = 1:90
    for j = 1:90
        distance(i,j) = sqrt(sum((input_matrix(i,: ) - input_matrix(j,: )).^2));
    end
end
end
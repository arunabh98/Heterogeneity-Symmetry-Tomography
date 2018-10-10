function [clustered_projections, clustered_angles] = ...
    cluster_projections(projections, num_clusters, original_theta)
    % Cluster the projections.
    [idx, C, ~, ~] = kmeans(projections', num_clusters,...
        'distance', 'cityblock');
    clustered_angles = zeros(1, num_clusters);
    for i=1:num_clusters
        clustered_angles(i) = mean(original_theta(find(idx == i)));
    end
    clustered_projections = C';
end
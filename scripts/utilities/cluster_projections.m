function [clustered_projections, clustered_angles, cluster_class, original_cluster_class] = ...
    cluster_projections(projections, num_clusters, original_theta, original_class)
    % Cluster the projections.
    [idx, C, ~, ~] = kmeans(projections', num_clusters,...
        'distance', 'cityblock');
    clustered_angles = zeros(1, num_clusters);
    for i=1:num_clusters
        clustered_angles(i) = mean(original_theta(find(idx == i)));
    end

    original_cluster_class = zeros(1, num_clusters);
    for i=1:num_clusters
         original_cluster_class = mode(original_class(find(idx == i)));
    end
    clustered_projections = C';

    cluster_class = randi(3, 1, num_clusters);
end
function [clustered_projections, clustered_angles, cluster_class, original_cluster_class] = ...
    cluster_projections(projections, num_clusters, original_theta, original_class)
    % Cluster the projections.
    [idx, C, ~, ~] = kmeans(projections', num_clusters,...
        'distance', 'cityblock');

    % Make the Distance Matrix
    D = zeros(size(projections,2), num_clusters);
    for i=1:size(projections,2)
        distVec = repmat(projections(:,i), [1 num_clusters]) - C';
        D(i,:) = sum(distVec.^2);
    end

    % Filtering
    percentile95 = prctile(min(D,[],2), 80);
    filteredY = projections(:, min(D,[],2) < percentile95);
    filteredIdx = idx(min(D,[],2) < percentile95);
    filteredAngles = original_theta(min(D,[],2) < percentile95);
    filteredClass = original_class(min(D,[],2) < percentile95);

    fprintf(1,'Number of Projections Filtered: %d, out of Total: %d',...
        size(projections,2)-size(filteredY,2),size(projections,2));
    clusteredProj = zeros(size(projections,1), num_clusters);

    % Averaging
    count = 1;
    for ang = 1:num_clusters
        clusterProjs = filteredY(:, filteredIdx == ang);
        % If some clusters are obtained
        if(~(size(clusterProjs,2) == 0))
            clusteredProj(:,count) = mean(clusterProjs,2);
            count = count+1;
        end
    end

    % Trim out the non-existant angles
    clustered_projections = clusteredProj(:,1:count-1);

    % Calculate the angles and class of projections in each cluster.
    new_num_clusters = size(clustered_projections, 2);
    clustered_angles = zeros(1, new_num_clusters);
    original_cluster_class = zeros(1, new_num_clusters);
    for i=1:new_num_clusters
        clustered_angles(i) = mean(filteredAngles(filteredIdx == ang));
        original_cluster_class(i) = mode(filteredClass(filteredIdx == ang));
    end

    % Randomly initialize the clusters.
    cluster_class = randi(3, 1, new_num_clusters);

    % Analyze cluster purity.
    original_cluster_purity = zeros(1, num_clusters);
    filtered_cluster_purity = zeros(1, num_clusters);

    for ang = 1:num_clusters
        classProjs = original_class(idx == ang);
        filteredClassProjs = filteredClass(filteredIdx == ang);

        frequent_original_class = mode(classProjs);
        original_cluster_purity(ang) =...
            (sum(classProjs == frequent_original_class)/size(classProjs, 2))*100;

        if(size(filteredClassProjs, 2) == 0)
            filtered_cluster_purity(ang) = 0;
            break;
        end

        frequent_filtered_class = mode(filteredClassProjs);
        filtered_cluster_purity(ang) =...
            (sum(filteredClassProjs == frequent_filtered_class)/size(filteredClassProjs, 2))*100;
    end

    fprintf(1,'--------------------------------------------------------\n');
    fprintf(1,'Projections Clustered\n');
    fprintf(1,'--------------------------------------------------------\n');
end
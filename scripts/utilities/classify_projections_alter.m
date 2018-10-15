function [projection_classes] = ...
    classify_projections_alter(projections, original_theta, original_class,...
    sigmaNoise, output_size)

    % Agglomerative clustering.
    Z = linkage(projections', 'average', 'euclidean');
    num_clusters = inf;
    for c=1.14:0.001:1.16
        idx = cluster(Z, 'cutoff', c);
        % idx = cluster(Z, 'Maxclust', 3);
        
        num_clus = size(unique(idx), 1);
        
        if (num_clus < num_clusters) && (num_clus > 300)
            num_clusters = num_clus;
            cutoff = c;
            curr_idx = idx;
        end
    end
    
    C = zeros(size(projections, 1), num_clusters);
    class_clustered = zeros(1, num_clusters);
    theta_clustered = zeros(1, num_clusters);
    std_theta = zeros(1, num_clusters);
    std_class = zeros(1, num_clusters);
    idx = curr_idx;
    parfor i=1:num_clusters
        projections_in_cluster = projections(:, idx == i);
        C(:, i) = mean(projections_in_cluster, 2);
        class_clustered(i) = mode(original_class(idx == i));
        std_class(i) = std(original_class(idx == i));
        theta_clustered(i) = mean(original_theta(idx == i));
        std_theta(i) = std(original_theta(idx == i));
    end
    clustered_projections = C;

    % Denoise the projections.
%     sigmaNoise = sigmaNoise*num_clusters/size(original_theta, 2);
    clustered_projections = denoise(clustered_projections, sigmaNoise, 145, 100);
    clustered_projections = max(0, clustered_projections);

    clustered_projections = [clustered_projections flipud(clustered_projections)];

    % Constants
    num_angles = num_clusters;
    epsilon = 3e3;

    % Define the weight matrix.
    W = pdist2(clustered_projections', clustered_projections');
    W = W.^2;
    W = exp(-W/(2*epsilon));

    % Define diagonal matrix D.
    D = diag(sum(W, 2));
    W_tilda = D\W/D;
    D_tilda = diag(sum(W_tilda, 2));

    L = D_tilda - W_tilda;

    [V, ~] = eig(L, D_tilda);

    phi1 = -V(:, 2);
    phi2 = -V(:, 3);
    phi3 = -V(:, 4);
    coeff = [phi1 phi2 phi3];
    
    % Extract the three main eigenvectors.
    coeff = coeff(1:num_clusters, :);
    phi1 = phi1(1:num_clusters);
    phi2 = phi2(1:num_clusters);
    phi3 = phi3(1:num_clusters);
    
    % Visualize the clusters.
    h(1) = figure; scatter3(phi1, phi2, phi3, 10, class_clustered');
    
    % Cluster using single linkage hierachical clustering
    Z = linkage(coeff, 'single', 'euclidean');
    c = cluster(Z, 'Maxclust', 3);
    h(2) = figure; scatter3(phi1, phi2, phi3, 10, c);

    % Save the figures
    savefig(h, 'TwoFiguresFile.fig');
    close(h);
    projection_classes = c;
end
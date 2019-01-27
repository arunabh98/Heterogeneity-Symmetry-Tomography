function [projections, s_vec, original_class_of_projections] = ...
    get_projections(angles, testImStack, outlier_mode, outlier_percentage)

    % Initialize the classes and the projections.
    no_of_classes = size(testImStack, 3);
    if outlier_mode ~= 1
    	original_class_of_projections = randi(no_of_classes, 1, size(angles, 2));
    end
    	
    [sample_projection, s_vec] = radon(testImStack(:, :, 1), 0);
    projections = zeros(size(sample_projection, 1), size(angles, 2));

    for i=1:size(angles, 2)
        projections(:, i) = radon(testImStack(:,:, original_class_of_projections(i)), angles(i));
    end
end

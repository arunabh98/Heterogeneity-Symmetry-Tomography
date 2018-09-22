function [projections, s_vec, original_class_of_projections] = ...
    get_projections(angles, testIm1, testIm2, testIm3)

    % Initialize the classes and the projections.
    original_class_of_projections = randi(3, 1, size(angles, 2));
    [sample_projection, s_vec] = radon(testIm1, 0);
    projections = zeros(size(sample_projection, 1), size(angles, 2));

    for i=1:size(angles, 2)
        if original_class_of_projections(i) == 1
            projections(:, i) = radon(testIm1, angles(i));
        elseif original_class_of_projections(i) == 2
            projections(:, i) = radon(testIm2, angles(i));
        else
            projections(:, i) = radon(testIm3, angles(i));
        end 
    end
end

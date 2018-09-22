% Increase the number of parpool workers.
parpool('local', 14)
% warning('off', 'MATLAB:rankDeficientMatrix');

% Include the moment based estimation scripts and noise scripts.
addpath(genpath('../data'));
addpath(genpath('moment_based_estimation'));
addpath(genpath('noise_scripts'));
addpath(genpath('utilities'));
addpath(genpath('symmetry_gradient'));
addpath(genpath('single_axis_symmetry'));
addpath(genpath('horizontal_symmetry_utilities'));

% Get the images
image_size = 100;
P1 = read_process_image('refs_008.png', image_size);
P2 = read_process_image('refs_009.png', image_size);
P3 = read_process_image('refs_010.png', image_size);

% Constants.
non_uniform_distribution = 0;
sigmaNoiseFraction = 0.05;
if non_uniform_distribution == 0
    filename = ...
        strcat('../results/heterogeneity/', num2str(sigmaNoiseFraction*100), '_percent_noise/');
else
    filename = ...
        strcat('../results/heterogeneity/', num2str(sigmaNoiseFraction*100), '_percent_noise/non_uniform_distribution/');
end
output_size = max(size(P1));

% Experimentatal conditions.
max_shift_amplitude = 0;
symmetry_prior = 1;
noisy_orientations = 1;
symmetry_method = 4;
include_clustering = 0;
num_clusters = 180;
num_theta = 270;

% Create the folder to hold the results of the experiment.
mkdir(strcat(filename, num2str(num_theta), '/all_variables/'));
if include_clustering ~= 1
    theta_to_write = zeros(10, num_theta);
end

% The file which contains all the errors.
fileID = fopen(strcat(filename,...
    num2str(num_theta), '/result.txt'),'w');

% Write the original images.
imwrite(P1, strcat(filename,...
    num2str(num_theta), '/original_image_1.png'));
imwrite(P2, strcat(filename,...
    num2str(num_theta), '/original_image_2.png'));
imwrite(P3, strcat(filename,...
    num2str(num_theta), '/original_image_3.png'));

% Define ground truth angles and take the tomographic projection.
if non_uniform_distribution == 0
    theta = datasample(0:0.005:179, num_theta);
elseif non_uniform_distribution == 1
    theta = [datasample(0:0.005:20, num_theta/5)...
    datasample(40:0.005:60, num_theta/5)...
    datasample(80:0.005:120, 2*num_theta/5)...
    datasample(140:0.005:160, num_theta/5)];
end

% Generate the projections including outliers of class 1 and 2.
[projections, svector, original_class_of_projections] = ...
    get_projections(theta, P1, P2, P3);

% [projections, svector] = radon(P,theta);
original_projections = projections;
original_shifts = zeros(size(theta));

% Shift each projection by an unknown amount.
parfor i=1:size(projections, 2)
    original_shifts(i) = ...
        randi([-max_shift_amplitude, max_shift_amplitude]);
    projections(:, i) = circshift(projections(:, i), original_shifts(i)); 
end

if include_clustering ~= 1
    % Ignore shifts if clustering.
    theta_to_write(6, :) = original_shifts;
end

% Normalize s to a unit circle
smax = max(abs(svector));
svector = svector / smax;

% Define the length of the projection.
projection_length = size(projections, 1);

% Add noise to projections.
[measured_projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);

% Initial error between the projections.
disp('**** L2-norm error between the original projections and measured projections ****');
disp(norm(measured_projections - original_projections, 'fro'));
disp('');

estimated_class_of_projections = randi(3, 1, num_theta);
if include_clustering == 1
    disp('**** Initial - Cluster the projections ****');
    [clustered_projections, clustered_angles, cluster_class, original_class_of_projections] =...
        cluster_projections(measured_projections, num_clusters, theta, original_class_of_projections);

    % Save the original projections and mark the clustered projections as measured
    % projections.
    originally_measured_projections = measured_projections;
    measured_projections = clustered_projections;
    original_theta = theta;
    theta = clustered_angles;

    % Calculate the new variance of noise in the projections.
    sigmaNoise = sigmaNoise*num_clusters/num_theta;

    % Update the number of clusters after filtering.
    num_clusters = size(clustered_projections, 2);
    theta_to_write = zeros(10, num_clusters);

    % First estimate of the class of the projections.
    estimated_class_of_projections = cluster_class;
end

theta_to_write(1, :) = theta;
theta_to_write(2, :) = original_class_of_projections;

% If orientations are noisy or completely unknown.
if noisy_orientations == 1
    if include_clustering == 1 
        initial_theta = theta + randi([-1 1], 1, num_clusters);
    else 
        initial_theta = theta + randi([-1 1], 1, num_theta);
    end
else
    if include_clustering == 1
        initial_theta = randi([1 179], num_clusters, 1);
    else 
        initial_theta = randi([1 179], num_theta, 1);
    end
end

% The shifts estimated using the center of mass theorem.
estimated_shifts = estimate_shifts(measured_projections, max_shift_amplitude); 

if include_clustering == 1
    % Ignore shifts if we are clustering
    original_shifts = estimated_shifts;
end

% Error after shift correction.
disp('**** L1-norm error between estimated shifts and actual shifts ****');
disp(norm(estimated_shifts - original_shifts, 1));
disp('');
theta_to_write(7, :) = estimated_shifts;     

% Predict the angles and shifts using moment angle estimation.
disp('**** Moment based estimation ****');
[refined_projections, noisy_theta, projection_shifts, refined_classes] = ...
    SHARPord_cluster(measured_projections, svector, sigmaNoise, max_shift_amplitude,...
    estimated_shifts, initial_theta, noisy_orientations, estimated_class_of_projections);
disp('');
projection_shifts = projection_shifts';
noisy_theta = noisy_theta';

theta_to_write(3, :) = refined_classes;

% % Error after noise removal.
% % Only makes sense if clustering is not there.
% if include_clustering ~= 1
%     disp('**** L2-norm error between the original projections and denoised and shift corrected projections ****');
%     disp(norm(refined_projections - original_projections, 'fro'));
%     disp('');
% end
% disp('**** L1-norm error between the moment estimated shifts and the actual shifts ****');
% disp(norm(projection_shifts - original_shifts, 1));
% disp('');

% % Start iteration.
% projection_shifts = projection_shifts';
% noisy_theta = noisy_theta';

% better_theta = noisy_theta;
% better_shift = projection_shifts;

% % Define the gradient descent parameters.
% resolution_angle = 0.005;
% resolution_shift = 1;
% angle_amplitude = 3;
% shift_amplitude = 0;
% errors = [];
% alpha_rate =  0.001;
% beta_rate = 0.001;
% error_delta = 10;

% % Calculate the estimate of the image based on moment based solver.
% reconstructed_image = ...
%     reconstruct_image(refined_projections, better_theta,...
%         better_shift, output_size);

% % If symmetric prior is enabled, estimate the axis of symmetry.
% if symmetry_prior == 1
%     % delta is the angle by which the image should be rotated
%     % so that the axis of symmetry is horizontal.
%     delta = estimate_axis_symmetry_alter(...
%         reconstructed_image, output_size, 90, 0, symmetry_method);
% end

% % Record the results given by the moment based estimation.
% theta_to_write(2, :) = better_theta;
% theta_to_write(8, :) = -better_shift;
% formatSpec = 'Estimated image relative error: %4.2f \r\n';
% fprintf(fileID, formatSpec, calculate_rmse_error(reconstructed_image, P));
% fclose(fileID);


% % Store the estimated image after moment based estimation.
% imwrite(reconstructed_image, ...
%     strcat(filename, num2str(num_theta), '/estimated_image.png'));

% % Since we also have to run the optimization another time to show the
% % improvement by using symmetric priors we store the initial values.
% initial_image = reconstructed_image;
% initial_iteration_theta = better_theta;
% initial_iteration_shift = better_shift;

% if symmetry_prior == 1
%     disp('**** Optimization error using symmetric prior ****');
%     fprintf('\nIteration Error:            \n');
%     for i=1:500

%         % Use the projection prior
%         gradient_vector = ...
%             reconstructed_image -...
%             reconstruct_image(refined_projections, better_theta, better_shift, output_size);

%         % Then use the symmetric prior.
%         % For calculating the symmetry gradient we first rotate the image
%         % to make the horizintal axis symmetrical.
%         cropped_image = extract_circular_patch(reconstructed_image);

%         rotated_image = imrotate(cropped_image, delta, 'bicubic', 'crop');
%         symmetry_gradient_vector = symmetry_gradient_matrix(rotated_image);
%         symmetry_gradient_vector = imrotate(symmetry_gradient_vector, -delta, 'bicubic', 'crop');

%         % Now finally update the image.
%         reconstructed_image = ...
%             reconstructed_image - alpha_rate*gradient_vector - beta_rate*symmetry_gradient_vector; 

%         % The optimization error for this iteration.
%         function_error = calculate_optimization_error_symmetry_single_axis(refined_projections,...
%             reconstructed_image, better_theta, better_shift, delta);

%         % Display the error for ths iteration.
%         fprintf('%6g \n', function_error); 
%         errors = [errors function_error];

%         % Do a brute force over all orientations and shifts.
%         % Store the old estimated before calculating the new ones.
%         old_theta = better_theta;
%         old_shift = better_shift;
%         [better_theta, better_shift] = best_orientation_estimate(...
%             refined_projections, reconstructed_image, better_theta, better_shift,...
%             angle_amplitude, shift_amplitude, resolution_angle,...
%             resolution_shift);

%         % Now update the symmetry axis.
%         delta = estimate_axis_symmetry_alter(...
%             reconstructed_image, output_size, error_delta, delta, symmetry_method);

%         error_delta = max(1, error_delta - 0.5);
%     end
%     fprintf('\n');
%     disp('');
    
%     imwrite(reconstructed_image, ...
%         strcat(filename, num2str(num_theta), '/reconstructed_symmetry.png'));

%      % Plot the function error.
%     figure; plot(errors);
%     saveas(gcf, ...
%         strcat(filename, num2str(num_theta), '/error_symmetric.png'));

%     % Save the image formed due to symmetric priors.
%     symmetric_reconstructed_image = reconstructed_image;
% end

% % Now reinitialize the parameters for the second iteration without
% % using symmetry.
% reconstructed_image = initial_image;
% better_theta = initial_iteration_theta;
% better_shift = initial_iteration_shift;
% errors = [];
% alpha_rate = 0.001;

% disp('**** Optimization error without using symmetric prior ****');
% fprintf('\nIteration Error:            \n');
% for i=1:500

%     % Use the projection prior
%     gradient_vector = ...
%         reconstructed_image -...
%         reconstruct_image(refined_projections, better_theta, better_shift, output_size);

%     % Now finally update the image.
%     reconstructed_image = ...
%         reconstructed_image - alpha_rate*gradient_vector; 

%     % The optimization error for this iteration.
%     function_error = calculate_optimization_error(refined_projections,...
%         reconstructed_image, better_theta, better_shift);

%     % Display the error for ths iteration.
%     fprintf('%6g \n', function_error); 
%     errors = [errors function_error];

%     % Do a brute force over all orientations and shifts.
%     [better_theta, better_shift] = best_orientation_estimate(...
%         refined_projections, reconstructed_image, better_theta, better_shift,...
%         angle_amplitude, shift_amplitude, resolution_angle,...
%         resolution_shift);

%     % % Decrease the step size as the iteration progresses.
%     % if mod(i, 5) == 0
%     %     alpha_rate = alpha_rate/3;
%     % end
% end
% fprintf('\n');
% disp('');

% imwrite(reconstructed_image, ...
%     strcat(filename, num2str(num_theta), '/reconstructed_sans_symmetry.png'));

% % Reopen the file.
% fileID = fopen(strcat(filename,...
%     num2str(num_theta), '/result.txt'),'at');

% % Calculate the mse errors in both cases.
% if symmetry_prior == 1
%     formatSpec = 'Final symmetric image rmse error: %4.2f \r\n';
%     fprintf(fileID, formatSpec, calculate_rmse_error(symmetric_reconstructed_image, P));
% end
% formatSpec = 'Final non-symmetric image rmse error: %4.2f \r\n';
% fprintf(fileID, formatSpec, calculate_rmse_error(reconstructed_image, P));

% % Record the refined shifts and thetas.
% theta_to_write(3, :) = better_theta';
% theta_to_write(9, :) = -better_shift';  

% % Plot the function error.
% figure; plot(errors);
% saveas(gcf, ...
%     strcat(filename, num2str(num_theta), '/error_non_symmetry.png'))

% Write the thetas to csv file.
csvwrite(strcat(filename,...
    num2str(num_theta), '/thetas.csv'), theta_to_write);
fclose(fileID);

% Save the important variables.
save(strcat(filename, num2str(num_theta), '/all_variables/all_variables.mat'),...
    '-regexp',...
    '^(?!(clustered_projections|flipped_P1|flipped_P2|measured_projections|refined_projections|original_projections|originally_measured_projections|projections|first_quadrant|second_quadrant|third_quadrant|fourth_quadrant)$).');


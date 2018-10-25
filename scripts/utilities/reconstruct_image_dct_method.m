function image_estimate_dct = reconstruct_image_dct_method(...
	projections, theta_estimate, shift_estimate, output_size, num_freq)

	shifted_projections = correct_projection_shifts(projections, shift_estimate);
	image_estimate = iradon(shifted_projections, theta_estimate, output_size, 'Cosine');
	D = dctmtx(output_size);

	dct_matrix = zeros(size(projections, 1), num_freq^2);
    
    coeff_mat = zeros(output_size, output_size);
    for k=1:size(projections, 2)
        c = 1;
        for i=1:num_freq
            for j=1:num_freq
                basis_matrix = D(i,:)' * D(j,:);
                proj = radon(basis_matrix, theta_estimate(k));
                dct_matrix(:, c)= proj;
            	c = c + 1;
            end
        end

        coeff = pinv(dct_matrix)*projections(:, k);
        disp(mean(coeff));
        c = 1;
        for i=1:num_freq
            for j=1:num_freq
                coeff_mat(i, j) = coeff_mat(i, j) + coeff(c);
                c = c + 1;
            end
        end
    end
    coeff_mat = coeff_mat./size(projections, 2);

    image_estimate_dct = zeros(output_size, output_size);
	for i = 1:output_size
	    for j = 1:output_size
            basis_matrix = D(i,:)' * D(j,:);
	        image_estimate_dct = image_estimate_dct + coeff_mat(i, j)*basis_matrix;
	    end
	end
end

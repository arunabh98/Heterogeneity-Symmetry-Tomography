function image_estimate = reconstruct_image_dct_method(...
	projections, theta_estimate, shift_estimate, output_size, num_freq)

	shifted_projections = correct_projection_shifts(projections, shift_estimate);
	image_estimate = iradon(shifted_projections, theta_estimate, output_size, 'Cosine');

	coeff = zeros(output_size, output_size);
	for i=1:num_freq
		for j=1:num_freq
			basis_matrix = D(i,:)' * D(j,:);
			left_side = image_estimate*basis_matrix';
			right_side = basis_matrix*basis_matrix';
			coeff(i, j) = mean(left_side(:))/mean(right_side(:));
		end
	end

	image_estimate = zeros(output_size, output_size);
	for i = 1:image_size
	    for j = 1:image_size
	    	basis_matrix = D(i,:)' * D(j,:);
	        image_estimate = image_estimate + coeff(i, j)*basis_matrix;
	    end
	end
end
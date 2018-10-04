function proj = calc_prob_for_proj(projections, Im1, Im2, Im3, theta)
    proj = zeros(3, size(projections, 2));

    parfor i=1:size(projections, 2)
        estimated_projection1 = radon(Im1, theta(i));
        estimated_projection2 = radon(Im2, theta(i));
        estimated_projection3 = radon(Im3, theta(i));

        dist1 = norm(estimated_projection1 - projections(:, i)).^2;
        dist2 = norm(estimated_projection2 - projections(:, i)).^2;
        dist3 = norm(estimated_projection3 - projections(:, i)).^2;

        dist1 = exp(-dist1/2500);
        dist2 = exp(-dist2/2500);
        dist3 = exp(-dist3/2500);

        p = [(dist1/(dist1 + dist2 + dist3));
             (dist2/(dist1 + dist2 + dist3));
             (dist3/(dist1 + dist2 + dist3))];
         
        proj(:, i) = p;
    end
end
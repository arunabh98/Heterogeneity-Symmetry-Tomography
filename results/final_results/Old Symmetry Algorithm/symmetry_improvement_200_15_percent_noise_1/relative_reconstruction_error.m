a = im2double(imread('reconstructed_symmetry.png'));
b = im2double(imread('original_image.png'));
c = im2double(imread('reconstructed_sans_symmetry.png'));
error1 = calculate_mse_error(a, b);
error2 = calculate_mse_error(c, b);

re1 = 100*(error1/(norm(b)));
re2 = 100*(error2/(norm(b)));
% MATLAB code for image compression
%% Clear section 
clear 
clc
close all
%% Read the original image and get its size
image = imread('image1.bmp'); 
rows = size(image, 1); 
columns = size(image, 2);
%% Storing the corresponding color plane
% im2double is used to convert the image to double precision.
% Red component
R_component = im2double(image(:, :, 1)); 
% Green component
G_component = im2double(image(:, :, 2));
% Blue component
B_component = im2double(image(:, :, 3));
%% Displaying the color components
zero_matrix = zeros(size(R_component));
% cat is used to concatenate the red, green and blue components to 
form a RGB image.
% To display a color component use cat to concatenate the color matrix 
with zero matricies for the other colors.
redImage = cat(3, R_component, zero_matrix, zero_matrix);
greenImage = cat(3, zero_matrix, G_component, zero_matrix);
blueImage = cat(3, zero_matrix, zero_matrix, B_component);
subplot(2,2,1)
imshow(redImage)
title("Red Component")
subplot(2,2,2)
imshow(greenImage)
title("Green Component")
subplot(2,1,2)
imshow(blueImage)
title("Blue Component")
%% Process each color component in blocks of 8 × 8 pixels and get DCT
% dctmtx(n) returns the n-by-n DCT matrix, which you can use to 
perform a 2-D DCT on an image.
Temp = dctmtx(8);
% dct is used to get the 2-D DCT of a block
dct = @(block_struct) Temp * block_struct.data * Temp';
% B = blockproc(A,[m n],fun) processes the image A by applying the 
function fun to each distinct block of size [m n] and concatenating 
the results into the output image, B. 
% Red Component
Red_DCT_Blocks = blockproc(R_component, [8 8], dct);
% Green Component
Green_DCT_Blocks = blockproc(G_component, [8 8], dct);
% Blue Component
Blue_DCT_Blocks = blockproc(B_component, [8 8], dct);
%% Compressing and Decompressing 
% imfinfo returns a structure whose fields contain information about 
an image
info = imfinfo('image1.bmp');
% info.FileSize returns the size of the image in bytes, so to get it 
in MB divide by 2^20
image_size = info.FileSize / 2^20;
% fprintf is used to display the text in the command window
fprintf("The size of the original image is %f MB\n", image_size)
figure
PSNR = 1 : 4;
for m = 1 : 4
 % Masking
 % Compressor_Matrix is the matrix which I used to let the values
of a m*m block of each 8*8 block and change the other values to
0.So, its size is 8*8 and has ones in the first m*m block and
zeros otherwise.
 Compressor_Matrix = ones(m, m); 
 Compressor_Matrix(8-end+m, 8-end+m) = 0;
 % Compressor_Function is the function used to multiply each 8*8
block and the compressor matrix by multiplying corresponding
elements to clear the unwanted elements.
 Compressor_Function = @(block_struct) Compressor_Matrix.*
block_struct.data;
 Compressed_Red_DCT_Blocks = blockproc(Red_DCT_Blocks, [8 8],
Compressor_Function); 
 Compressed_Green_DCT_Blocks = blockproc(Green_DCT_Blocks, [8 8],
Compressor_Function); 
 Compressed_Blue_DCT_Blocks = blockproc(Blue_DCT_Blocks, [8 8],
Compressor_Function);
 % Compressing
 Compressed_Red_Component = zeros (rows/8*m, columns/8*m);
 Compressed_Green_Component = zeros (rows/8*m, columns/8*m);
 Compressed_Blue_Component = zeros (rows/8*m, columns/8*m);
 % This loop is used to decrease the size of the compressed matrix
by clearing the zero elements which are out of the first m*m
block of each 8*8 block
 k = 1;
 for i= 1:8:rows
 l = 1;
 for j = 1:8:columns
 Compressed_Red_Component(k:k+m-1, l:l+m-1) =
Compressed_Red_DCT_Blocks(i:i+m-1, j:j+m-1);
Compressed_Green_Component(k:k+m-1, l:l+m-1) = 
Compressed_Green_DCT_Blocks(i:i+m-1, j:j+m-1);
Compressed_Blue_Component(k:k+m-1, l:l+m-1) = 
Compressed_Blue_DCT_Blocks(i:i+m-1, j:j+m-1);
 l = l + m;
 end
 k = k + m;
 end
compressed_Image = cat(3, Compressed_Red_Component,
Compressed_Green_Component, Compressed_Blue_Component);
 % saving the compressed images and printing their sizes
 file_name = "compressed_" + m + ".bmp";
 imwrite(compressed_Image, file_name, 'bmp')
 info = imfinfo(file_name);
 image_size = info.FileSize/2^20;
 fprintf("The size of the compressed image at m = %d is %f MB\n",
m, image_size)
 % Decompressing
 % invdct is used to get the inverse DCT of a block
 invdct = @(block_struct) Temp' * block_struct.data * Temp;
 Decompressed_Red_Component = blockproc(Compressed_Red_DCT_Blocks ,
[8 8], invdct);
 Decompressed_Green_Component =
blockproc(Compressed_Green_DCT_Blocks, [8 8], invdct);
 Decompressed_Blue_Component =
blockproc(Compressed_Blue_DCT_Blocks, [8 8], invdct);
 Decompressed_Image = cat(3, Decompressed_Red_Component,
Decompressed_Green_Component, Decompressed_Blue_Component);
 % Displaying the decompressed images
 subplot(2,2,m)
 imshow(Decompressed_Image);
 file_name = "Decompressed_" + m + ".bmp";
 % imwrite is used to save the pictures
 imwrite(Decompressed_Image, file_name, 'bmp')
 title("Decompressed Image at m = " + m)
 % immse calculates the mean-squared error (MSE) between two
matrices. 
 mse = immse(im2uint8(image), im2uint8(Decompressed_Image));
 % Calculating the Peak signal-to-noise ratio (PSNR) for each value
of m
 PSNR(m) = 10 * (log10((255^2) / mse));
end
% Plotting PSNR against m
M = 1:4;
figure
plot(M, PSNR);
xlabel('m');
ylabel('PSNR');
title("Relation between PSNR and m");

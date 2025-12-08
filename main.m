% LPSD_v1.0 test demo
% Written by liujy

clc; clear;
close all;
warning('off');
addpath('utils', 'LPSD', genpath('data'), '-begin');

start_image = 1;
end_image = 36;
for image_idx = start_image:end_image
    fprintf("test image %d... ", image_idx);
    image_name1 = sprintf('pair%d_1.png', image_idx);
    image_name2 = sprintf('pair%d_2.png', image_idx);
    try
        test_img1 = im2uint8(imread(image_name1));
        test_img2 = im2uint8(imread(image_name2));
    catch ME
        fprintf('\nWarning: Failed to load images for test %d: %s\n', image_idx, ME.message);
        continue;
    end

    tic
    [H, ~, matched_points1, matched_points2, params] = LPSD_main(test_img1, test_img2);
    fprintf("%fs\n", toc);

    close all
    image_fusion(test_img1, test_img2, double(inv(H)));
    showMatchedFeatures(test_img1,test_img2,matched_points1,matched_points2,"montage");
end

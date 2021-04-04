%%

% mound_volume_calculation - Script to read host nuclei images, binarize
% them and calculate the resulting volume of the mounds using alpha shapes
%
% 1. Read the images of the nuclei, perform background and flatfield
% correction
% 2. Binarize the images using nested segmentation
% 3. Draw an alpha shape around each binarized slice and calculate area  
% 4. Calculate resulting volume 
% Last modified: Prathima Radhakrishnan : 2021-02-18


% all confocal stacks to be analyzed are of the form "1nemo_1", "1nemo_2", "1nemo_3", etc. - this step is a shorthand to read in these image files from the directory
serum_root = '1control_'; 

% read in the background image
dirname2= '/Volumes/homes/Prathima/2020_01_28_CNT_NEMO_RAPA_60x_opt1_02um/MED_60X_1X_background_3_MMStack_Pos0_1.ome.tif'; 

%convert image type to double
flat=double(imread(dirname2)); 

%read in the flatfield image
dirname3= '/Volumes/homes/Prathima/2020_01_28_CNT_NEMO_RAPA_60x_opt1_02um/MED_Concatenated Stacks.tif'; 

%convert image type to double
flat1=double(imread(dirname3)); 

% creating a uniform matrix of similar intensity to the flatfield image
% 65355 is the max intensity for 16bit image
flat2= double(ones(1024,1024)).*65355; 


%Calling forth the first confocal stack of images for this condition 
for k = 1; 
    %The following code applies to the first slice of the image
      for j=1
        mydata{k} = double(imread([serum_root num2str(k) '.tif'],j));%read in the first slice of this image
        C{k} = mydata{k} - flat;%background correction
        D{k} = C{k}./flat1;%flatfield correction
        E{k} = uint16(D{k}*65355);%This conversion corrects for the fact that the flatfield image, which is of a concentrated dye slide, 
        %is much brighter than the nucleus image. The number 65355 is used
        %here because it is the maximum value for a 16 bit image. 
        %However, any appropriate value can be used to make this
        %conversion.
        F = mat2gray(E{k});% Converting the matrix to a greyscale image
        F(F<=0.15*max(max(F))) = 0;%this is a form of background subtraction, in which all pixels below a certain intensity are converted to zero. 
        %I use this instead of the traditional background subtraction because later on, when I apply a mask to delete certain parts of the field of view, 
        %they turn completely black, which is often darker than the background of the original image and causes segmentation problems.
        %If necessary, vary the number here to modulate the amount of background subtracted.
        level = graythresh(F);%setting a threshhold to allow for the binarization of the image. 
        %A manual rather than automatic threshhold can be applied to
        %account for the fact that nuclei become dimmer father away from
        %the coverslip. Altering the amount of background subtracted is an
        %alternative method to accomplish the same task.
        bw = im2bw(F, level);%converting to a binary image
        bw2 = imfill(bw,'holes');%filling holes
        bw3 = imclose(bw2, strel('disk',4));%closing step to perfect segmentation
        bw4 = imopen(bw3, strel('disk',2));%opening step to perfect segmentation
        [x,y] = find(bw4);%noting the co-ordinates of the nuclei in the binarized image
        l = boundary(x,y);%drawing a boundary around the nuclei
        G = size(bw4);%noting the size of the image
        H{k,j} = poly2mask(y(l),x(l),G(1),G(2));%drawing a polygon around the image
        J{k} = bw4;
        imwrite(J{k},strcat('tile',num2str(k),'.TIFF'), 'WriteMode', 'append',  'Compression','none');%binarized image stack is formed
      end
    for j=2:7%the following code applies from the 2nd to the 7th slice of the image
        mydata{k} = double(imread([serum_root num2str(k) '.tif'],j));%read in 29 slices of this image
        C{k} = mydata{k} - flat;%background correction
        D{k} = C{k}./flat1;%flatfield correction
        E{k} = uint16(D{k}*65355);%This conversion corrects for the fact that the flatfield image, which is of a concentrated dye slide, 
        %is much brighter than the nucleus image.
        F = mat2gray(E{k});% Converting the matrix to a greyscale image
        F(F<=0.15*max(max(F))) = 0;%this is a form of background subtraction, in which all pixels below a certain intensity are converted to zero. 
        %I use this instead of the traditional background subtraction because later on, when I apply a mask to delete certain parts of the field of view, 
        %they turn completely black, which is often darker than the background of the original image and causes segmentation problems.
        %If necessary, vary the number here to modulate the amount of background subtracted.
        I{k} = H{k,j-1}.*F;%By multiplying with the mask of the polygon formed by the boundary of the previous z-slice, 
        %the empty area not occupied by cells in the previous z-slice is
        %not considered while analyzing the current z-slice. 
        level = graythresh(I{k});%setting a threshhold to allow for the binarization of the image. 
        %A manual rather than automatic threshhold can be applied to
        %account for the fact that nuclei become dimmer father away from
        %the coverslip. Altering the amount of background subtracted is an
        %alternative method to accomplish the same task.
        bw = im2bw(F, level);%converting to a binary image
        bw2 = imfill(bw,'holes');%filling holes
        bw3 = imclose(bw2, strel('disk',4));%closing step to perfect segmentation
        bw4 = imopen(bw3, strel('disk',2));%opening step to perfect segmentation
        [x,y] = find(bw4);%noting the co-ordinates of the nuclei in the binarized image
        l = boundary(x,y);%drawing a boundary around the nuclei
        G = size(bw4);%noting the size of the image
        H{k,j} = poly2mask(y(l),x(l),G(1),G(2));%drawing a polygon around the image
        J{k} = bw4;
        imwrite(J{k},strcat('tile',num2str(k),'.TIFF'), 'WriteMode', 'append',  'Compression','none');%binarized image stack is formed
    end
      for j=8:52%the following code applies from the 8th to the 52nd slice of the image
        mydata{k} = double(imread([serum_root num2str(k) '.tif'],j));%read in the first 19 slices of this image
        C{k} = mydata{k} - flat;%background correction
        D{k} = C{k}./flat2;%rather than applying a flatfield correction at these higher values of z, where there are fewer nuclei at the periphery, I am dividing by a uniform matrix of similar intensity to the flatfield image.
        E{k} = uint16(D{k}*65355);%This conversion corrects for the fact that the flatfield image, which is of a concentrated dye slide, is much brighter than the nucleus image
        F = mat2gray(E{k});% Converting the matrix to a greyscale image
        F(F<=0.10*max(max(F))) = 0;%this is a form of background subtraction, in which all pixels below a certain intensity are converted to zero. 
        %I use this instead of the traditional background subtraction because later on, when I apply a mask to delete certain parts of the field of view, 
        %they turn completely black, which is often darker than the background of the original image and causes segmentation problems.
        %If necessary, vary the number here to modulate the amount of background subtracted.
        I{k} = H{k,j-1}.*F;%By multiplying with the mask of the polygon formed by the boundary of the previous z-slice, 
        %the empty area not occupied by cells in the previous z-slice is
        %not considered while analyzing the current z-slice. 
        level = graythresh(I{k});%setting a threshhold to allow for the binarization of the image. 
        %A manual rather than automatic threshhold can be applied to
        %account for the fact that nuclei become dimmer father away from
        %the coverslip. Altering the amount of background subtracted is an
        %alternative method to accomplish the same task.
        bw = im2bw(F, level);%converting to a binary image
        bw2 = imfill(bw,'holes');%filling holes
        bw3 = imclose(bw2, strel('disk',4));%closing step to perfect segmentation
        bw4 = imopen(bw3, strel('disk',2));%opening step to perfect segmentation
        [x,y] = find(bw4);%noting the co-ordinates of the nuclei in the binarized image
        l = boundary(x,y);%drawing a boundary around the nuclei
        G = size(bw4);%noting the size of the image
        H{k,j} = poly2mask(y(l),x(l),G(1),G(2));%drawing a polygon around the image
        J{k} = bw4;
        imwrite(J{k},strcat('tile',num2str(k),'.TIFF'), 'WriteMode', 'append',  'Compression','none');%binarized image stack is formed
      end
      
      for j=53:57%the following code applies from the 53rd to the 57th slice of the image
        mydata{k} = double(imread([serum_root num2str(k) '.tif'],j));%read in the first 19 slices of this image
        C{k} = mydata{k} - flat;%background correction
        D{k} = C{k}./flat2;%rather than applying a flatfield correction at these higher values of z, where there are fewer nuclei at the periphery, I am dividing by a uniform matrix of similar intensity to the flatfield image.
        E{k} = uint16(D{k}*65355);%This conversion corrects for the fact that the flatfield image, which is of a concentrated dye slide, is much brighter than the nucleus image
        F = mat2gray(E{k});% Converting the matrix to a greyscale image
        F(F<=0.15*max(max(F))) = 0;%this is a form of background subtraction, in which all pixels below a certain intensity are converted to zero. 
        %I use this instead of the traditional background subtraction because later on, when I apply a mask to delete certain parts of the field of view, 
        %they turn completely black, which is often darker than the background of the original image and causes segmentation problems.
        %If necessary, vary the number here to modulate the amount of background subtracted.
        I{k} = H{k,j-1}.*F;%By multiplying with the mask of the polygon formed by the boundary of the previous z-slice, 
        %the empty area not occupied by cells in the previous z-slice is
        %not considered while analyzing the current z-slice. 
        level = graythresh(I{k});%setting a threshhold to allow for the binarization of the image. 
        %A manual rather than automatic threshhold can be applied to
        %account for the fact that nuclei become dimmer father away from
        %the coverslip. Altering the amount of background subtracted is an
        %alternative method to accomplish the same task.
        bw = im2bw(F, level);%converting to a binary image
        bw2 = imfill(bw,'holes');%filling holes
        bw3 = imclose(bw2, strel('disk',4));%closing step to perfect segmentation
        bw4 = imopen(bw3, strel('disk',2));%opening step to perfect segmentation
        [x,y] = find(bw4);%noting the co-ordinates of the nuclei in the binarized image
        l = boundary(x,y);%drawing a boundary around the nuclei
        G = size(bw4);%noting the size of the image
        H{k,j} = poly2mask(y(l),x(l),G(1),G(2));%drawing a polygon around the image
        J{k} = bw4;
        imwrite(J{k},strcat('tile',num2str(k),'.TIFF'), 'WriteMode', 'append',  'Compression','none');%binarized image stack is formed
      end
      
      for j=58:105%the following code applies from the 58th to the 105th slice of the image
        mydata{k} = double(imread([serum_root num2str(k) '.tif'],j));%read in the first 19 slices of this image
        C{k} = mydata{k} - flat;%background correction
        D{k} = C{k}./flat2;%rather than applying a flatfield correction at these higher values of z, where there are fewer nuclei at the periphery, I am dividing by a uniform matrix of similar intensity to the flatfield image.
        E{k} = uint16(D{k}*65355);%This conversion corrects for the fact that the flatfield image, which is of a concentrated dye slide, is much brighter than the nucleus image
        F = mat2gray(E{k});% Converting the matrix to a greyscale image
        F(F<=0.2*max(max(F))) = 0;%this is a form of background subtraction, in which all pixels below a certain intensity are converted to zero. 
        %I use this instead of the traditional background subtraction because later on, when I apply a mask to delete certain parts of the field of view, 
        %they turn completely black, which is often darker than the background of the original image and causes segmentation problems.
        %If necessary, vary the number here to modulate the amount of background subtracted.
        I{k} = H{k,j-1}.*F;%By multiplying with the mask of the polygon formed by the boundary of the previous z-slice, 
        %the empty area not occupied by cells in the previous z-slice is
        %not considered while analyzing the current z-slice. 
        level = graythresh(I{k});%setting a threshhold to allow for the binarization of the image. 
        %A manual rather than automatic threshhold can be applied to
        %account for the fact that nuclei become dimmer father away from
        %the coverslip. Altering the amount of background subtracted is an
        %alternative method to accomplish the same task.
        bw = im2bw(F, level);%converting to a binary image
        bw2 = imfill(bw,'holes');%filling holes
        bw3 = imclose(bw2, strel('disk',4));%closing step to perfect segmentation
        bw4 = imopen(bw3, strel('disk',2));%opening step to perfect segmentation
        [x,y] = find(bw4);%noting the co-ordinates of the nuclei in the binarized image
        l = boundary(x,y);%drawing a boundary around the nuclei
        G = size(bw4);%noting the size of the image
        H{k,j} = poly2mask(y(l),x(l),G(1),G(2));%drawing a polygon around the image
        J{k} = bw4;
        imwrite(J{k},strcat('tile',num2str(k),'.TIFF'), 'WriteMode', 'append',  'Compression','none');%binarized image stack is formed
      end
      
      for j=106:120%the following code applies from the 31st to the 50th slice of the image
        mydata{k} = double(imread([serum_root num2str(k) '.tif'],j));%read in the first 19 slices of this image
        C{k} = mydata{k} - flat;%background correction
        D{k} = C{k}./flat2;%rather than applying a flatfield correction at these higher values of z, where there are fewer nuclei at the periphery, I am dividing by a uniform matrix of similar intensity to the flatfield image.
        E{k} = uint16(D{k}*65355);%This conversion corrects for the fact that the flatfield image, which is of a concentrated dye slide, is much brighter than the nucleus image
        F = mat2gray(E{k});% Converting the matrix to a greyscale image
        F(F<=0.3*max(max(F))) = 0;%this is a form of background subtraction, in which all pixels below a certain intensity are converted to zero. 
        %I use this instead of the traditional background subtraction because later on, when I apply a mask to delete certain parts of the field of view, 
        %they turn completely black, which is often darker than the background of the original image and causes segmentation problems.
        %If necessary, vary the number here to modulate the amount of background subtracted.
        I{k} = H{k,j-1}.*F;%By multiplying with the mask of the polygon formed by the boundary of the previous z-slice, 
        %the empty area not occupied by cells in the previous z-slice is
        %not considered while analyzing the current z-slice. 
        level = graythresh(I{k});%setting a threshhold to allow for the binarization of the image. 
        %A manual rather than automatic threshhold can be applied to
        %account for the fact that nuclei become dimmer father away from
        %the coverslip. Altering the amount of background subtracted is an
        %alternative method to accomplish the same task.
        bw = im2bw(F, level);%converting to a binary image
        bw2 = imfill(bw,'holes');%filling holes
        bw3 = imclose(bw2, strel('disk',4));%closing step to perfect segmentation
        bw4 = imopen(bw3, strel('disk',2));%opening step to perfect segmentation
        [x,y] = find(bw4);%noting the co-ordinates of the nuclei in the binarized image
        l = boundary(x,y);%drawing a boundary around the nuclei
        G = size(bw4);%noting the size of the image
        H{k,j} = poly2mask(y(l),x(l),G(1),G(2));%drawing a polygon around the image
        J{k} = bw4;
        imwrite(J{k},strcat('tile',num2str(k),'.TIFF'), 'WriteMode', 'append',  'Compression','none');%binarized image stack is formed
      end
      
       for j=121:127%the following code applies from the 121st to the 127th slice of the image
        mydata{k} = double(imread([serum_root num2str(k) '.tif'],j));%read in the first 19 slices of this image
        C{k} = mydata{k} - flat;%background correction
        D{k} = C{k}./flat2;%rather than applying a flatfield correction at these higher values of z, where there are fewer nuclei at the periphery, I am dividing by a uniform matrix of similar intensity to the flatfield image.
        E{k} = uint16(D{k}*65355);%This conversion corrects for the fact that the flatfield image, which is of a concentrated dye slide, is much brighter than the nucleus image
        F = mat2gray(E{k});% Converting the matrix to a greyscale image
        F(F<=0.35*max(max(F))) = 0;%this is a form of background subtraction, in which all pixels below a certain intensity are converted to zero. 
        %I use this instead of the traditional background subtraction because later on, when I apply a mask to delete certain parts of the field of view, 
        %they turn completely black, which is often darker than the background of the original image and causes segmentation problems.
        %If necessary, vary the number here to modulate the amount of background subtracted.
        I{k} = H{k,j-1}.*F;%By multiplying with the mask of the polygon formed by the boundary of the previous z-slice, 
        %the empty area not occupied by cells in the previous z-slice is
        %not considered while analyzing the current z-slice. 
        level = graythresh(I{k});%setting a threshhold to allow for the binarization of the image. 
        %A manual rather than automatic threshhold can be applied to
        %account for the fact that nuclei become dimmer father away from
        %the coverslip. Altering the amount of background subtracted is an
        %alternative method to accomplish the same task.
        bw = im2bw(F, level);%converting to a binary image
        bw2 = imfill(bw,'holes');%filling holes
        bw3 = imclose(bw2, strel('disk',4));%closing step to perfect segmentation
        bw4 = imopen(bw3, strel('disk',2));%opening step to perfect segmentation
        [x,y] = find(bw4);%noting the co-ordinates of the nuclei in the binarized image
        l = boundary(x,y);%drawing a boundary around the nuclei
        G = size(bw4);%noting the size of the image
        H{k,j} = poly2mask(y(l),x(l),G(1),G(2));%drawing a polygon around the image
        J{k} = bw4;
        imwrite(J{k},strcat('tile',num2str(k),'.TIFF'), 'WriteMode', 'append',  'Compression','none');%binarized image stack is formed
      end
end


%%

dirname= '/Volumes/homes/Prathima/2019_10_30_siEcad_60x_opt1/tile1.TIFF'; %the stack of binarized nuclei


for j=1:50%number of slices
    nuclei{j}=double(imread(dirname,j)); %reading in the binarized stack and converting image type to double
    [x,y] = find(nuclei{j});
    shp = alphaShape(x,y,50);%forming an alphashape around the nuclei
    A = area(shp);%calculating the area of the alpha shape
    B{j} = A;
end

Volume = sum(cell2mat(B)).*0.7;% volume is the sum of all of the alpha shapes multiplied by the increment between z slices - make sure to multiply this value by the conversion to microns cubed
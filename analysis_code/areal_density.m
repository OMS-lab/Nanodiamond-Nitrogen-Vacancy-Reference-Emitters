%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used to calculated the areal density of nanodiamonds from a
% dark-field microscope image of the sample
% Simply "run" to calculate.
% note: If a different image is selected, then image filtering
% parameters may need adjustment to remove anomalous feature detections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load RGB microscope image
if ~exist("nds_rgb","var")
    nds_rgb = imread("D:\20230109_NPL_Nanodiamonds\microscope_images\middle_darkfield_singleImage_B.jpg");
end
nds_gray = rgb2gray(nds_rgb); %convert to grayscale
nds_gray = imcrop(nds_gray,[3.5 0.5 2461 2057]); %crop out metadata
res = 20/size(imcrop(nds_gray,[2094.5 1935.5 292 22]),2); %um/px (scale bar is 20 µm) %calculate resolution
nds_thresh = nds_gray > 50; %threshold grayscale image

%plot grayscale image
figure(1);clf
imagesc(nds_gray);

%feature detection via regionprops
stats = regionprops("table",nds_thresh,"Centroid", ...
        "MajorAxisLength","MinorAxisLength","BoundingBox","ConvexImage", ...
        "Area");
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;

f = diameters < 100 & diameters > 3; %filter out anomalously detected objects

%circle detected objects
hold on
viscircles(stats.Centroid(f,:),radii(f,:),'color','yellow','LineWidth',0.1);
hold off
colormap("bone")
axis image

%print the calculated density to console
disp(strcat("Density is ",num2str(nnz(f)/(res*size(nds_gray,1)*res*size(nds_gray,1))),"µm-2"))
function [pslgDomain] = findBoundaryPolygons(filename,variables)
%% Image editing
% Read file
domain_image=imread(filename);

% Variables
treshold_value = variables.treshold_value;
smallest_cluster_in_px = variables.smallest_cluster_in_px;
image_contains_outer_boundary = variables.image_contains_outer_boundary;
Lx = variables.Lx;
Ly = variables.Ly;
smooth_param = variables.smooth_param;
skip_param = variables.skip_param;

% Convert to grayscale and create binary image:
if size(domain_image,3)==3
    domain_image = rgb2gray(domain_image);
end
img_bw = domain_image>treshold_value;
img_bw = 1-img_bw;

% divide into clusters:
[img_clusters,number_of_clusters] = bwlabel(img_bw,4);

% Delete smallest clusters
stats = regionprops(img_clusters,'Area');
clusters_to_keep = find([stats.Area] > smallest_cluster_in_px);
img = ismember(img_clusters,clusters_to_keep);
% Remove regions inside and small structures
img = bwmorph(img,'erode');

img = 1-img;
img(round(end/2)+1:end,:) = imfill(img(round(end/2)+1:end,:),'holes');
img(1:round(end/2),:) = imfill(img(1:round(end/2),:),'holes');
without_border = imclearborder(img);
img_boundary = (img-without_border);
img = without_border;
img = imfill(img,'holes');

img = bwmorph(img,'clean');
img = bwmorph(img,'hbreak');
img = bwmorph(img,'spur');
%img = bwmorph(img,'majority');

% Delete smallest clusters again
img = bwlabel(img,4);
stats = regionprops(img,'Area');
clusters_to_keep = find([stats.Area] > smallest_cluster_in_px);
img = ismember(img,clusters_to_keep);

%% MESH
% Find holes:
B_all = bwboundaries(img,8);

img = bwlabel(img);
holes = regionprops(img,'pixelList');

% Extract boundaries and scale to domain (-1,1)x(-1,1)
x=[];
y=[];
pslgInner.holes=[];
pslgInner.segments=[];
pslgInner.segmentmarkers=[];

for i = 1:length(B_all)
    B = B_all{i};
    [~,ind] = unique(B,'rows');
    B = B(sort(ind),:);
    B = B(1:skip_param:end,:);
    B(:,1)=2*B(:,1)/size(img,1)-1;
    B(:,2)=2*B(:,2)/size(img,2)-1;
    B_per_smooth = [B(end-smooth_param:end,:);B(1:end,:);B(1:smooth_param,:)];
    x_polygon = smooth(B_per_smooth(:,2),smooth_param)'+1e-9*rand(1,length(B_per_smooth)); x_polygon = x_polygon(smooth_param+2:end-smooth_param);
    y_polygon = -smooth(B_per_smooth(:,1),smooth_param)'+1e-9*rand(1,length(B_per_smooth)); y_polygon = y_polygon(smooth_param+2:end-smooth_param);
    x = [x x_polygon];
    y = [y y_polygon];
    size(holes)
    size(B_all)
    
    xrand = 2*holes(i).PixelList(:,1)'/size(img,2)-1;
    yrand = -2*holes(i).PixelList(:,2)'/size(img,1)+1;
    for j = 1:length(xrand)
        if inpolygon(xrand(j),yrand(j),[x_polygon x_polygon(1)],[y_polygon y_polygon(1)])
            pslgInner.holes = [pslgInner.holes [xrand(j);yrand(j)]];
            break;
        end
    end
    pslgInner.segments = [pslgInner.segments length(pslgInner.segments)+[ (1:length(x_polygon));[(2:length(x_polygon)),1]]];
    pslgInner.segmentmarkers = [pslgInner.segmentmarkers ones(1,length(x_polygon))];
end

pslgInner.points = [x;y];

if image_contains_outer_boundary
    %Extract outer walls:
    B = bwboundaries(img_boundary,4,'holes');
    [cell_length] = cellfun(@length, B);
    [~,max_ind] = max(cell_length);
    B = B{max_ind};
    [~,ind] = unique(B,'rows');
    B = B(sort(ind),:);
    B = B(1:skip_param:end,:);
    B(:,1)=2*B(:,1)/size(img,1)-1;
    B(:,2)=2*B(:,2)/size(img,2)-1;
    x_polygon = smooth(B(:,2),smooth_param)'+1e-9*rand(1,length(B));
    y_polygon = -smooth(B(:,1),smooth_param)'+1e-9*rand(1,length(B));
    
    %Delete edge points for sharp boundary:
    ind_right = find(x_polygon>0.9);
    x_polygon(ind_right)=[];
    y_polygon(ind_right) =[];
    x_polygon = [x_polygon(1:(min(ind_right)-1)) 1 1 x_polygon((min(ind_right)):end)];
    y_polygon = [y_polygon(1:(min(ind_right)-1)) 1 -1 y_polygon((min(ind_right)):end)];
    ind_left = find(x_polygon<-0.9);
    x_polygon(ind_left)=[];
    y_polygon(ind_left) =[];
    x_polygon = [x_polygon(1:(min(ind_left)-1)) -1 -1 x_polygon((min(ind_left)):end)];
    y_polygon = [y_polygon(1:(min(ind_left)-1)) -1 1 y_polygon((min(ind_left)):end)];
    ind_left = find(x_polygon<-(1-1e-9));
    ind_right = find(x_polygon>(1-1e-9));
    pslgOuter.points = [x_polygon(:)';y_polygon(:)'];
    tmp = [1:size(pslgOuter.points,2)-1; 2:size(pslgOuter.points,2)];
    tmp(end)=1;
    pslgOuter.segments = size(pslgInner.points,2) + tmp;
    tmp_markers = ones(1,size(pslgOuter.segments,2));
    tmp_markers(ind_left(1))=2; tmp_markers(ind_right(1))=3;
    if size(pslgInner.points,2)>0
        pslgOuter.segmentmarkers = max(pslgInner.segmentmarkers)+tmp_markers;
    else
        pslgOuter.segmentmarkers = tmp_markers;
    end
else
    x_wall_all{1} = [-1 1];
    x_wall_all{2} = 0;
    x_wall_all{3} = -1;
    x_wall_all{4} = -1;
    pslgOuter.points = [[x_wall_all{3}       1+0*x_wall_all{1}       fliplr(x_wall_all{4})   -1+0*x_wall_all{2}];
        [1+0*x_wall_all{3}   fliplr(x_wall_all{1})   -1+0*x_wall_all{4}      x_wall_all{2}]];
    pslgOuter = fliplr(pslgOuter);
    tmp = [1:size(pslgOuter.points,2)-1; 2:size(pslgOuter.points,2)];
    tmp(end)=1;
    pslgOuter.segments = size(pslgInner.points,2) + tmp;
    pslgOuter.segmentmarkers = max(pslgInner.segmentmarkers)+ones(1,size(pslgOuter.segments,2));
end

% Combine:
if size(pslgInner.points,2)>0
    pslgDomain.points = [pslgInner.points pslgOuter.points];
    pslgDomain.holes = pslgInner.holes;
    pslgDomain.segments = [pslgInner.segments pslgOuter.segments];
    pslgDomain.segmentmarkers = [pslgInner.segmentmarkers pslgOuter.segmentmarkers];
else
    pslgDomain = pslgOuter;
end

%Scale to domain size
pslgDomain.points(1,:) = pslgDomain.points(1,:)*Lx;
pslgDomain.points(2,:) = pslgDomain.points(2,:)*Ly;
if size(pslgInner.points,2)>0
    pslgDomain.holes(1,:) = pslgDomain.holes(1,:)*Lx;
    pslgDomain.holes(2,:) = pslgDomain.holes(2,:)*Ly;
else
   pslgDomain.holes=[]; 
end


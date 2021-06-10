%Example script to create vectorized plot from an image.
%


%%  Set Variables
var.treshold_value = 50; % threshold value for creating binary image from color image. Transforms image to grayscale and then applies the criterion. (1-255)
var.smallest_cluster_in_px = 30; %Smallest allowed cluster in pixels that is inluded in the vectorized plot
var.image_contains_outer_boundary = false; % Extract outer walls if they exist. If set to false, a box with dimensions [-Lx,Lx],[-Ly,Ly] is returned
var.Lx = 2; % X-axis is scaled to [-Lx,Lx]
var.Ly = 1; % Y-axis is scaled to [-Ly,Ly]
var.smooth_param = 5; % Apply smoothing to the boundaries to avoid pixelized look (smooth(x,smooth_param) is used)
var.skip_param = 2; % Skip points to avoid pixelized look.

%% Extract points and segments
pslg = findBoundaryPolygons('fracture_network.png',var);

%% plot the segments
x = [pslg.points(1,pslg.segments(1,:));pslg.points(1,pslg.segments(2,:))];
y = [pslg.points(2,pslg.segments(1,:));pslg.points(2,pslg.segments(2,:))];
plot(x,y,'-k')
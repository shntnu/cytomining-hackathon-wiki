function SPADE_contour2D(x,y,contour_levels,outlier_percentage)
% FlowJo_contour2D(x,y, contour_levels, outlier_percentage)
% outlier_percentage = 1, 2, 5, 10, ...
% example:
%       x = randn(1,50000) + [zeros(1,20000), ones(1,30000)]; 
%       y = randn(1,50000)+ [3*ones(1,20000),zeros(1,30000)];
%       outlier_percentage = 5;
%       FlowJo_contour2D(x,y, [0:15:60,70:10:90,93:3:99])

if exist('contour_levels')~=1
    contour_levels = [0:15:60,70:10:90,93:3:99];
end

if exist('outlier_percentage')~=1
    outlier_percentage=5;
end

if length(x)>200000
    ind = randsample(1:length(x),200000);
    x_downsample = x(ind); 
    y_downsample = y(ind);
else
    x_downsample = x;      
    y_downsample = y;
end
display('computing ksdensity grids for each variable ...')
[F,XI,hx] = ksdensity(x_downsample); 
[F,YI,hy] = ksdensity(y_downsample); 
Z = zeros(length(XI),length(YI));

display('computing the gaussian kernel density of each point on the grid ...')
dist_XI_x = exp(-(XI(:)*ones(1,length(x_downsample)) - ones(length(XI),1)*(x_downsample(:)')).^2/(2*hx^2));
dist_YI_y = exp(-(YI(:)*ones(1,length(y_downsample)) - ones(length(YI),1)*(y_downsample(:)')).^2/(2*hy^2));
Z = (dist_XI_x*dist_YI_y')';
Z = Z./(sqrt(2*pi*hx^2*2*pi*hy)*length(x_downsample));


Y = sort(Z(:));
[dummy,ind] = min(abs(cumsum(Y)-sum(Y)*outlier_percentage/100));
outer_ring = Y(ind);

display('drawing the contour plot ...')
contour(XI,YI,Z,prctile(Y(ind:end),contour_levels));
axis_tmp = axis;

% display('drawing outliers ...')
% 
% delta_XI = XI(2)-XI(1);
% delta_YI = YI(2)-YI(1);
% flag_dense_points = zeros(1,length(x));
% for i=1:size(Z,1)
%     start_cell=[];end_cell=[];
%     for j=1:size(Z,2)
%         if Z(i,j)<outer_ring && isempty(start_cell)
%             continue;
%         end
%         if Z(i,j)>=outer_ring && isempty(start_cell)
%             start_cell = j; end_cell=j;
%         end
%         if Z(i,j)>=outer_ring && ~isempty(start_cell)
%             end_cell=j;
%         end
%         if Z(i,j)<outer_ring && ~isempty(start_cell)
%             flag_dense_points(x>=XI(start_cell)-delta_XI/2 & x<=XI(end_cell)+delta_XI/2 & abs(y-YI(i))<=delta_YI/2)=1;
%             start_cell=[]; end_cell=[];
%         end
%     end
% end
% hold on; plot(x(flag_dense_points==0),y(flag_dense_points==0),'.','color',[0.3 0.3 0.3],'markersize',4); hold off
% 

return
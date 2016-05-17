function draw_colorbar(position,colors,xticklabels)

x = position(1); y = position(2);
len = position(3); height = position(4);
% patch([x, x+len, x+len, x],[y, y, y+height, y+height],[1 1 1]); hold on; 
delta_x = len/size(colors,1);
for i=1:size(colors,1)
    color_tmp = colors(i,:);
    patch([i-1,i,i,i-1]*delta_x+x,[0,0,height,height]+y,color_tmp,'edgecolor',color_tmp); 
end
line([x, x+len, x+len, x, x],[y, y, y+height, y+height, y],'color',[0 0 0])
text(x-10,       0+y-7, xticklabels{1},'fontsize',10)
text(x+len/2, 0+y-7, xticklabels{2},'fontsize',10) 
text(x-10+len,   0+y-7, xticklabels{3},'fontsize',10) 

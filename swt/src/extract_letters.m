f =  fspecial('gaussian', 2*ceil(1*2.5)+1, 1);
% g_ = imfilter(g,f);
g_ = g;
[a,b] = SWT_c(g_, 1, 500);


THRESHOLD_X = 8;
THRESHOLD_Y = 8;
%%
no_of_components = max(max(b));
bounding_boxes = cell(no_of_components,1);
component_regions = cell(no_of_components,1);
is_text = ones(no_of_components,1);
%%
for i = 1:no_of_components
% for i = 1:5
    component_regions{i} = find(b == i);
    [cr_x,cr_y] = ind2sub(size(g), component_regions{i});
    patch = a(cr_x, cr_y);
    [min_x, min_x_posn] = min(cr_x);
    [max_x, max_x_posn] = max(cr_x);
    [min_y, min_y_posn] = min(cr_y);
    [max_y, max_y_posn] = max(cr_y);
    
    if(max_x - min_x < THRESHOLD_X || max_y - min_y < THRESHOLD_Y)
        is_text(i) = 0;
    end
    
    angle1 = abs(atan2(cr_y(max_x_posn) - cr_y(min_x_posn) , max_x - min_x));
    angle2 = abs(atan2(cr_x(max_y_posn) - cr_x(min_y_posn) , max_y - min_y));
    if(~(angle1 > 10 * 180 /pi || angle1 < 80 * 180 /pi || ...
        angle2 > 10 * 180 /pi || angle2 < 80 * 180 /pi))
        is_text(i) = 0;
    end

    bounding_boxes{i} = g(min_x:max_x, min_y:max_y);
end

%%
figure
j = 1;
k = 75;
% bbr = bounding_boxes{is_text};
for i=k:k+25
    subplot(5,5, j);
    j = j+1;
    if(is_text(i) == 1)
        imshow(bounding_boxes{i});
    end
    hold on
end
%%
imwrite(imcomplement(bounding_boxes{2}), 'hopeful.jpg', 'jpg')

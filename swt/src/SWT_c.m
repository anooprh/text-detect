function [swt swtcc] = SWT_c( img, dol, maxWidth )

if size( img, 3 ) == 3
    img = rgb2gray(img);
end
img = im2single(img);

edgeMap = single( edge( img, 'canny', .15 ) ); 
img = imfilter( img, fspecial('gauss',[5 5], 0.3*(2.5-1)+.8) );
gx = imfilter( img, fspecial('prewitt')' ); %//'
gy = imfilter( img, fspecial('prewitt') );
gx = single(medfilt2( gx, [3 3] ));
gy = single(medfilt2( gy, [3 3] ));

[swt swtcc] = Source( edgeMap.', gx.', gy.', dol, maxWidth ); %//'

swt = swt'; %//'
swtcc = double(swtcc'); %//'
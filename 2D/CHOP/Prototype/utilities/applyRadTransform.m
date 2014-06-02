%====================================================================
%> @brief This function applies radial transform to a 16-bit image,
%> effectively reducing it to a 8-bit image.
%>
%> @author anonymous
%>
%> De�i�iklik                                  Tarih       Yapan
%> Version 1.0                          02/01/2011 Anonymous
%> Yorumlar eklendi                 11/08/2011 Umit Rusen Aktas
%>
%> @param im16: 16-bit image.
%>
%> @retval imEnh: The enhanced 8-bit image.
%>
%> @example
%> img = imread('523232_liman1.tif');
%> imEnh = applyRadTransform(img);
%> figure,imshow(imEnh(:,:,1:3));
%====================================================================
function [imEnh] = applyRadTransform(im16)
    %% Read band number.
    [~, ~, k] = size(im16);
    %% Allocate space
    im = zeros(size(im16));
    %% Apply transform to each band.
    for i=1:k
        maxVal = double(max(max(im16(:,:,i))));
        im(:,:,i) = round(double(im16(:,:,i))./maxVal*255);
    end
    im = uint8(im);
    imEnh = uint8(im);
    for i=1:k
        imEnh(:,:,i) = imadjust(im(:,:,i),stretchlim(im(:,:,i),[0.05 0.95]));
    end
end
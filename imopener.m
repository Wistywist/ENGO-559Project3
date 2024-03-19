function im1 = imopener(impath)
im1 = fopen(impath);
im1 = fread(im1,[960,720]);
im1 = uint8(im1);

end

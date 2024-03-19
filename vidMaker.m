function v = vidMaker(initial_vid,frames,vidName)

v = VideoWriter(vidName);
open(v)

for frame = 1:size(frames,3)
    writeVideo(v,frames(:,:,frame));
end

close(v)
end


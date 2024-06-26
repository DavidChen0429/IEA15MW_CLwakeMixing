imageDir = 'Data/Figures/LiDAR/1D_circle_Uni8_300s';
imageFiles = dir(fullfile(imageDir, 'figure_*.png'));

% Extract the indices from the filenames
indices = zeros(1, length(imageFiles));
for i = 1:length(imageFiles)
    indices(i) = str2double(regexp(imageFiles(i).name, '\d+', 'match'));
end

% Sort the indices and corresponding filenames
[~, sortedIndices] = sort(indices);
sortedFiles = imageFiles(sortedIndices);

%% Create a VideoWriter object
Video = VideoWriter(fullfile(imageDir, 'LiDAR_sampling.avi'));
Video.FrameRate = 24;
open(Video);

for i = 1:length(sortedFiles)
    img = imread(fullfile(imageDir, sortedFiles(i).name));
    writeVideo(Video, img);
end
close(Video);

%% Test
LiDARAvi = VideoReader(fullfile(imageDir,"LiDAR_sampling.avi"));
i = 1;
while hasFrame(LiDARAvi)
   mov(i) = im2frame(readFrame(LiDARAvi));
   i = i + 1;
end
vf = figure(Position=[0 0 LiDARAvi.Width LiDARAvi.Height]);
imshow(mov(1).cdata,Border="tight")
movie(vf,mov,1,LiDARAvi.FrameRate)
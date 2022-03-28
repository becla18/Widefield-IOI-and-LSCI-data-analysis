function signal = extractROISignal(file)
%EXTRACTROISIGNAL Code designed to extract differential movement signal from avi webcam video
v = VideoReader(file);
signal = zeros(v.NumFrames,1);
frame = readFrame(v);
ROI = ROIselect(frame); %Binary mask
npixelROI = sum(sum(ROI)); %no. of pixels in ROI
for i = 1:v.NumFrames
    frame = read(v,i);
    frame = double(rgb2gray(frame)); %Convert to proper gray scale format
    if i == 1
        lastFrame = frame;
    end
    diffFrame = abs(frame-lastFrame);
    signal(i) = sum(sum(diffFrame.*ROI))/npixelROI; %mean diff signal in ROI
    lastFrame = frame;
end
end


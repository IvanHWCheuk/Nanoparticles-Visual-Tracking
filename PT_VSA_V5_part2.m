% Hoi Wa Cheuk (Ivan)
% Video Processing Research for NUANCE TEM Lab
% Code for Part 2

particleTraced = 5;

video = VideoReader('Part 2 Video.avi');
frame = 1:video.NumFrames;
second = frame./30;
areas = zeros(1,length(frame));
pos = zeros(length(frame),2);
increment = 80; % in nm of the video's scope

% Analyse the first frame
vidFrame = readFrame(video); 
Ffill = processImage(vidFrame);
[stat,Bound] = locateObjects(Ffill,1500);

% Find the scope of the video, assuming the box is longest
scopeStat = regionprops(Ffill,'MajorAxisLength');
scope = 0;
for i = 1:length(scopeStat)
    if scopeStat(i).MajorAxisLength > scope
        scope = scopeStat(i).MajorAxisLength;
    end
end

% Assigning lock for one specific polymer tracing
lock = stat(particleTraced).Centroid;
radius = sqrt(stat(particleTraced).Area / pi);

for ind = 1:video.NumFrames-1
    vidFrame = readFrame(video);
    
    % Slow motion gap
    pause(1/video.FrameRate);
    
    % Process image and locate objects
    Ffill = processImage(vidFrame);
    [stat,Bound] = locateObjects(Ffill,1500);
    
    for j = 1:length(stat)
        if sum(abs(stat(j).Centroid-lock)) < radius + radius*0.4
            newstat = stat(j);
            newBound = Bound(j);
        end
    end
    lock = newstat.Centroid;
    radius = sqrt(newstat.Area / pi);
    
    % Display each frame with label and outline on objects
    imshow(vidFrame); hold on
    for k = 1 : length(newBound) % Don't need loop here since only one
        b = newBound{k};
        c = newstat(k).Centroid;
        plot(b(:,2),b(:,1),'g','linewidth',2);
        text(c(1),c(2),num2str(k),'backgroundcolor','g');
    end 
    hold off
    
    % Record data
    areas(ind) = newstat.Area;
    pos(ind,:) = newstat.Centroid;
    
end

totvel = zeros(1,length(frame));
for i = 1:length(pos)-1
    totvel(i) = sqrt((pos(i+1,1)-pos(i,1))^2 + (pos(i+1,2)-pos(i,2))^2);
end

% Adjust in nm and seconds
areasNM = areas .* (increment/scope)^2;
posNM = pos .* (increment/scope);

totvelNM = zeros(1,length(frame));
for i = 1:length(posNM)-1
    totvelNM(i) = sqrt((posNM(i+1,1)-posNM(i,1))^2 + (posNM(i+1,2)-posNM(i,2))^2);
end

figure(2)
subplot(1,3,1); plot(second(1:end-2),areas(1:end-2)); 
xlabel('Time (s)'); ylabel('Area (pixels)');
subplot(1,3,2); plot(pos(1:end-2,1),pos(1:end-2,2)); 
xlabel('x-position (pixels)'); ylabel('y-position (pixels)');
subplot(1,3,3); plot(second(1:end-2),totvel(1:end-2)); 
xlabel('Time (s)'); ylabel('Total Velocity (pixels/frame)');




function Ffin = processImage(image)
    % Adjust image to black and white given automated gray threshold
    level = graythresh(image);
    Fgray = rgb2gray(image);
    Fbw = imbinarize(Fgray,level-0.07);
    
    % Fill holes within object
    Ffin = imfill(Fbw,'holes');
end

function [outStat, outBound] = locateObjects(image, areaThreshold)
    outBound = bwboundaries(image);
    outStat = regionprops(image,'Centroid','Area');
    
    % Clean the image from noise by deleting small areas
    for i = length(outStat):-1:1
        if outStat(i).Area < areaThreshold
            outStat(i) = [];
            outBound(i,:) = [];
        end
    end
end

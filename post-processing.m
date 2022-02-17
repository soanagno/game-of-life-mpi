% --- GAME OF LIFE POST-PROCESSING --- %
clc; close all; clear all;

% Just change the number of cores used for compliling 

cores = 18;  % Number of processors / cores

color = 4;
head = 2;  % Lines to skip in the text files which contain processor info

% Initialisations
A = zeros(cores, 6);  % Create 2D matrix A to store info about each core
for i = 1:cores
    fileID = fopen(string(i - 1) + '_pixels.txt', 'r');
    A(i, :) = fscanf(fileID, '%d', [1 6]);
    fclose(fileID);
end
m_max = A(1, 4); n_max = A(1, 5);  % Size [m, n] of largest core
cm_max = A(1, 2); cn_max = A(1, 3);  % Size [cm, cn] of full grid
cores = A(1, 1); iter = A(1, 6);
data = zeros((iter + 1) * m_max, n_max, cores);
corner_x = zeros(1, cores + 1); corner_y = zeros(1, cores + 1);
C = cell(1, cores);

% Extract data from each processor (core) file
for i = 1:cores
    m = A(i, 4); n = A(i, 5);
    data(1:m*iter, 1:n, i) = dlmread(string(i - 1) + '_pixels.txt', ' ', [head 0 m*iter+1 n - 1]);
    i;
end

% Extract corner coordinates from each core
for i = 1:cores
    corner_x(1, i) = dlmread(string(i - 1) + '_pixels.txt', ' ', [1 0 1 0]);
    corner_y(1, i) = dlmread(string(i - 1) + '_pixels.txt', ' ', [1 1 1 1]);
end

% Keep unique corner points for grid plot
corner_x = unique(corner_x) + 0.5;
corner_y = unique(corner_y) + 0.5;

% Initialise video
myVideo = VideoWriter('GameofLife'); % open video file
myVideo.FrameRate = 20; % can be adjusted, 5 - 10 works well
open(myVideo)

for i = 1:iter
    
    for j = 1:cores
        % Get the data for each core
        cn = A(j, 2); cm = A(j, 3); m = A(j, 4); n  = A(j, 5);
        temp = color * data((i-1) * m + 1:m * i, 1:n, j);
        % Merge them into one cell array C
        C(1, j) = {temp'};
    end
    
    % Create the full grid b reshaping C in the initial dimensions cm, cn
    full_grid = cell2mat(reshape(C', [cm, cn]));
    
    % PLOTS
    image(full_grid')  % plot image
    pbaspect([n/m 1 1])
    title('Game of Life')
    set(gca, 'XAxisLocation','top');  % Place x axis on top of the plot
    colormap(bone(5))  % Set color theme
    
    % Set custom grid based on individual processor area using the
    % coordinates of their corners, [corner_x, corner_y].
    grid on;
    set(gca,'xtick', corner_y,'TickLength',[0 0], 'GridColor', 'red', 'LineWidth', 2, 'GridAlpha', 0.6);
    set(gca,'ytick', corner_x,'TickLength',[0 0], 'GridColor', 'red', 'LineWidth', 2, 'GridAlpha', 0.6);
    xticklabels({0:cn_max}); yticklabels({0:cm_max});
    
    frame = getframe(gcf);  % get frame
    % writeVideo(myVideo, frame);  % save frame to video
    pause(0.01);  % Pause for animation plot
end
close(myVideo)
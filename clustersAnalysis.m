%% Temporal
k = 10;
idx = kmeans(HbOMatrix,k);
movie = HbOMovie.data;
mean_frames = zeros(size(movie,1),size(movie,2),k);
for i = 1:k
    time_points = find(idx == i);
    mean_frames(:,:,i) = mean(movie(:,:,time_points),3);   
end

figure;
t = tiledlayout(2,5);
for i = 1:k
    nexttile;
    imagesc(mean_frames(:,:,i),[-5 10]);
    title(['Cluster' num2str(i)]);
end
%% Spatial
k = 10;
idx = kmeans(HbOMatrix',k);
tmp_mask = zeros(size(mask));
tmp_mask(logical(mask(:))) = idx;
figure;
imagesc(tmp_mask);
colorbar;




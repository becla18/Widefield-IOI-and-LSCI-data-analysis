function showMovie(movie,colorLimits)
%showMovie displays the 3D matrix movie as a sequence of colormapped images using imagesc. Can
%optionnally specify range for the colormap by entering a two-element vector colorLimits =
%[lowerLimit,upperLimit]
nframes = size(movie,3);
figure
ax = gca;
colormap hot
imagesc(movie(:,:,1));
c = colorbar;
if nargin == 1
    colorLimit1 = c.Limits(1);
    colorLimit2 = c.Limits(2);
else
    colorLimit1 = colorLimits(1);
    colorLimit2 = colorLimits(2);
end
ax.Visible = 'off';
for i = 2:nframes
    imagesc(movie(:,:,i),[colorLimit1 colorLimit2]);
    colorbar;
    ax.Visible = 'off';
    pause(0.0001);
end

function [HbOMovie,HbRMovie] = convertToHb(GreenMovie,RedMovie)
%convertToHb converts green and red relative reflectance data to changes in Hbr and HbO concentration.
%Data must be in 3D array format (2D images by time)

%Parameters
whichCurve = 'Ma';
rescaling_factor = 1e6;
lambda1=450;
lambda2=700;
npoints=1000;
baseline_hbt = 100;
baseline_hbo = 60;
baseline_hbr = 40;
filter = 0; %change to 1 if fluorescence imaging is performed during IOI

eps_pathlength = ioi_epsilon_pathlength(lambda1,lambda2,npoints,whichCurve,baseline_hbt,baseline_hbo,baseline_hbr,filter);
Ainv=single(rescaling_factor*pinv(eps_pathlength)); %Pseudoinverse matrix

LogGreenAvg = -log(reshape(GreenMovie,1,[]));
LogRedAvg = -log(reshape(RedMovie,1,[]));

LogMat = [LogRedAvg;LogGreenAvg];
Hbs = Ainv*LogMat;
HbOMovie = Hbs(1,:);
HbRMovie = Hbs(2,:);

HbOMovie = reshape(HbOMovie, size(GreenMovie));
HbRMovie = reshape(HbRMovie, size(GreenMovie));

%Protection against aberrant data points
HbOMovie(isnan(HbOMovie)) = 0;
HbOMovie(isinf(HbOMovie)) = 0;
HbRMovie(isnan(HbRMovie)) = 0;
HbRMovie(isinf(HbRMovie)) = 0;

end


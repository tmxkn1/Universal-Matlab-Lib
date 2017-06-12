indices = 1:1000;
data = zeros(size(indices));

% some regions of data
region1 = [50:100 200:340 450:500 670:980];
region2 = setdiff(indices, region1);

% generating random data
data(region1) = rand(size(region1)) + 1;
data(region2) = rand(size(region2));


figure(1);
cla(gca);
hold on;
plot(region1, data(region1));
plot(region2, data(region2));
hold off;
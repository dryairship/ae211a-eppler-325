data       = load('Data.txt');
X1          = data(:,1);
Y1          = data(:,2);
X=X1/1000;
Y=Y1/1000;
subplot(1,1,1);
p = plot(X,Y,'b');
p.LineWidth=5;
grid;
axis([-0.01 0.11 -0.030 0.030]);
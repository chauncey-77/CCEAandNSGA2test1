dat = load('outdata.txt');
y1 = dat(:,1)';
y2 = dat(:,2)';
xlabel('f1'),ylabel('f2');
plot(y1,y2,'.');
axis([0 1 0 2]);
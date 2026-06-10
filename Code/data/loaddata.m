clear;

currdir=pwd;

[DATA,TXT] = xlsread('data.xlsx','matlab');
 
Time = datenum(TXT(2:end,1),'dd/mm/yyyy');

Dates = cellstr(datestr(Time));

VarNames = TXT(1,2:end);  

Y = DATA;

y = [100*log(Y(:,1:8)), Y(:,9:12), 100*log(Y(:,13:15)), Y(:,16), 100*log(Y(:,17)), Y(:,18)];
 
save data y Time VarNames Dates Y;
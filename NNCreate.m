clear; clc;

load('reference.mat');

% net=newff(minmax(P),[points_number,3,3],{'logsig' 'logsig' 'purelin'},'trainbfg');
% net=newff(P,T,[21 14]);
net = feedforwardnet([50 100]);

net.layers{1}.transferFcn = 'radbas'; % 'radbas'; % 'tansig'; % poslin
net.layers{2}.transferFcn = 'radbas'; % 'poslin'; % 'tansig'; % 'radbas'; % 
net.layers{3}.transferFcn = 'purelin'; % 'tansig'; % 'radbas'; % 
net.trainFcn = 'trainlm'; % 'trainbfg'; % 
net.performFcn= 'mse'; % 'sse'; % 
net.trainParam.goal=0.001;
net.trainParam.epochs=2000;
net.sampleTime = Tq;

[net,tr]=train(net,P,T');

% gensim(net)

save nn3 net;
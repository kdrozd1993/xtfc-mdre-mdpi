% Aquiring sensitivity results

clear all
clc

load('R:\xtfc-mrde\Trials\Fifth_Order_Fluid_Catalytic_Reactor\Results_MDRE\Sensitivity_Tanh_n160m120_20211222.mat');

% I = 1;
% N = 8;
% cnt = 0;

% for i = I:N
%     
%     cnt = cnt+1;
%     
%     mat(1,cnt) = mean(output.results(i).compTime.inv);
%     mat(2,cnt) = std(output.results(i).compTime.inv);
%     mat(3,cnt) = mean(output.results(i).compTime.iter);
%     mat(4,cnt) = std(output.results(i).compTime.iter);
%     mat(5,cnt) = sqrt(sum(output.results(i).lossNorm.interp.^2));
%     mat(6,cnt) = sqrt(sum(output.results(i).residualNorm.PinnKE.^2));
%     mat(7,cnt) = sqrt(sum(output.results(i).residualNorm.PinnLyap.^2));
%     mat(8,cnt) = sqrt(sum(output.results(i).residualNorm.LyapKE.^2));
%     
% end

for i = 1:length(output.results)
    
    mat(1,i) = mean(output.results(i).compTime.inv);
    mat(2,i) = std(output.results(i).compTime.inv);
    mat(3,i) = mean(output.results(i).compTime.iter);
    mat(4,i) = std(output.results(i).compTime.iter);
    mat(5,i) = sqrt(sum(output.results(i).lossNorm.interp.^2));
    mat(6,i) = sqrt(sum(output.results(i).residualNorm.PinnKE.^2));
    mat(7,i) = sqrt(sum(output.results(i).residualNorm.PinnLyap.^2));
    mat(8,i) = sqrt(sum(output.results(i).residualNorm.LyapKE.^2));
    
end

a = mat(:,1:3:length(output.results));
b = mat(:,2:3:length(output.results));
c = mat(:,3:3:length(output.results));
% d = mat(:,4:6:length(output.results));
% e = mat(:,5:6:length(output.results));
% f = mat(:,6:6:length(output.results));
% g = mat(:,7:12:length(output.results));
% h = mat(:,8:12:length(output.results));
% i = mat(:,9:12:length(output.results));
% j = mat(:,10:12:length(output.results));
% k = mat(:,11:12:length(output.results));
% l = mat(:,12:12:length(output.results));

mmat = [a,b,c];



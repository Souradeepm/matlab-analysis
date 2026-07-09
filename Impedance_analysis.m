clc;
clear all;
close all;
repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
[data text raw]=xlsread(input_path);
text=char(text(1,:));
[p1,q]=size(text);
temperature=zeros(p1,1);
for p=1:p1
temp=text(p,:);
temp(temp=='K')=[];
temp=str2double(temp);
temperature(p)=temp;
end;
temperature=temperature(~isnan(temperature));
drt_input_analysis_matlab2011(repo_root,input_path,temperature,data);
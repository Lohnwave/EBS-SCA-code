
clear
clc;

mex  cec13_func.cpp -DWINDOWS 
D=30;
M=10;
Xmin=-100;
Xmax=100;
pop_size=50;
iter_max=6000;
max_fes=300000;
warning off;
runs=30;
func_number=28;

algorithm_Name={'EBS_SCA'};
fhd=str2func('cec13_func');
Data1=zeros(func_number,runs,size(algorithm_Name,2));   % Store Error
Data2=zeros(func_number,runs,size(algorithm_Name,2));   % Store the Fitness value of the Global optimum solution
Data3=zeros(runs,iter_max,size(algorithm_Name,2));      % Store the fitness value of the best individual in the search process
Data4=zeros(runs,D,func_number,size(algorithm_Name,2)); % Store the Global optimum solution
t=zeros(func_number,runs,size(algorithm_Name,2));
f_mean=zeros(func_number,size(algorithm_Name,2));
t_mean=zeros(func_number,size(algorithm_Name,2));
E_mean=zeros(func_number,size(algorithm_Name,2));
for k=1:size(algorithm_Name,2)
    fprintf('Algorithm =\t %d\n',k);
    for i=1:func_number
        fprintf('Problem =\t %d\n',i);
        for j=1:runs
            fprintf('run =\t %d\n',j);
            BSOFUNC=algorithm_Name{k};
            tic;
            [Error,gbest,gbestval,allgbestval,fitcount]=feval(BSOFUNC,fhd,...
                iter_max,max_fes,pop_size,D,M,Xmin,Xmax,i);
            t(i,j,k)=toc;
            Data1(i,j,k)=Error;
            Data2(i,j,k)=gbestval;
            Data3(j,:,i,k)=allgbestval;
            Data4(j,:,i,k)=gbest;
        end
        f_mean(i,k)=mean(Data2(i,:,k),2);
        t_mean(i,k)=mean(t(i,:,k),2);
        E_mean(i,k)=mean(Data1(i,:,k),2);
    end
    file_name= [algorithm_Name{k},'_M10_30D_30runs.mat'];
    save (file_name);
end











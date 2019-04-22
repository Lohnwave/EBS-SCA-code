%% Information
% Last Revised: 26th July 2018

% Enhanced Brain Storm SCA (EBS-SCA) (Standard Version)

% Developer: CHUNQUAN LI, ZU LUO, ZHENSHOU SONG, FENG YANG, JINGHUI FAN, and PETER X. LIU. 

% Note that in order to obey the copy-right rules, please cite our published paper properly.

% INPUTS:

% max_iteration:                Maximum number of iterations
% max_FES:                      Maximum number of fitness evaluations
% N:                            Population size
% D:                            Dimension
% M:                            Number of clusters
% lb:                           Lower bound of a problem
% ub:                           Upper bound of a problem

% OUTPUTS:

% Error:                        Error between the fitness of the Global optimum solution and target value
% X_gbest:                      Global optimum solution 
% gbestval:                     Fitness value of the Global optimum solution
% allgbestval:                  Store the fitness value of the best individual in the search process
% fitcount:                     Number of fitness evaluations
function[Error,X_gbest,gbestval,allgbestval,fitcount] = EBS_SCA(fhd,max_iteration,max_FES,N,D,M,lb,ub,varargin) 
		%Load Target value
        load fbias_data;
        % *************** Initialization *******************
        if size(lb,2)==1
            X = lb + (ub - lb) * rand(N,D); % initialize the population of individuals        
            popu_sorted  = lb + (ub - lb) * rand(N,D); % initialize the  population of individuals sorted according to clusters        
        else
            X = repmat(lb,N,1) + repmat((ub - lb),N,1).* rand(N,D);
            popu_sorted  = repmat(lb,N,1) + repmat((ub - lb),N,1).* rand(N,D);
        end
        n_iteration = 1; % current iteration number
        % initialize cluster probability to be zeros
        prob = zeros(M,1); % n_c: number of clusters        
        best = zeros(M,1); % index of best individual in each cluster
        if size(lb,2)==1
            centers = lb + (ub - lb) * rand(M,D);  % initialize best individual in each cluster
            centers_copy = lb + (ub - lb) * rand(M,D);  % initialize best individual-COPY in each cluster FOR the purpose of introduce random best
        else
            centers = repmat(lb,M,1) + repmat((ub - lb),M,1).* rand(M,D);  % initialize best individual in each cluster
            centers_copy = repmat(lb,M,1) + repmat((ub - lb),M,1).* rand(M,D);
        end
        best_fitness = 1000000*ones(max_iteration,1);%initialize the best fitness in every iteration
        fitness_X = 1000000*ones(N,1);  % store fitness value for each individual
        fitness_X_sorted = 1000000*ones(N,1);  % store  fitness value for each sorted individual        
        M_rang_l=ones(1,D).*lb;
        M_rang_r=ones(1,D).*ub;       
        indi_temp = zeros(1,D);  % store temperary individual        
        fitcount=0;
        for i=1:N
            fitness_X(i,:) = (feval(fhd,X(i,:)',varargin{:}))';
        end
        fitcount=fitcount+N;
        [gbestval,gbestid]=min(fitness_X);
        best_fitness(n_iteration,1)=gbestval;
        X_gbest=X(gbestid,:);
        pbestpos=X;
        pbestval=fitness_X;
        gbestpos_o=X_gbest(1,:);
        gbestval_o=best_fitness(1,1);
        % Parameter seting
        pr11=0.7;
        pr12=0.7;
        pr1_low=0.2;
        pr1_high=0.7;
        pr1=pr1_low+pr1_high*(1:max_iteration)*(1/max_iteration);%0.2-0.9
		w=0.8;		% the constant scale factor for distinguishing whether to use IIU-I or IIU-II
        a=1;        % corresponding to the parameter ¦Ñ in the article
        sg=5;       % corresponding to the parameter ¦« in the article
        hnum=N/M;   % corresponding to the parameter ¦Ä in the article
        
        while n_iteration<=max_iteration && fitcount<=max_FES                
            n_iteration=n_iteration+1;
            st=0;
            if n_iteration <= w.*max_iteration
                stop=0;
                st=1;
            end            
            if stop >= sg || st==1;%stop >=5 || n_iteration <= Gap_ratio.*max_iteration
            % **************** IUS-I & IUS-II ******************
            % if st==1 i.e., n_iteration <= Gap_ratio.*max_iteration,
            % execute IUS-I, otherwise execute IUS-II.
                stop=0;
                % ***************** RG ******************
                %Step 2:Clustering
                a1=randperm(N);
                for i=1:M
                    cluster(a1((i-1)*hnum+1:i*hnum),:)=i;
                end
				% initialized the fitness for the center of each cluster and assigning a initial big fitness value as best fitness for each cluster in minimization problems				
                fit_values = inf*ones(M,1);
				% initialize 0 individual in each cluster				
                number_in_cluster = zeros(M,1);  
                % find the best individual in each clusters
                for idx = 1:N  
					% number_in_cluster is a matrix initializedwith five rows and one column, which is used to described the population size for five clusters.
                    number_in_cluster(cluster(idx,1),1)= number_in_cluster(cluster(idx,1),1) + 1; 
					% minimization
                    if fit_values(cluster(idx,1),1) > pbestval(idx,1)  
                        fit_values(cluster(idx,1),1) = pbestval(idx,1);
						% store the index of the best individual in each clusters that is the center as each cluster.
                        best(cluster(idx,1),1) = idx;
                    end
                end
                % form population sorted according to clusters
				% initialize cluster counter to be 0
                counter_cluster = zeros(M,1);  
				% initialize accumulated number of individuals in previous clusters
                acculate_num_cluster = zeros(M,1);  
            
                for idx =2:M
                    acculate_num_cluster(idx,1) = acculate_num_cluster((idx-1),1) + number_in_cluster((idx-1),1);%
                end
				% idx is the index of each individual after sorting
                for idx = 1:N 
					%number_in_cluster is a matrix initializedwith five rows and one column
                    counter_cluster(cluster(idx,1),1) = counter_cluster(cluster(idx,1),1) + 1 ;
					% temIdx is the index of each individual after sorting					
                    temIdx = acculate_num_cluster(cluster(idx,1),1) +  counter_cluster(cluster(idx,1),1); 
					% change the individual into new position in popu_sorted varian
                    popu_sorted(temIdx,:) = pbestpos(idx,:); 
					% change the fitness of individual into the new position
                    fitness_X_sorted(temIdx,1) = pbestval(idx,1); 
                end
                %record the best individual in each cluster as the center of each cluster
                for idx = 1:M
                    centers(idx,:) = pbestpos(best(idx,1),:);
                end
                centers_copy = centers;  % make a copy   
                
                %************** MIGM *****************
                for idx = 1:N 
                    r_1 = rand();
                    if r_1 < pr1(n_iteration) % select one cluster
                        r_11 =rand();
                        if   r_11 < pr11
                            % A new idea is generated by dimensionally
                            % randomly selecting a center.
                            % Generate a new idea
                            for dim=1:D
                                % Randomly select a center
                                cluster_1 = ceil(rand() * M);   
                                indi_temp(1,dim)= centers(cluster_1,dim);
                            end
                        else
                            % A new idea is generated by randomly selecting
                            % two ideas from one group.
                            indi_1 = acculate_num_cluster(cluster(idx),1) + ceil(rand()*number_in_cluster(cluster(idx),1));
                            indi_2 = acculate_num_cluster(cluster(idx),1) + ceil(rand()*number_in_cluster(cluster(idx),1));
                            r=rand(1,D);
                            % Generate a new idea
                            indi_temp(1,:) = r.*popu_sorted(indi_1,:)+(1-r).*popu_sorted(indi_2,:);
                        end
                    else
                        r_12 = rand();
                        if r_12 < pr12 
                            % A new idea is generated by randomly selecting
                            % two ideas from two different groups.
                            tem=rand(1,D);
                            cluster_1 = ceil(rand() * M);
                            cluster_2 = ceil(rand() * M);
                            indi_1 = acculate_num_cluster(cluster_1,1) + ceil(rand() * number_in_cluster(cluster_1,1));
                            indi_2 = acculate_num_cluster(cluster_2,1) + ceil(rand() * number_in_cluster(cluster_2,1));
                            % Generate a new idea
                            indi_temp(1,:) = tem.* popu_sorted(indi_1,:) + (1-tem).* popu_sorted(indi_2,:);
                        else   
                            % A new idea is generated by dimensionally randomly
                            % selecting one idea.
                            % Generate a new idea
                            for dim=1:D
                                cluster_1 = ceil(rand() * M);
                                indi_1 = acculate_num_cluster(cluster_1,1) + ceil(rand() * number_in_cluster(cluster_1,1));
                                indi_temp(1,dim)= popu_sorted(indi_1,dim) ;
                            end
                        end
                    end
                    
                    %************* Updating individuals ***************
                    r2=(2*pi)*rand(1,D);                    % corresponding to the parameter ¦Ë2 in the article
                    r3=2*rand(1,D);                         % corresponding to the parameter ¦Ë3 in the article
                    r1=a-n_iteration*((a)/max_iteration);   % corresponding to the parameter ¦Ë1 in the article
                    
                    if n_iteration <= w.*max_iteration
                    %*************** IUS-I: Eq. (14)************************
                        if rand()<=0.5 %Eq. (14)
                            indi_temp(1,:) = indi_temp(1,:)+(r1.*sin(r2).*abs(r3.*(indi_temp(1,:)-pbestpos(idx,:))));
                        else
                            indi_temp(1,:) = indi_temp(1,:)+(r1.*cos(r2).*abs(r3.*(indi_temp(1,:)-pbestpos(idx,:))));
                        end
                    else
                    %*************** IUS-II: Eq. (1)************************
                       if rand()<=0.5  
                            indi_temp(1,:) = indi_temp(1,:)+(r1.*sin(r2).*abs(r3.*X_gbest(1,:)-indi_temp(1,:)));
                       else
                            indi_temp(1,:) = indi_temp(1,:)+(r1.*sin(r2).*abs(r3.*X_gbest(1,:)-indi_temp(1,:)));
                        end
                    end
					% Check boundary
                    indi_temp=(indi_temp<M_rang_l).*M_rang_l+(indi_temp>M_rang_r).*M_rang_r...
                        +((indi_temp>=M_rang_l)&(indi_temp<=M_rang_r)).*indi_temp;
					% Evaluate fitness value
                    fv = feval(fhd,indi_temp(1,:)',varargin{:});
                    fitcount=fitcount+1;
					% better than the previous one, replace.
                    if fv < pbestval(idx,1)  
                        pbestval(idx,1)=fv;
                        pbestpos(idx,:)=indi_temp(1,:);
                        X(idx,:) = indi_temp(1,:);
                    end
                end
            else
                % ****************** IUS-II ************************
                for idx=1:N
                    r2=(2*pi)*rand(1,D);
                    r3=2*rand(1,D);
                    r1=a-n_iteration*((a)/max_iteration);
                    %*************** IUS-II: Eq. (1)************************
                    if rand()<=0.5 
					   X(idx,:) = X(idx,:)+(r1.*sin(r2).*abs(r3.*X_gbest(1,:)-X(idx,:)));
                    else
					   X(idx,:) = X(idx,:)+(r1.*sin(r2).*abs(r3.*X_gbest(1,:)-X(idx,:)));
                    end
					% Check boundary
                    X(idx,:)=(X(idx,:)<M_rang_l).*M_rang_l+(X(idx,:)>M_rang_r).*M_rang_r...
                        +((X(idx,:)>=M_rang_l)&(X(idx,:)<=M_rang_r)).*X(idx,:);
					% Evaluate fitness value
                    fv = feval(fhd,X(idx,:)',varargin{:});
                    fitcount=fitcount+1;
					% Store the pbest
                    if fv < pbestval(idx,1)  
                        pbestval(idx,1)=fv;
                        pbestpos(idx,:)=X(idx,:);
                    end
                end
            end

            [gbestval,gbestid]=min(pbestval);
            if gbestval< best_fitness(n_iteration-1,1)
                best_fitness(n_iteration,1)=gbestval;
                X_gbest=pbestpos(gbestid,:);
            else
                best_fitness(n_iteration,1)=best_fitness(n_iteration-1,1);
            end
            % Count the number of global optimal stop improvements.
            if best_fitness(n_iteration,1) < gbestval_o
               gbestpos_o=X_gbest;
               gbestval_o=best_fitness(n_iteration,1);
               stop=0;
            else
               stop=stop+1;
            end
            if fitcount>=max_FES
                break;
            end

            if (n_iteration==max_iteration)&&(fitcount<max_FES)
                n_iteration=n_iteration-1;
            end 
        end
		
		% OUTPUTS
		
		allgbestval=best_fitness;
		gbestval=best_fitness(n_iteration,1);
		X_gbest;
		Error=gbestval-f_bias(varargin{:});
		fitcount;
end
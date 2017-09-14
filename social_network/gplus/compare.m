p=genpath(pwd);
addpath(p);

load qualified;

d_within1=[];
d_within2=[];
d_between1=[];
d_between2=[];

 %choice = 'normal_simple';
% choice = 'normal_ias';
% choice = 'unnorm_simple';
% choice = 'unnorm_ias';

choices = {'normal_simple','normal_ias'};

count = zeros(length(choices),2); %(i_choice, p, Id) 2 is for p=1;p=2
count_bad = 0;

for i_choice = 1:2
	choice = choices{i_choice};
	 for Id = 1:length(qualified)
	%------------------------------load data-------------------------------------
		fileID= qualified(Id);
		fprintf(1,'\n-------------------------------');
		fprintf(1,'\nworking on file           %s\n',char(fileID));
		load(sprintf('matlab_data/new_feat/%s_feat',char(fileID)),'feat');
		load(sprintf('matlab_data/new_ias/%s_ias',char(fileID)),'ias');
		load(sprintf('matlab_data/new_circles/%s_circles',char(fileID)),'circles');
		load(sprintf('matlab_data/disjoint/%s_disjoint',char(fileID)),'disjoint');
		for p =1:2
			d_between=[];
			%first compute the distances between the circles
			for i=1:length(disjoint)-1
				for k =i+1:length(disjoint)
					i1 = disjoint(i);%index for a circle
					i2 = disjoint(k);%index for another circle				
					circle1=circles{i1}; %indices for the members in a circle
					circle2=circles{i2};			
					%compute distance matrices between circle1 and circle2
					for j=1:length(circle1)
						ID1 = circle1(j);
						ID1 = find(feat(:,1) == ID1 );
						feat1=feat(ID1,2:end);
						
						if (sum(feat1)==0)
							continue;
						end
						
						for l=1:length(circle2)
							ID2 = circle2(l);
							ID2 = find(feat(:,1) == ID2);
							feat2=feat(ID2,2:end);
							
							if (sum(feat2)==0)
								continue;
							end
							
							switch choice
								case 'normal_ias'
									d_between=[d_between,d_normal(feat(ID1,2:end),feat(ID2,2:end),p,ias)]; %feat(:,1) stores the ids of people
								case 'normal_simple'
									d_between=[d_between,d_normal_k(feat(ID1,2:end),feat(ID2,2:end),p)]; %feat(:,1) stores the ids of people
								case 'unnorm_simple'
									d_between=[d_between,d_unnorm_k(feat(ID1,2:end),feat(ID2,2:end),p)]; %feat(:,1) stores the ids of people
								case 'unnorm_ias'
									d_between=[d_between,d_unnorm(feat(ID1,2:end),feat(ID2,2:end),p,ias)]; %feat(:,1) stores the ids of people
							end %end switch
							
						end %end l
					end %end j
				end %end k
			end %end i

			%now compute the distances within the circles
			d_within=[];			
			for i=1:length(disjoint)-1
				i1 = disjoint(i);%index for a circle
				circle=circles{i1}; %indices for the members in a circle
				for j=1:length(circle)
					ID1 = circle(j);
					ID1 = find(feat(:,1) == ID1 );
					feat1= feat(ID1, 2:end);
					
					if ( sum(feat1)==0)
						continue;
					end
					
					for l=i+1:length(circle)
						ID2 = circle(l);
						ID2 = find(feat(:,1) == ID2 );
						feat2 = feat(ID2,2:end);
						
						if (sum(feat2)==0)
							continue;
						end
		
						switch choice
							case 'normal_ias'
								d_temp=d_normal(feat(ID1,2:end),feat(ID2,2:end),p,ias);
							case 'normal_simple'
								d_temp=d_normal_k(feat(ID1,2:end),feat(ID2,2:end),p);
							case 'unnorm_simple'
								d_temp=d_unnorm_k(feat(ID1,2:end),feat(ID2,2:end),p);
							case 'unnorm_ias'
								d_temp=d_unnorm(feat(ID1,2:end),feat(ID2,2:end),p,ias);
						end %end switch
						
						d_within=[d_within,d_temp];
					end %end l
				end %end j
			end	 %end i

			len1=length(d_between);
			len2=length(d_within);
			distance = horzcat(d_between,d_within);
			
			fprintf(1,'len1=%d , len2 = %d\n',len1,len2);
			
			%check if the samples are too small
			if(len1<=2 | len2 <=2)
				fprintf(1,'sample size is too small\n');
				if(p==1)
					count_bad = count_bad +1;
				end
				continue
			elseif (sum(d_between)==0 | sum(d_within)==0)
				fprintf(1,'At least one of d_between and d_within is all zeros\n');
				if(p==1)
					count_bad = count_bad +1;
					printf(1,'count_bad = %d\n',count_bad);
				end
				continue
			end
			
			%get the count
			if median(d_between)>median(d_within)
				count(i_choice,p) = count(i_choice,p) +1;
			end
		
			a=zeros(1,len1);
			b=ones(1,len2);
			indicators = horzcat(a,b);

			corr = corrcoef(distance,indicators);
			
			fprintf(1,'p=%d.       correlation: %f\n',p,corr(1,2));
		    % fprintf(1,'p=%d choice = %s. median comparison: %d(between) vs. %d (within)\n', p,choice,count(i_choice,p),Id- count(i_choice, p));
			
		end%for p
		
		    fprintf(1,'p=%d choice = %s. median comparison: %d(between) vs. %d (within)\n', 1,choices{i_choice},count(i_choice,1),Id-count_bad- count(i_choice, 1));
		    fprintf(1,'p=%d choice = %s. median comparison: %d(between) vs. %d (within)\n', 2,choices{i_choice},count(i_choice,2),Id-count_bad- count(i_choice, 2));
	end%for Id
end%for i_choice
	
	
for i_choice = 1:2
	for p =1:2
		fprintf(1,'\n p=%d choice = %s\nmedian comparison: %d(between) vs. %d (within) \n', p,choices{i_choice},count(i_choice,p),Id-count_bad -count(i_choice,p));
	end
end
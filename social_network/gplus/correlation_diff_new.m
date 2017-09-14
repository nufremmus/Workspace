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

difference = zeros(2,2,length(qualified)); %(i_choice, p, Id)

for i_choice = 1:2
	choice = choices{i_choice};
	% for Id = 1:length(qualified)
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
			
			%get the count
			if median(d_between)>median(d_within)
				difference(i_choice, p,Id) = median(d_between)- median(d_within);
			end
		
			a=zeros(1,len1);
			b=ones(1,len2);
			indicators = horzcat(a,b);
		
			if(a*b > 0)
				corr = corrcoef(distance,indicators);
			elseif(a*b ==0)
				corr = NaN;
			end
			fprintf(1,'p=%d.       correlation: %f\n',p,corr(1,2));
			% fprintf(1,'p=%d. median comparison: %d(between) vs. %d (within)\n', p,sum(count(p,:)),Id- sum(count(p,:)));
			
			% boxplot(distance,indicators);
			% saveas(gca,sprintf('boxplot_%s_%s_p%d.eps',char(fileID),choice,p),'epsc');
		end%for p
		
		% clear disjoint feat circles
		%disjoint feat circles get rewritten for next for loop for ID
	end%for Id
end%for i_choice

save difference difference
 
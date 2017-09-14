%generate the top 10 boxplots from gplus and facebook combined
%put the paired boxplots in one plot using subplot
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

load difference_new;
[r,c,v] = ind2sub(size(difference),find(difference > 0.0331 & difference<= 0.1668));
OK = [r,c,v];

count = 1; % count the current plots

%for subplotting
figure

for k=1:length(OK)

%get the parameters
	Id = v(k);
	% if(ismember(Id,deletedId))
		% continue;
	% end

	i_choice = r(k);
	choice = choices{i_choice};
	
	p=c(k);

%---------------------------start-----------------------------------------
	fileID= qualified(Id);
	fprintf(1,'\n-------------------------------');
	fprintf(1,'\nworking on file           %s\n',char(fileID));
	load(sprintf('matlab_data/feat/%s_feat',char(fileID)),'feat');
	load(sprintf('matlab_data/ias/%s_ias',char(fileID)),'ias');
	load(sprintf('matlab_data/circles/%s_circles',char(fileID)),'circles');
	load(sprintf('matlab_data/disjoint/%s_disjoint',char(fileID)),'disjoint');

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
	
	a=zeros(1,len1);
	b=ones(1,len2);
	indicators = horzcat(a,b);

	corr = corrcoef(distance,indicators);

	fprintf(1,'p=%d.       correlation: %f\n',p,corr(1,2));
	% fprintf(1,'p=%d. median comparison: %d(between) vs. %d (within)\n', p,sum(count(p,:)),Id- sum(count(p,:)));
	
	%-----------subplot------------
	subplot(2,5,count)
	count = count+1;
	 boxplot(distance,indicators,'labels',{'B','W'});
	 % saveas(gca,sprintf('boxplot_%s_%s_p%d.eps',char(fileID),choice,p),'epsc');

	% clear disjoint feat circles
	%disjoint feat circles get rewritten for next for loop for ID
end%end k

  saveas(gca,'boxplot_top10.eps','epsc');
% fprintf(1,'Median comparison: %d(between) vs. %d (within) for p=%d\n', sum(count(1,:)),length(qualified),1);
% fprintf(1,'Median comparison: %d(between) vs. %d (within) for p=%d\n', sum(count(2,:)),length(qualified),2);
 
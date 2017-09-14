function flag=InCircles(x,y,Circles)
%x and y are id's of people
%Circles is a cell file that saves different circles as a vector
[m,n] = size(Circles);
flag =0;
for i=1:m
	if find(Circles{i}==x)
		if find(Circles{i}==y)
			flag = 1;
			break;
		end
	end
end
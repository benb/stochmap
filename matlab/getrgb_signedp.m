function rgbvals=getrgb_signedp(values,signs)
%get rgb colours for a vector of signed p-values
%SYNTAX
%rgbvals=getrgb_signedp(values,signs)
%INPUTS
%values is a vector of p-values
%signs is a vector of signs for the test scores
%OUTPUTS
%rgbvals is an array with one row per original value, and three columns:
%red, green and blue, each on a scale of 0 to 255.
%negative values go from blue (0) to black (-1)
%positive values go from red (0) to black (1)

rgbvals=zeros(length(values),3);

rgbvals(:,2)=0;%constant green component
ind=(signs==-1);

%red
rgbvals(ind,1)=0;
rgbvals(~ind,1)=round(255*(1-values(~ind)));

%blue
rgbvals(~ind,3)=0;
rgbvals(ind,3)=round(255*(1-abs(values(ind))));

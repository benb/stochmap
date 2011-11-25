function rgbvals=getrgb(values,range)
%get rgb colours for a vector of values
%SYNTAX
%rgbvals=getrgb(values,range)
%INPUTS
%values is a vector of values
%range (optional) is the limits of the range. If not specified, use min and
%max from data)
%OUTPUTS
%rgbvals is an array with one row per original value, and three columns:
%red, green and blue, each on a scale of 0 to 255. Red increases, green
%decreases, and blue is constant, as values increase (cf the predefined
%'cool' colormap)

if nargin<2 | isempty(range)
    mn=min(values);
    mx=max(values);
else
    mn=range(1);
    mx=range(2);
end
rgbvals=zeros(length(values),3);

%truncate at the specified range
values(values<mn)=mn;
values(values>mx)=mx;

rgbvals(:,3)=255;%constant blue component
rgbvals(:,1)=round(255*(values-mn)/(mx-mn));
rgbvals(:,2)=round(255*(1-(values-mn)/(mx-mn)));
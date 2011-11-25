function colortreeedges2(edgelengths,values,treefile,outfile,range,sign,exclude,colorchoice,doFDR,threshold)
%colour edges on a Dendroscope tree using values. Should really write this
%in perl and pay more attention to error conditions e.g. edge lengths don't
%match exactly
%This version: supply two different thresholds, giving an output file with
%three different edge widths
%SYNTAX
%colortreeedges2(edgelengths,values,treefile,outfile,range,sign,exclude,colorchoice,doFDR,threshold)
%INPUTS
%edgelengths is a vector of the lengths of edges on the tree
%values is a vector of values by which to colour the edges (in the same
%order)
%treefile is a dendroscope file containing the tree structure, saved after
%using Show Edge Weights (we look for edge lengths in the tree file that
%match the ones in edgelengths in order to put the colours in the right
%places)
%outfile is the file to write to, which will contain the edge colours
%range (optional) is the range of values we use to calculate colours. If
%any values lie outside the range, they are truncated. If not specified,
%use min and max of values
%sign (optional) is signs for the colours (e.g. to distinguish positive and
%negative scores when using p-values)
%exclude (optional) is a logical vector with 1 for edges to exclude (colour
%black)
%colorchoice (optional) is zero if we want edges coloured according to values
%in a color map similar to the predefined matlab 'cool' color map
%or one if we want to emphasize values close to zero (redder if positive, bluer if
%negative): use for p-values
%doFDR (optional) is 1 (default) if we want to do False Discovery Rate
%calculation on p-values (threshold 0.1), and set edges to be thicker if significant under
%this criterion
%threshold (optional, default [0.05,0.1]) is pair of thresholds for FDR
%OUTPUTS
%A file readable by dendroscope, with edges coloured, wider for edges significant at FDR threshold 2,
%widest for threshold 1
%USES
%getrgb.m
%getrgb_signedp.m
%testfdr.m

if nargin<5
    range=[min(values) max(values)];
end
if nargin<6
    sign=[];
end
if nargin<7
    exclude=[];
end
if nargin<8
    colorchoice=0;
end
if nargin<9
    doFDR=1;
end
if nargin<10
    threshold=[0.05 0.1];
end
if doFDR
   FDRvec1=testfdr(values,threshold(1));%edges significant under FDR of threshold 1
   FDRvec2=testfdr(values,threshold(2));%edges significant under FDR of threshold 2
end

if ~isempty(sign)
   values=values.*sign;
end
if ~colorchoice
    rgbvals=getrgb(values,range);%calculate the RGB values
else
    rgbvals=getrgb_signedp(values,sign);%calculate the RGB values: emphasize low values
end
rgbvals(edgelengths==0,:)=round(0.65*255);%grey if edge length zero
if ~isempty(exclude)
    rgbvals(exclude,:)=round(0.65*255);%colour excluded edges grey
end
nedges=length(edgelengths);
fid=fopen(treefile,'r');
fout=fopen(outfile,'w');
matched=0;
while ~feof(fid)%scan the file until we find the line 'edges'
    line=fgetl(fid);
    if strcmpi(line,'edges')
        fprintf(fout,[line,'\n']);
        while matched<nedges && ~feof(fid)
            line=fgetl(fid);
            [tok,rem]=strtok(line,'''');
            [tok,rem]=strtok(rem,'''');
            treel=textscan(tok,'%f');%match the edge length
            treel=treel{1};
            ind=find(edgelengths==treel);
            if length(ind)>1
                disp(['warning: edge length ',num2str(treel),' not unique']);
            elseif length(ind)==1
                matched=matched+1;%write the colour
                
                %is there an existing colour statement?
                oldcol=textscan(line,'%s fg=%d %d %d %s');
                if ~isempty(oldcol{2})
                   line=strrep(line,['fg=',num2str(oldcol{2}),' ',num2str(oldcol{3}),' ',num2str(oldcol{4})],'');%delete existing colour
                end
                if doFDR%use edge thickness for FDR
                    %is there an existing edge thickness statement?
                    oldw=textscan(line,'%s w=%d %s');
                    if ~isempty(oldw{2})
                       line=strrep(line,['w=',num2str(oldw{2})],'');%delete existing edge thickness
                    end
                    mywidth=2+2*FDRvec1(ind)+2*FDRvec2(ind);%lines get thicker as we get significant at lower FDR
                    myline=strrep(line,':',[': fg=',num2str(rgbvals(ind,:)),' w=',num2str(mywidth)]);%write in new colour and edge thickness
                else
                    myline=strrep(line,':',[': fg=',num2str(rgbvals(ind,:))]);%write in new colour
                end
                fprintf(fout,[myline,'\n']);
            else
                disp(['warning: edge length ',num2str(treel),' not matched']); 
            end
        end
    else
        fprintf(fout,[line,'\n']);
    end
end
fclose(fid);
fclose(fout);
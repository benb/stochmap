function [rawdata,branchdata,sitedata,score]=read_stochmap_out2(fname)
%read output from stochmap
%SYNTAX
%[rawdata,branchdata,sitedata,score]=read_stochmap_out2(fname)
%INPUTS
%fname is the name of the .out file
%OUTPUTS
%rawdata is a matrix of raw site x branch data with columns branch,
%process, site, conditional expected number of changes, prior expected
%number of changes, prior variance, z-score ((branch
%conditional-branch prior)/branch prior st dev)
%branchdata is similar, but for events summed over branches, with cols
%branch, conditional, prior, zscore, branch length
%sitedata is similar, but for events summed over sites, with cols site,
%conditional, prior, zscore
%score is scalar overall score (sum of site conditional - prior)

%first we'll find out how many sites and branches we have
fid=fopen(fname,'r');%open file
sf=0;
while (~feof(fid) && ~sf)
   line=fgetl(fid);
   [tok,rem]=strtok(line);
   if strcmp(tok,'Branch')
       sf=1;
   end
end
done=0;
while (~feof(fid) && ~done)
    line=fgetl(fid);%read raw site x branch data: branch, process, site, conditional, prior, variance, zscore
    if isempty(line)
        done=1;
    else%read the data
        A=sscanf(line,'%d %d %d');
    end
end
fclose(fid);%close file

%set branches and sites to max values read so far
nbranch=A(1)+1;
nsite=A(3)+1;

%now read the data in a second pass
fid=fopen(fname,'r');%open file
sf=0;
while (~feof(fid) && ~sf)
   line=fgetl(fid);
   [tok,rem]=strtok(line);
   if strcmp(tok,'Branch')
    sf=1;
   end
end
done=0;
rawdata=zeros(nbranch*nsite,7);
c=1;
while (~feof(fid) && ~done)
    line=fgetl(fid);%read raw site x branch data: branch, process, site, conditional, prior, variance, zscore
    if isempty(line)
        done=1;
    else%read the data
        rawdata(c,:)=sscanf(line,'%f %f %f %f %f %f %f');
        c=c+1;
    end
end

done=0;
branchdata=zeros(nbranch,5);
c=1;
fgetl(fid);
%read branch sums
while (~feof(fid) && ~done)
    line=fgetl(fid);%read branch, conditional, prior, zscore, branch length
    if isempty(line)
        done=1;
    else%read the data
        branchdata(c,:)=sscanf(line,'%f %f %f %f %f');
        c=c+1;
    end
end

done=0;
sitedata=zeros(nsite,4);
c=1;
fgetl(fid);
%read site sums
while (~feof(fid) && ~done)
    line=fgetl(fid);%read site, conditional, prior, zscore
    if isempty(line)
        done=1;
    else%read the data
        sitedata(c,:)=sscanf(line,'%f %f %f %f');
        c=c+1;
    end
end

score=sum(sitedata(:,2)-sitedata(:,3));%sum of conditional - prior over sites

fclose(fid);%close file
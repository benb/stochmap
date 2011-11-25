function [gapcvectplus,gapcvectminus,seq,gapsigns]=addgaps2(fname,target,cvect,gapchar,signs)
%add gaps (NaN) to colour vector for structure representation, where the
%structure comes from a different sequence to the colour vector data
%This version: two separate colour vectors, one for positives and one for
%negatives
%WARNING: assumes there are no gaps in the aligned sequence from which the
%structure came (which will probably be the case, but should check)
%SYNTAX
%[gapcvectplus,gapcvectminus,seq,gapsigns]=addgaps2(fname,target,cvect,gapchar,signs)
%INPUTS
%fname is a .aln file from a profile alignment in clustalx, in which the
%sequence corresponding to the structure was added to a pre-existing
%alignment which was used to obtain the colour data
%target is the name of the target species from which we will add gaps to
%the colour data
%cvect is a vector of ungapped colour data (e.g. p values, range [0,1])
%gapchar is the gap character in the alignment
%signs is a vector of ungapped signs for colour data (-1 and +1)
%OUTPUTS
%gapcvectplus is a vector of gapped colour data [range 0,1] for positive signs, with NaN where there are gaps
%in the alignment for the target species
%gapcvectminus is a vector of gapped colour data [range 0,1] for negative signs, with NaN where there are gaps
%in the alignment for the target species
%seq is the target sequence with gaps
%gapsigns is the gapped signs

missval=99;%value we use to code gaps and opposite-sign values

seq=[];
%open the file
fid=fopen(fname,'r');
while ~feof(fid)
    %read lines from target species
    line=fgetl(fid);
    [seqname,sl]=strtok(line);
    if strcmpi(seqname,target)%line from our target species
       seq=[seq,strtrim(sl)];%append to sequence, after trimming leading and trailing white space 
    end
end
fclose(fid);

%create vector to store colour values
gapcvectplus=missval*ones(size(seq'));

%put NaN where there are gaps in target species
notgaps=(seq~=gapchar);
gapcvectplus(notgaps)=cvect;
gapcvectminus=gapcvectplus;

%put NaN where sign is opposite to what we want
gapsigns=NaN*ones(size(seq'));
gapsigns(notgaps)=signs;%put gaps into signs
gapcvectplus(gapsigns==-1)=missval;
gapcvectminus(gapsigns==1)=missval;


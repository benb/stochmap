reinitialize

#load the structure
cmd.load("1IJF.pdb")

hide all
show cartoon

#read in the colour values for positive scores
stored.newBplus=[]
inFile=open("siteshiftplus.txt",'r')
for line in inFile.readlines(): stored.newBplus.append(float(line))
inFile.close()

alter 1IJF,b=0.0
alter 1IJF and n. CA and chain A, b=stored.newBplus.pop(0)
spectrum b

set_color lgray= [0.90 , 0.90 , 0.90]
set_color dgray= [0.50 , 0.50 , 0.50]


color dgray, chain A #so that we see dark grey for residues with no shift info

#colour all positive sites pink
select pos, chain A and name CA and b<99.0
color lightpink, pos

#read in the colour values for negative scores
stored.newBminus=[]
inFile=open("siteshiftminus.txt",'r')
for line in inFile.readlines(): stored.newBminus.append(float(line))
inFile.close()

alter 1IJF,b=0.0
alter 1IJF and n. CA and chain A, b=stored.newBminus.pop(0)

#colour all negative sites cyan
select neg, chain A and name CA and b<99.0
color palecyan, neg

#distinguish the EF1-beta in yellow
color paleyellow, chain B

#rotate the view
set_view (\
    -0.633794725,   -0.418469369,    0.650525153,\
    -0.580507398,   -0.298458308,   -0.757577658,\
     0.511179328,   -0.857783258,   -0.053767774,\
     0.000000000,    0.000000000, -255.919601440,\
    20.852466583,   34.344417572,   27.699907303,\
   201.768966675,  310.070251465,    0.000000000 )


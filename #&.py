# -*- coding: cp1253 -*-
from Bio import SeqIO
handle = open("C:\Users\Αναστασία\Desktop\Daily_work\Raw_materials_for_samples\disprot.fasta","rU")
disordered = list(SeqIO.parse(handle,"fasta"))
handle.close()


#find the value of the Net Charge of each amino acid

def NetCharge(sequence):
    charge=[]
    for k in sequence:
        if k=="K" or k=="R": charge.append(-1)
        elif k=="D" or k=="E": charge.append(1)
        else: charge.append(0)

    return charge


#find the hydrophobicity for an amino acid, based on K/D scale

kd={'L': 3.8, 'S': -0.8,'E': -3.5, 'K': -3.9, 'Q': -3.5, 'C': 2.5, 'F': 2.8, 'T': -0.7, 'H': -3.2, 'I': 4.5, 'M': 1.9, 'Y': -1.3, 'D': -3.5, 'A': 1.8, 'V': 4.2, 'W': -0.9, 'G': -0.4, 'P': -1.6, 'N': -3.5, 'R': -4.5}

def Hydrophobicity(sequence):
    h=[]

    for residue in sequence:
        h.append(kd[residue])

    return h




#finds the sequence complexity for a window of N amino acids, based on Shannon's Entropy
import math

def Shannon(sequence, N): #N is the given window
    shannon=[]
                            #finds the frequency of each amino acid, for a window of N
    length=len(sequence)    #used for Shannon entropy computation
    freq=[]                 
    position=0
    for aa in sequence:
        freq.append(frequency_of_each_aminoacid(aa, position, sequence, length, N))
        position=position + 1

        
    length=len(sequence)
    
    for position in range(length):
        start=0
        end= length - 1
        if position-N>start: start= position-N
        if position+N<end: end= position+N

        K=0
        for k in freq[start:end]:
            if k!=0:
               K=K-(k*math.log(k,2))
        K = round(K,3) #returns only the first 3 digits after "."
        shannon.append(K)
        
    return shannon


def frequency_of_each_aminoacid(aa, position, sequence, length, N):
    f=0
    start=0
    end = length - 1

    if position-N>start: start= position-N
    if position+N<end: end= position+N

    for k in range(start,end+1):
        if sequence[k]== aa:
            f=f+1

    aminoacid_number= end-start+1
    return round(1.0 * f / aminoacid_number , 3) #returns only the first 3 digits after "."





#finds the frequency of all 20 amino acids in a window of N for the given sequence 

def Frequency(sequence , N):
    f={'L': [], 'S': [],'E': [], 'K': [], 'Q': [], 'C': [], 'F': [], 'T': [], 'H': [], 'I': [], 'M': [], 'Y': [], 'D': [], 'A': [], 'V': [], 'W': [], 'G': [], 'P': [], 'N': [], 'R': []}
    length=len(sequence)
    position=0
    
    for aa in sequence:
        start=0
        end = length - 1
        if position-N>start: start= position-N
        if position+N<end: end= position+N

        for i in f:
            f[i].append(window(i, sequence, start, end))

        position=position + 1

    return f

def window(aa, sequence, start, end):

    s=0
    for k in range(start,end+1):
        if sequence[k]== aa:
            s=s+1
    return round(1.0 * s / (end-start+1) , 3)


	
	
	
#hydrophobic cluster 1st version
#f.e.   ... 0 0 0 0.125 0.25 0.5 1 1 1 1 0.5 0.25 0.125 0 0 0 ...

def Cluster(sequence):
    c=[]
    for k in sequence:
        if k=="V" or k == "I" or k == "L" or k == "F" or k == "M" or k == "Y" or k == "W":
            c.append(1)
        elif k=="P":
            c.append(2)
        else: c.append(0)
               
    start=[]
    end=[]
    start_check=False

    for position in range (1,len(c)):
        if  start_check==False and (c[position-1]==2 or(position>3 and c[position-1]==0 and c[position-2]==0 and c[position-3]==0 and c[position-4]==0)):
            if position<len(c)-5 and c[position]!=2 and c[position+1]!=2 and c[position+2]!=2 and c[position+3]!=2 and (c[position]==1 or c[position+1]==1 or c[position+2]==1 or c[position+3]==1): 
                start.append(position)
                start_check=True

        if  start_check==True and ((position<len(c)-1 and c[position+1]==2) or(position<len(c)-5 and c[position+1]==0 and c[position+2]==0 and c[position+3]==0 and c[position+4]==0)):
            end.append(position)
            start_check= False

    if len(end)!= len(start):
        end.append(len(c)-1)
  
    l=[]
    for k in range (len(c)):
        l.append("X")
    for k in range(len(start)):
        if end[k]-start[k]>2:
            for j in range(start[k],end[k]+1):
                l[j]=1

    for k in range(len(l)-3):
        if l[k]==1 and l[k+1]!=1:
            l[k+1]=0.5
            if l[k+2]=='X' or l[k+2]==0.125:
                l[k+2]=0.25
                if l[k+3]=='X':
                    l[k+3]=0.125
    
    for k in range (3,len(l)):
        if l[k]==1 and l[k-1]!=1:
            l[k-1]=0.5
            if l[k-2]=='X' or l[k-2]==0.125:
                l[k-2]=0.25
                if l[k-3]=='X':
                    l[k-3]=0.125

    for k in range(len(l)):
        if l[k]=='X':
            l[k]=0        
    return l



#hydrophobic cluster 2st version
#f.e. ... 0 0 0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 ...
def Cluster2(sequence):
    c=[]
    for k in sequence:
        if k=="V" or k == "I" or k == "L" or k == "F" or k == "M" or k == "Y" or k == "W":
            c.append(1)
        elif k=="P":
            c.append(2)
        else: c.append(0)
               
    start=[]
    end=[]
    start_check=False

    for position in range (1,len(c)):
        if  start_check==False and (c[position-1]==2 or(position>3 and c[position-1]==0 and c[position-2]==0 and c[position-3]==0 and c[position-4]==0)):
            if position<len(c)-5 and c[position]!=2 and c[position+1]!=2 and c[position+2]!=2 and c[position+3]!=2 and (c[position]==1 or c[position+1]==1 or c[position+2]==1 or c[position+3]==1): 
                start.append(position)
                start_check=True

        if  start_check==True and ((position<len(c)-1 and c[position+1]==2) or(position<len(c)-5 and c[position+1]==0 and c[position+2]==0 and c[position+3]==0 and c[position+4]==0)):
            end.append(position)
            start_check= False

    if len(end)!= len(start):
        end.append(len(c)-1)
  
    l=[]
    for k in range (len(c)):
        l.append(0)
    for k in range(len(start)):
        if end[k]-start[k]>2:
            for j in range(start[k],end[k]+1):
                l[j]=1
                
    return l







# uses the symbols # for disordered and & for ordered regions
# from the disordered.description
# for the classification of each aminoacid
def classify(Name):
    for k in disordered:
        if k.name==Name:
            print Name###
            d=k.description
            d=d.split(' ')
            break
    clas=['X']*len(k.seq)
    
    for j in d:
        if j[0]=='&':
            temp=j.replace('&','')
            temp=temp.split('-')
            for t in range(int(temp[0])-1,int(temp[1])):
                clas[t]=0                                
        elif j[0]=='#':
            temp=j.replace('#','')
            temp=temp.split('-')
            for t in range(int(temp[0])-1,int(temp[1])):
                clas[t]=1
    return clas




for k in range (len(disordered)):
    disordered[k].name= (disordered[k].name).replace('|','_')
    disordered[k].name= (disordered[k].name).replace('-','_')
    disordered[k].name= (disordered[k].name).replace('/','_')


#check what exists in the NetSurf file
#for each aminoacid in each sequence, compute the values from my functions
c= open('C:\Users\Αναστασία\Desktop\Daily_work\Raw_materials_for_samples\surf.txt', 'r')
output1=open('C:\Users\Αναστασία\Desktop\Daily_work\Raw_materials_for_samples\surf3.txt', 'w+')

#output1.write('Protein_Position'+ '\t' +  'Amino Acid' + '\t' + 'Net Charge' + '\t' + 'Hydrophobicity'  + '\t' + 'Shannon Entropy' + \
#    '\t' + 'Frequency A'  + '\t' + 'Frequency C'  + '\t' + 'Frequency D'  + '\t' + 'Frequency E'  + '\t' + \
#      'Frequency F'  + '\t' + 'Frequency G'  + '\t' + 'Frequency H'  + '\t' + 'Frequency I'  + '\t' + 'Frequency K'  + '\t' + 'Frequency L'  + \
#    '\t' +'Frequency M'  + '\t' + 'Frequency N'  + '\t' + 'Frequency P'  + '\t' + 'Frequency Q'  + '\t' + 'Frequency R'  + '\t' + 'Frequency S'  + \
#   '\t' +'Frequency T'  + '\t' + 'Frequency V'  + '\t' + 'Frequency W'  + '\t' + 'Frequency Y'+ \
#   '\t' +'Cluster 1' + '\t' +'Cluster 2' + '\t' +'B=0/E=1' + '\t' + 'RSA' +  '\t' + 'ASA' + '\t' + \
#      'Alpha-Helix' +  '\t' + 'Beta-strand' +  '\t' 'Coil' +  '\t' 'Ordered=0/Disordered=1' + '\n\n')

Buried_Exposed=[]
Aminoacid=[]
start=False

Position=[]
RSA=[]
ASA=[]

Alpha_Helix=[]
Beta_Strand=[]
Coil=[]

windowN1=20
windowN2=25
windowN3=30
windowN4=35
windowN5=40

for line in c:
    word=line.split(' ')
    b=filter(None,word)

    if start==False:
        Name=b[2]
        start=True

    
    if b[2]!=Name:
        net_charge=NetCharge(Aminoacid)
        hydrophobicity=Hydrophobicity(Aminoacid)
        
        shannon1=Shannon(Aminoacid, windowN1)
        shannon2=Shannon(Aminoacid, windowN2)
        shannon3=Shannon(Aminoacid, windowN3)
        shannon4=Shannon(Aminoacid, windowN4)
        shannon5=Shannon(Aminoacid, windowN5)
        
        frequency1=Frequency(Aminoacid, windowN1)
        frequency2=Frequency(Aminoacid, windowN2)
        frequency3=Frequency(Aminoacid, windowN3)
        frequency4=Frequency(Aminoacid, windowN4)
        frequency5=Frequency(Aminoacid, windowN5)
        
        cluster=Cluster(Aminoacid)
        cluster2=Cluster2(Aminoacid)

        classification=classify(Name)
        for index in range (len(net_charge)):
            if classification[index]!='X':
                space=5-len(Position[index])# I want to print the values aligned!!!
                output1.write(Name +"_"+ str(Position[index])+ (" ")*space + '\t'+ Aminoacid[index]  +'\t'+ str(net_charge[index]) +'\t' +str(hydrophobicity[index]) +'\t' + str(cluster[index]) +'\t' + str(cluster2[index]) +'\t' + str(Buried_Exposed[index]) +'\t' + \
                              str(shannon1[index]) +'\t' + str(shannon2[index]) +'\t' + str(shannon3[index]) +'\t' + str(shannon4[index]) +'\t' + str(shannon5[index]) +'\t' + \
                              str(frequency1['A'][index]) +'\t' + str(frequency2['A'][index]) +'\t' + str(frequency3['A'][index]) +'\t' + str(frequency4['A'][index]) +'\t' + str(frequency5['A'][index]) +'\t' + \
                              str(frequency1['C'][index]) +'\t' + str(frequency2['C'][index]) +'\t' + str(frequency3['C'][index]) +'\t' + str(frequency4['C'][index]) +'\t' + str(frequency5['C'][index]) +'\t' + \
                              str(frequency1['D'][index]) +'\t' + str(frequency2['D'][index]) +'\t' + str(frequency3['D'][index]) +'\t' + str(frequency4['D'][index]) +'\t' + str(frequency5['D'][index]) +'\t' + \
                              str(frequency1['E'][index]) +'\t' + str(frequency2['E'][index]) +'\t' + str(frequency3['E'][index]) +'\t' + str(frequency4['E'][index]) +'\t' + str(frequency5['E'][index]) +'\t' + \
                              str(frequency1['F'][index]) +'\t' + str(frequency2['F'][index]) +'\t' + str(frequency3['F'][index]) +'\t' + str(frequency4['F'][index]) +'\t' + str(frequency5['F'][index]) +'\t' + \
                              str(frequency1['G'][index]) +'\t' + str(frequency2['G'][index]) +'\t' + str(frequency3['G'][index]) +'\t' + str(frequency4['G'][index]) +'\t' + str(frequency5['G'][index]) +'\t' + \
                              str(frequency1['H'][index]) +'\t' + str(frequency2['H'][index]) +'\t' + str(frequency3['H'][index]) +'\t' + str(frequency4['H'][index]) +'\t' + str(frequency5['H'][index]) +'\t' + \
                              str(frequency1['I'][index]) +'\t' + str(frequency2['I'][index]) +'\t' + str(frequency3['I'][index]) +'\t' + str(frequency4['I'][index]) +'\t' + str(frequency5['I'][index]) +'\t' + \
                              str(frequency1['K'][index]) +'\t' + str(frequency2['K'][index]) +'\t' + str(frequency3['K'][index]) +'\t' + str(frequency4['K'][index]) +'\t' + str(frequency5['K'][index]) +'\t' + \
                              str(frequency1['L'][index]) +'\t' + str(frequency2['L'][index]) +'\t' + str(frequency3['L'][index]) +'\t' + str(frequency4['L'][index]) +'\t' + str(frequency5['L'][index]) +'\t' + \
                              str(frequency1['M'][index]) +'\t' + str(frequency2['M'][index]) +'\t' + str(frequency3['M'][index]) +'\t' + str(frequency4['M'][index]) +'\t' + str(frequency5['M'][index]) +'\t' + \
                              str(frequency1['N'][index]) +'\t' + str(frequency2['N'][index]) +'\t' + str(frequency3['N'][index]) +'\t' + str(frequency4['N'][index]) +'\t' + str(frequency5['N'][index]) +'\t' + \
                              str(frequency1['P'][index]) +'\t' + str(frequency2['P'][index]) +'\t' + str(frequency3['P'][index]) +'\t' + str(frequency4['P'][index]) +'\t' + str(frequency5['P'][index]) +'\t' + \
                              str(frequency1['Q'][index]) +'\t' + str(frequency2['Q'][index]) +'\t' + str(frequency3['Q'][index]) +'\t' + str(frequency4['Q'][index]) +'\t' + str(frequency5['Q'][index]) +'\t' + \
                              str(frequency1['R'][index]) +'\t' + str(frequency2['R'][index]) +'\t' + str(frequency3['R'][index]) +'\t' + str(frequency4['R'][index]) +'\t' + str(frequency5['R'][index]) +'\t' + \
                              str(frequency1['S'][index]) +'\t' + str(frequency2['S'][index]) +'\t' + str(frequency3['S'][index]) +'\t' + str(frequency4['S'][index]) +'\t' + str(frequency5['S'][index]) +'\t' + \
                              str(frequency1['T'][index]) +'\t' + str(frequency2['T'][index]) +'\t' + str(frequency3['T'][index]) +'\t' + str(frequency4['T'][index]) +'\t' + str(frequency5['T'][index]) +'\t' + \
                              str(frequency1['V'][index]) +'\t' + str(frequency2['V'][index]) +'\t' + str(frequency3['V'][index]) +'\t' + str(frequency4['V'][index]) +'\t' + str(frequency5['V'][index]) +'\t' + \
                              str(frequency1['W'][index]) +'\t' + str(frequency2['W'][index]) +'\t' + str(frequency3['W'][index]) +'\t' + str(frequency4['W'][index]) +'\t' + str(frequency5['W'][index]) +'\t' + \
                              str(frequency1['Y'][index]) +'\t' + str(frequency2['Y'][index]) +'\t' + str(frequency3['Y'][index]) +'\t' + str(frequency4['Y'][index]) +'\t' + str(frequency5['Y'][index]) +'\t' + \
                              str(RSA[index]) +'\t' + str(ASA[index]) +'\t' + str(Alpha_Helix[index]) +'\t' + str(Beta_Strand[index]) +'\t' +str(Coil[index]) +'\t' +str(classification[index]) + '\n')

            
           
        Buried_Exposed=[]
        Aminoacid=[]
        Position=[]
        RSA=[]
        ASA=[]
        Alpha_Helix=[]
        Beta_Strand=[]
        Coil=[]
        Name=b[2]
      
    if b[0]=='B':
        Buried_Exposed.append(0)
    else:
        Buried_Exposed.append(1)
    
    Aminoacid.append(b[1])

    Position.append(b[3])
    RSA.append(b[4])
    ASA.append(b[5])
    Alpha_Helix.append(b[7])
    Beta_Strand.append(b[8])
    Coil.append(b[9].replace('\n',''))
    del b
net_charge=NetCharge(Aminoacid)
hydrophobicity=Hydrophobicity(Aminoacid)

shannon1=Shannon(Aminoacid, windowN1)
shannon2=Shannon(Aminoacid, windowN2)
shannon3=Shannon(Aminoacid, windowN3)
shannon4=Shannon(Aminoacid, windowN4)
shannon5=Shannon(Aminoacid, windowN5)
        
frequency1=Frequency(Aminoacid, windowN1)
frequency2=Frequency(Aminoacid, windowN2)
frequency3=Frequency(Aminoacid, windowN3)
frequency4=Frequency(Aminoacid, windowN4)
frequency5=Frequency(Aminoacid, windowN5)
        
cluster=Cluster(Aminoacid)
cluster2=Cluster2(Aminoacid)

classification=classify(Name)
for index in range (len(net_charge)):
    if classification[index]!='X':
        output1.write(Name +"_"+ str(Position[index])+ (" ")*space + '\t'+ Aminoacid[index]  +'\t'+ str(net_charge[index]) +'\t' +str(hydrophobicity[index]) +'\t' + str(cluster[index]) +'\t' + str(cluster2[index]) +'\t' + str(Buried_Exposed[index]) +'\t' + \
                      str(shannon1[index]) +'\t' + str(shannon2[index]) +'\t' + str(shannon3[index]) +'\t' + str(shannon4[index]) +'\t' + str(shannon5[index]) +'\t' + \
                      str(frequency1['A'][index]) +'\t' + str(frequency2['A'][index]) +'\t' + str(frequency3['A'][index]) +'\t' + str(frequency4['A'][index]) +'\t' + str(frequency5['A'][index]) +'\t' + \
                      str(frequency1['C'][index]) +'\t' + str(frequency2['C'][index]) +'\t' + str(frequency3['C'][index]) +'\t' + str(frequency4['C'][index]) +'\t' + str(frequency5['C'][index]) +'\t' + \
                      str(frequency1['D'][index]) +'\t' + str(frequency2['D'][index]) +'\t' + str(frequency3['D'][index]) +'\t' + str(frequency4['D'][index]) +'\t' + str(frequency5['D'][index]) +'\t' + \
                      str(frequency1['E'][index]) +'\t' + str(frequency2['E'][index]) +'\t' + str(frequency3['E'][index]) +'\t' + str(frequency4['E'][index]) +'\t' + str(frequency5['E'][index]) +'\t' + \
                      str(frequency1['F'][index]) +'\t' + str(frequency2['F'][index]) +'\t' + str(frequency3['F'][index]) +'\t' + str(frequency4['F'][index]) +'\t' + str(frequency5['F'][index]) +'\t' + \
                      str(frequency1['G'][index]) +'\t' + str(frequency2['G'][index]) +'\t' + str(frequency3['G'][index]) +'\t' + str(frequency4['G'][index]) +'\t' + str(frequency5['G'][index]) +'\t' + \
                      str(frequency1['H'][index]) +'\t' + str(frequency2['H'][index]) +'\t' + str(frequency3['H'][index]) +'\t' + str(frequency4['H'][index]) +'\t' + str(frequency5['H'][index]) +'\t' + \
                      str(frequency1['I'][index]) +'\t' + str(frequency2['I'][index]) +'\t' + str(frequency3['I'][index]) +'\t' + str(frequency4['I'][index]) +'\t' + str(frequency5['I'][index]) +'\t' + \
                      str(frequency1['K'][index]) +'\t' + str(frequency2['K'][index]) +'\t' + str(frequency3['K'][index]) +'\t' + str(frequency4['K'][index]) +'\t' + str(frequency5['K'][index]) +'\t' + \
                      str(frequency1['L'][index]) +'\t' + str(frequency2['L'][index]) +'\t' + str(frequency3['L'][index]) +'\t' + str(frequency4['L'][index]) +'\t' + str(frequency5['L'][index]) +'\t' + \
                      str(frequency1['M'][index]) +'\t' + str(frequency2['M'][index]) +'\t' + str(frequency3['M'][index]) +'\t' + str(frequency4['M'][index]) +'\t' + str(frequency5['M'][index]) +'\t' + \
                      str(frequency1['N'][index]) +'\t' + str(frequency2['N'][index]) +'\t' + str(frequency3['N'][index]) +'\t' + str(frequency4['N'][index]) +'\t' + str(frequency5['N'][index]) +'\t' + \
                      str(frequency1['P'][index]) +'\t' + str(frequency2['P'][index]) +'\t' + str(frequency3['P'][index]) +'\t' + str(frequency4['P'][index]) +'\t' + str(frequency5['P'][index]) +'\t' + \
                      str(frequency1['Q'][index]) +'\t' + str(frequency2['Q'][index]) +'\t' + str(frequency3['Q'][index]) +'\t' + str(frequency4['Q'][index]) +'\t' + str(frequency5['Q'][index]) +'\t' + \
                      str(frequency1['R'][index]) +'\t' + str(frequency2['R'][index]) +'\t' + str(frequency3['R'][index]) +'\t' + str(frequency4['R'][index]) +'\t' + str(frequency5['R'][index]) +'\t' + \
                      str(frequency1['S'][index]) +'\t' + str(frequency2['S'][index]) +'\t' + str(frequency3['S'][index]) +'\t' + str(frequency4['S'][index]) +'\t' + str(frequency5['S'][index]) +'\t' + \
                      str(frequency1['T'][index]) +'\t' + str(frequency2['T'][index]) +'\t' + str(frequency3['T'][index]) +'\t' + str(frequency4['T'][index]) +'\t' + str(frequency5['T'][index]) +'\t' + \
                      str(frequency1['V'][index]) +'\t' + str(frequency2['V'][index]) +'\t' + str(frequency3['V'][index]) +'\t' + str(frequency4['V'][index]) +'\t' + str(frequency5['V'][index]) +'\t' + \
                      str(frequency1['W'][index]) +'\t' + str(frequency2['W'][index]) +'\t' + str(frequency3['W'][index]) +'\t' + str(frequency4['W'][index]) +'\t' + str(frequency5['W'][index]) +'\t' + \
                      str(frequency1['Y'][index]) +'\t' + str(frequency2['Y'][index]) +'\t' + str(frequency3['Y'][index]) +'\t' + str(frequency4['Y'][index]) +'\t' + str(frequency5['Y'][index]) +'\t' + \
                      str(RSA[index]) +'\t' + str(ASA[index]) +'\t' + str(Alpha_Helix[index]) +'\t' + str(Beta_Strand[index]) +'\t' +str(Coil[index]) +'\t' +str(classification[index]) + '\n')

c.close()
output1.close()

#!/usr/bin/python3
#-*- coding : utf-8 -*-
#CHANGEMENT POUR TESTER COMMIT


__authors__ = ("Mathilde Chatain", "Lucien Maurau")
__contact__ = ("mathilde.chatain@etu.umontpellier.fr","lucien.maurau@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "11/24/2025"
__licence__ ="This program is free software: you can redistribute it and/or modifyit under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."


     
    ### OPTION LIST:
        ##-h or --help : help information
        ##-i or --input: input file (.sam)
        ##-o or --output: output name files (.txt)

    #Synopsis:
        ##SamReader.py -h or --help # launch the help.
        ##SamReader.py -i or --input <file> # Launch SamReader to analyze a samtools file (.sam) and print the result in the terminal
        ##SamReader.py -i or --input <file> -o or --output <name> # Launch SamReader to analyze a samtools file (.sam) and print the result in the file called <name>
  
############### FUNCTIONS TO :

from itertools import islice
import os;
import sys;
import re


#Function that checks the integrity of the input file
def check_file(file):
    if not os.path.exists(file): #check if the file exists
        print(f"Error : the path '{file}' does not exist.")
        return False

    if not os.path.isfile(file): #check if the path is a file
        print(f"Error : the path '{file}' is not a file.")
        return False
    size = os.path.getsize(file)
    if size == 0: #check if the file is empty
        print(f"Error : the file '{file}' is empty.")
        return False

    print(f"OK: the file '{file}' exists and is not empty (size: {size} bytes)") #confirm the file is valid
    return True

## 2/ Read, 
#Function that reads the file and separates the lines according to a minimum quality value
def SamRead(file,minval):
    with open(file,"r") as Str:  
        L=list(islice(enumerate(Str),2,None))
    return [t[1] for t in L if int(t[1].split("\t")[4]) > minval] #We take 1 to not have the line number
    
## 3/ Store,
#Function that calls check_file and samread and stores the results of unmapped and partiallymapped
def Store(file,minval):
    if check_file(file): #We check if the file is not empty and exists
        Qual=ReadQuality(SamRead2(file, -1))
        sam_lines=SamRead(file,minval) #We obtain a list of lines 
        unmapped_count,un=unmapped(sam_lines) #We obtain the number of unmapped reads and a dictionary containing them
        partially_mapped_count,par=partiallyMapped(sam_lines) #We obtain the number of partiallyMapped reads and a dictionary containing them
        Pm=PerfectMapped(sam_lines) #We obtain a dictionary containing the IDs of perfectly mapped reads (we keep only one occurrence per pair since we are not interested in the line content)
        idd=Ids(file,minval) # We build a dictionary of read IDs with their occurrence counts (to verify that each read is properly paired)

        return unmapped_count,partially_mapped_count,Pm,un,par,idd,Qual
    else:
        return None #Change for an exit
    
def Analyse(Pm,Un,Par,idd,file,minval): # Function that analyzes read pairs and calculates the global percentage of each mutation based on the CIGAR string
    MapUnMap={} #Dictionary that will store mapped–unmapped read pairs
    MapParMap={} #Dictionary that will store Mapped–PartiallyMapped read pairs
    IdLines=SamRead2(file,minval) #Dictionary that contains all lines, with the read ID as the key and the rest of the line as the value
    for s in idd.keys(): #We iterate over each read ID
        if s in Pm and s in Un: #We check if the pair is present in perfectMapped and Unmapped
            MapUnMap[s]=Un[s] #We add the first read of the pair
            MapUnMap[s+"-1"]=Pm[s] #We add the second read of the pair
        if s in Pm and s in Par:  #We check if the pair is present in perfectMapped and partiallyMapped
            MapParMap[s]=Par[s] #We add the first read of the pair
            MapParMap[s+"-1"]=Pm[s] #We add the second read of the pair
        with open("../Results/outpuTable_cigar.txt","a+") as outpuTable: #We will write the mutation percentages for each line
            if idd[s]==2: #Verify that the ID is associated with a read pair
                D1,D2=readCigar(IdLines[s][4]),readCigar(IdLines[s+"-1"][4]) #Obtain the mutation percentages for both reads in the pair
                LTXT=percentMutation(D1)+";"+percentMutation(D2) #Concatenate into a single line
                outpuTable.write(LTXT+"\n") #Write to a TXT file that will contain the results
    return MapUnMap,MapParMap 
    
def ReadQuality(sam_line):
    D={}
    for line in sam_line.values():
        value=int(line[3])
        if value not in D.keys():
            D[value]=1
        else :
            D[value]+=1
    with open("../Results/QualityDistribution.tsv","a+") as outpuTable:
        outpuTable.write("Quality"+",\t"+"Number \n")
        for k in D.keys():
            outpuTable.write(str(k)+"\t"+str(D[k])+"\n")
    return D
            

#Converti les flag en binaire, adapter la lecture en fonction 
def flagBinary(flag):
    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB


#Take the unmapped lines and write a FASTA file containing them, as well as a summary file reporting the number of unmapped reads
def unmapped(sam_line):
    D={} #Dictionnary that will contains unmapped reads
    unmapped_count = 0
    with open ("../Results/only_unmapped.fasta", "a+") as unmapped_fasta, open("../Results/summary_unmapped.txt", "w") as summary_file:
        for line in sam_line: 
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])
            if int(flag[-3]) == 1:
                D[col_line[0]]=col_line[1:]
                unmapped_count += 1
                unmapped_fasta.write(line) #retire "toStringOutput"
        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
        return unmapped_count,D

#Function that counts partially mapped reads by examining the CIGAR string
def partiallyMapped(sam_line):
    D={}
    partially_mapped_count = 0
    with open ("../Results/only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("../Results/summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1]) # We compute the same 
            if int(flag[-2]) == 1: 
                if col_line[5] != str(len(col_line[9]))+"M": #à spécifier 
                    D[col_line[0]]=col_line[1:]
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(line) #retired "toStringOutput"
        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count,D

def PerfectMapped(sam_line):
    D={}
    for line in sam_line:
        col_line = line.split("\t")
        if int(flagBinary(col_line[1])[-3])==0 and col_line[5] == str(len(col_line[9]))+"M":
            D[col_line[0]]=col_line[1:]
        #if not isPartMapped(col_line) and not isunmapped(col_line):
        #    D[col_line[0]]=col_line[1:]
    return D

def isPartMapped(col_line):
    return int(flagBinary(col_line[1])[-2]) ==1 and col_line[5] != "100M"

def isunmapped(col_line):
    return int(flagBinary(col_line[1])[-3]) == 1

def readCigar(cigar): 
    ext = re.findall('\w',cigar) # Split all characters
    key=[] 
    value=[]    
    val=""
    for i in range(0,len(ext)): # For each numeric values or alpha numeric
        if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
            key.append(ext[i]) #Previous if: if a letter is found, add it to the dictionary with the letter as key and the preceding number as value 
            value.append(val) #we add the value
            val = "" # Re-initialize

        else :
            val = "" + val + ext[i]  #Store values while iterating through the CIGAR string
    dico = {}
    n = 0
    for k in key:   #we make the dictionary         
        if k not in dico.keys():    #We do not want the same key twice, so we check each time if it is already present
            dico[k] = int(value[n])   #Assign the value 
            n += 1
        else:
            dico[k] += int(value[n])  #If the key already exists, add to the existing value
            n += 1
    return dico #Returns a dictionary giving the sum of operations, but not their order (which was recorded in the CIGAR)

#Function that analyzes the dictionary created by the readCigar function
def percentMutation(dico):
    totalValue = 0 #total value of mutations
    for v in dico :
        totalValue += dico[v] #Total number of operations recorded in the CIGAR string
    mutList = ['M','I','D','S','H','N','P','X','='] #Bizarre, pourquoi il y a le M qui devrait correspondre à "match", c'est pas une mutation ?? 
    res = "" 
    for mut in mutList : #For each type of mutation, iterate through the dictionary
        if mut in dico.keys() :
            res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";") #Add the percentage of the current operation
        else :
            res += ("0.00" + ";") #Add zero if the mutation type is not present
    return res[:-1] #on met zéro sinon, l'ordre spécifie l'opération (voir fonction suivante), donc pas besoin de nommer à chaque fois

#Function that writes a .txt summarizing the distribution of all CIGAR operations using a file containing the outputs of the previous function (iterated pair by pair?)
def globalPercentCigar():
    nbReads, M, I, D, S, H, N, P, X, Egal = [0 for _ in range(10)]
    with open ("../Results/outpuTable_cigar.txt","r") as outpuTable, open("../Results/Final_Cigar_table.txt", "w") as FinalCigar: #Explicitly open the .txt files for reading and writing
        for line in outpuTable : #Iterate over each line in the input table
            mutValues = line.split(";")
            nbReads += 2 #Total number of reads (2 per line)
            M += float(mutValues[0])+float(mutValues[9])
            I += float(mutValues[1])+float(mutValues[10])
            D += float(mutValues[2])+float(mutValues[11])
            S += float(mutValues[3])+float(mutValues[12])
            H += float(mutValues[4])+float(mutValues[13])
            N += float(mutValues[5])+float(mutValues[14])
            P += float(mutValues[6])+float(mutValues[15])
            X += float(mutValues[7])+float(mutValues[16])
            Egal += float(mutValues[8])+float(mutValues[17])
        FinalCigar.write("Global cigar mutation observed :"+"\n" #Write the results to the corresponding output table
                        +"Alignment Match : "+str(round(M/nbReads,4))+"% "+str(round(M,4))+" total"+"\n" #Calculate all percentages
                        +"Insertion : "+str(round(I/nbReads,4))+"% "+str(round(I,4))+" total"+"\n"
                        +"Deletion : "+str(round(D/nbReads,4))+"% "+str(round(D,4))+" total"+"\n"
                        +"Skipped region : "+str(round(S/nbReads,4))+"% "+str(round(S,4))+" total"+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,4))+"% "+str(round(H,4))+" total"+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,4))+"% "+str(round(N,4))+" total"+"\n"
                        +"Padding : "+str(round(P/nbReads,4))+"% "+str(round(P,4))+" total"+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,4))+"% "+str(round(Egal,4))+" total"+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,4))+"% "+str(round(X,4))+" total"+"\n")
        
def Ids(file,minval): #Function that stores the read IDs and their number of occurrences
    D={}
    with open(file,"r") as Str:  
        L=list(islice(enumerate(Str),2,None))
    List=[t[1] for t in L if int(t[1].split("\t")[4]) > minval]
    for t in List:
        t=t.split("\t")
        if t[0] in D.keys():
            D[t[0]]=2
        else:
            D[t[0]]=1
    return D    
    
def SamRead2(file,minval): # Function that stores reads in a dictionary with Key = ID and Value = the rest of the line
    D={}
    with open(file,"r") as Str:  
        L=list(islice(enumerate(Str),2,None))
    List=[t[1] for t in L if int(t[1].split("\t")[4]) > minval]
    for t in List:
        t=t.split("\t")
        if t[0] not in D.keys():
            D[t[0]]=t[1:]
        else :
            D[t[0]+"-1"]=t[1:]
    return D

def unmapped2(sam_line):
    unmapped_count = 0
    with open ("../Results/only_unmapped.fasta", "a+") as unmapped_fasta, open("../Results/summary_unmapped.txt", "w") as summary_file:
        for keys,line in sam_line.items():
            flag = flagBinary(line[0])
            if int(flag[-3]) == 1:
                unmapped_count += 1
                unmapped_fasta.write(keys+"\t"+str(line)[1:-1]+"\n") #retired "toStringOutput"
        summary_file.write("Total unmapped reads: " + str(unmapped_count)) 
        return unmapped_count

def partiallyMapped2(sam_line):
    partially_mapped_count = 0
    with open ("../Results/only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("../Results/summary_partially_mapped.txt", "w") as summary_file:
        for keys,line in sam_line.items():
            flag = flagBinary(line[0]) # We compute the same 
            if int(flag[-2]) == 1: 
                if line[4] != "100M":
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(keys+"\t"+str(line)[1:-1]+"\n") #retired "toStringOutput"
        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count
    
def Store2(file):
    if check_file(file):
        sam_lines=SamRead2(file)
        unmapped_count=unmapped2(sam_lines)
        partially_mapped_count=partiallyMapped2(sam_lines)
        return unmapped_count,partially_mapped_count
    else:
        return None 

#### Summarise the results ####

#def Summary(fileName):
    
   

#### Main function ####

def main(argv):
    if len(argv) < 2:
        print("Error, an argument is missing") #return error if an argument is missing ( 2 arguments required so len(argv)=2 )
        return None
    minval=int(argv[1]) #the second argument is the minimum quality value required
    input_file=argv[0] #the first argument is the path of the SAM file
    #check_file(input_file)
    unmapped_count,partially_mapped_count,Pm,un,par,idd,Qual=Store(input_file,minval)
    T,TT=Analyse(Pm,un,par,idd,input_file,minval)
    globalPercentCigar()
    

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])

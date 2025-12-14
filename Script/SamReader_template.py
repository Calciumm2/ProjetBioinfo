#!/usr/bin/python3
#-*- coding : utf-8 -*-
#CHANGEMENT POUR TESTER COMMIT


__authors__ = ("Mathilde Chatain", "Lucien Maurau")
__contact__ = ("mathilde.chatain@etu.umontpellier.fr","lucien.maurau@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "12/14/2021"
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


#Fonction qui vérifie si le lien vers "fichier" existe
def check_file(file):
    if not os.path.exists(file):
        print(f"Erreur : le chemin '{file}' n'existe pas.")
        return False

    if not os.path.isfile(file):
        print(f"Erreur : le chemin '{file}' n'est pas un fichier.")
        return False
    taille = os.path.getsize(file)
    if taille == 0:
        print(f"Erreur : le fichier '{file}' est vide.")
        return False

    print(f"OK : le fichier '{file}' existe et est non vide (taille : {taille} octets).")
    return True

## 2/ Read, 
#Fonction qui lis le fichier, sépare ligne 
def SamRead(file,minval):
    with open(file,"r") as Str:  
        L=list(islice(enumerate(Str),2,None))
    return [t[1] for t in L if int(t[1].split("\t")[4]) > minval] #On prend [1] pour ne pas avoir le num de la ligne
    
## 3/ Store,
#Fonction qui appelle check_file et samread et qui stocke les résultats de unmapped et de partiallymapped
def Store(file,minval):
    if check_file(file): #On check si le fichier est non vide
        Qual=ReadQuality(SamRead2(file, -1))
        sam_lines=SamRead(file,minval) #On obtient une liste des lignes
        unmapped_count,un=unmapped(sam_lines) #On obtient le nombre de read unmapped et un dico les contenant
        partially_mapped_count,par=partiallyMapped(sam_lines) #On obtient le nombre de read partiallymapped et un dico les contenant
        Pm=PerfectMapped(sam_lines) #On obtient un dico contenant les Id de read perfectMapped (on ne garde qu'une seule occurence pour les paires car on ne s'intéresse pas au contenu de la ligne)
        idd=Ids(file,minval) #On obtient un dico contenant les Id de chaque read, avec son nombre d'occurence (on veut être sûr que chaque read est paired)
        return unmapped_count,partially_mapped_count,Pm,un,par,idd,Qual
    else:
        return None #Changer pour un exit 
    
def Analyse(Pm,Un,Par,idd,file,minval): #Fonction qui va analyser les paires de read, et le pourcentage global de chaque mutation grace au cigar
    MapUnMap={} #Dico qui contiendra les paires de read Mapped-UnMapped
    MapParMap={} #Dico qui contiendra les paires de read Mapped-PartMapped
    IdLines=SamRead2(file,minval) #Dico qui contient toutes les lignes selon clé = Id du read, valeur = reste de la ligne
    for s in idd.keys(): #On parcourt pour chaque Id de read 
        if s in Pm and s in Un: #On regarde si la paire est présente dans perfectMapped et dans Unmapped
            MapUnMap[s]=Un[s] #On ajoute le premier read de la paire 
            MapUnMap[s+"-1"]=Pm[s] #/// le second
        if s in Pm and s in Par:  #On regarde si la paire est présente dans parMapped et dans Unmapped
            MapParMap[s]=Par[s] #Même que avant
            MapParMap[s+"-1"]=Pm[s]
        with open("../Results/outpuTable_cigar.txt","a+") as outpuTable: #On va écrire les pourcentages de mutations pour chaque ligne
            if idd[s]==2: #On vérifie que l'Id est associé à une paire de read
                D1,D2=readCigar(IdLines[s][4]),readCigar(IdLines[s+"-1"][4]) #On obtient les deux % de chaque mutations pour les deux read de la paire 
                LTXT=percentMutation(D1)+";"+percentMutation(D2) #On concatène en une ligne 
                outpuTable.write(LTXT+"\n") #On écrit dans un fichier Txt qui contiendra les résultats 
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


#Prend les lignes unmapped et écrit un fichier qui les contient en Fasta et un fichier résumé qui donne le nombre de reads unmapped
def unmapped(sam_line):
    D={} #Dico qui contiendra les read UnMapped
    unmapped_count = 0
    with open ("../Results/only_unmapped.fasta", "a+") as unmapped_fasta, open("../Results/summary_unmapped.txt", "w") as summary_file:
        for line in sam_line: 
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])
            if int(flag[-3]) == 1:
                D[col_line[0]]=col_line[1:]
                unmapped_count += 1
                unmapped_fasta.write(line) #retiré "toStringOutput"
        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
        return unmapped_count,D

#Fonction qui compte les read partially mapped en regardant le CIGAR
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
                    partillay_mapped_fasta.write(line) #retiré "toStringOutput"
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
    ext = re.findall('\w',cigar) #sépare tous les caractères
    key=[] 
    value=[]    
    val=""
    for i in range(0,len(ext)): # For each numeric values or alpha numeric
        if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
            key.append(ext[i]) #if précédent : si on trouve une lettre, on ajoute au dico, clé : la lettre, valeur : le numéro avant 
            value.append(val) #on ajoute la valeur
            val = "" #on rénétialise
        else :
            val = "" + val + ext[i]  #On stocke valeur en parcourant le cigar
    dico = {}
    n = 0
    for k in key:   #On construit le dico             
        if k not in dico.keys():    #On ne veut pas deux fois la même clé, on vérifie à chaque si elle est déjà dedans
            dico[k] = int(value[n])   #On assigne la valeur
            n += 1
        else:
            dico[k] += int(value[n])  #si la clé existe déjà, on ajoute
            n += 1
    return dico #On obtient un dico qui donne la somme des opérations, mais plus leur ordre (qui était consigné dans le CIGAR)

#Fonction qui analyse le dico crée par la fonction readCigar
def percentMutation(dico):
    totalValue = 0 #valeur totale mutations
    for v in dico :
        totalValue += dico[v] #nombre total d'opérations consignés dans le cigar
    mutList = ['M','I','D','S','H','N','P','X','='] #Bizarre, pourquoi il y a le M qui devrait correspondre à "match", c'est pas une mutation ?? 
    res = "" 
    for mut in mutList : # Pour chaque de mutation on itère dans le dico
        if mut in dico.keys() :
            res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";") #On ajoute le pourcentage de l'opération considérée
        else :
            res += ("0.00" + ";") 
    return res[:-1] #on met zéro sinon, l'ordre spécifie l'opération (voir fonction suivante), donc pas besoin de nommer à chaque fois

#Fonction qui écrit un .txt qui résume la distribution de tous les cigar à partir d'un cigar contenant les contenus de la fonction précédente (itérée paire par paire ?)
def globalPercentCigar():
    with open ("../Results/outpuTable_cigar.txt","r") as outpuTable, open("../Results/Final_Cigar_table.txt", "w") as FinalCigar: #On explicite le .txt qu'on lit et écrit
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)] #chaque valeur initiée à zéro
        for line in outpuTable : #on parcourt chaque ligne de la table donnée
            mutValues = line.split(";")
            nbReads += 2 #nombre de reads total (2 par ligne)
            M += float(mutValues[0])+float(mutValues[9])
            I += float(mutValues[1])+float(mutValues[10])
            D += float(mutValues[2])+float(mutValues[11])
            S += float(mutValues[3])+float(mutValues[12])
            H += float(mutValues[4])+float(mutValues[13])
            N += float(mutValues[5])+float(mutValues[14])
            P += float(mutValues[6])+float(mutValues[15])
            X += float(mutValues[7])+float(mutValues[16])
            Egal += float(mutValues[8])+float(mutValues[17])
        FinalCigar.write("Global cigar mutation observed :"+"\n" #On écrit les résultats dans la table correspondante 
                        +"Alignment Match : "+str(round(M/nbReads,4))+"% "+str(round(M,4))+" total"+"\n" #On calcule tous les pourcentages
                        +"Insertion : "+str(round(I/nbReads,4))+"% "+str(round(I,4))+" total"+"\n"
                        +"Deletion : "+str(round(D/nbReads,4))+"% "+str(round(D,4))+" total"+"\n"
                        +"Skipped region : "+str(round(S/nbReads,4))+"% "+str(round(S,4))+" total"+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,4))+"% "+str(round(H,4))+" total"+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,4))+"% "+str(round(N,4))+" total"+"\n"
                        +"Padding : "+str(round(P/nbReads,4))+"% "+str(round(P,4))+" total"+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,4))+"% "+str(round(Egal,4))+" total"+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,4))+"% "+str(round(X,4))+" total"+"\n")
        
def Ids(file,minval): #Fonction qui store les Id des read et leur nombre d'occurence 
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
    
def SamRead2(file,minval): #Fonction qui store les read dans un dico sous la forme Clé = Id, Valeur = reste de la ligne 
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
                unmapped_fasta.write(keys+"\t"+str(line)[1:-1]+"\n") #retiré "toStringOutput"
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
                    partillay_mapped_fasta.write(keys+"\t"+str(line)[1:-1]+"\n") #retiré "toStringOutput"
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
    if len(argv) < 1:
        print("erreur")
        return None
    minval=int(argv[1])
    input_file=argv[0]
    #check_file(input_file)
    unmapped_count,partially_mapped_count,Pm,un,par,idd,Qual=Store(input_file,minval)
    T,TT=Analyse(Pm,un,par,idd,input_file,minval)
    globalPercentCigar()
    

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main(sys.argv[1:])

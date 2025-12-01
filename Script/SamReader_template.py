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
def SamRead(file):
    with open(file,"r") as Str:  
        L=list(islice(enumerate(Str),2,None))
    return [t[1] for t in L]
    
## 3/ Store,
#Fonction qui appelle check_file et samread et qui stocke les résultats de unmapped et de partiallymapped
def Store(file):
    if check_file(file):
        sam_lines=SamRead(file)
        unmapped_count=unmapped(sam_lines)
        partially_mapped_count=partiallyMapped(sam_lines)
        return unmapped_count,partially_mapped_count,
    else:
        return None 

## 4/ Analyse 

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
    unmapped_count = 0
    with open ("../Results/only_unmapped.fasta", "a+") as unmapped_fasta, open("../Results/summary_unmapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])
            if int(flag[-3]) == 1:
                unmapped_count += 1
                unmapped_fasta.write(line) #retiré "toStringOutput"
        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
        return unmapped_count

#Fonction qui compte les read partially mapped en regardant le CIGAR
def partiallyMapped(sam_line):
    partially_mapped_count = 0
    with open ("../Results/only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("../Results/summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1]) # We compute the same 
            if int(flag[-2]) == 1: 
                if col_line[5] != "100M":
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(line) #retiré "toStringOutput"
        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count


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
    with open ("outpuTable_cigar.txt","r") as outpuTable, open("Final_Cigar_table.txt", "w") as FinalCigar: #On explicite le .txt qu'on lit et écrit
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)] #chaque valeur initiée à zéro
        for line in outpuTable : #on parcourt chaque ligne de la table donnée
            mutValues = line.split(";")
            nbReads += 2 #nombre de reads par lignes ?
            M += float(mutValues[2])+float(mutValues[12]) #théorie : la valeur "de gauche" (resp "de droite") vaut pour le M du premier (resp deuxième) read de la paire 
            I += float(mutValues[3])+float(mutValues[13])
            D += float(mutValues[4])+float(mutValues[14])
            S += float(mutValues[5])+float(mutValues[15])
            H += float(mutValues[6])+float(mutValues[16])
            N += float(mutValues[7])+float(mutValues[17])
            P += float(mutValues[8])+float(mutValues[18])
            X += float(mutValues[9])+float(mutValues[19])
            Egal += float(mutValues[10])+float(mutValues[20])
        FinalCigar.write("Global cigar mutation observed :"+"\n" #On écrit les résultats dans la table correspondante 
                        +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
                        +"Insertion : "+str(round(I/nbReads,2))+"\n"
                        +"Deletion : "+str(round(D/nbReads,2))+"\n"
                        +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                        +"Padding : "+str(round(P/nbReads,2))+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")
        
        
def test(line):
    nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)] 
    print(M)
    mutValues = line.split(";")
    nbReads += 2
    M += float(mutValues[0])#+float(mutValues[12])
    I += float(mutValues[1])#+float(mutValues[13])
    D += float(mutValues[2])#+float(mutValues[14])
    S += float(mutValues[3])#+float(mutValues[15])
    H += float(mutValues[4])#+float(mutValues[16])
    N += float(mutValues[5])#+float(mutValues[17])
    P += float(mutValues[6])#+float(mutValues[18])
    X += float(mutValues[7])#+float(mutValues[19])
    Egal += float(mutValues[8])#+float(mutValues[20])
    print("Global cigar mutation observed :"+"\n"
                +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
                +"Insertion : "+str(round(I/nbReads,2))+"\n"
                +"Deletion : "+str(round(D/nbReads,2))+"\n"
                +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                +"Padding : "+str(round(P/nbReads,2))+"\n"
                +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")
    
def SamRead2(file):
    D={}
    with open(file,"r") as Str:  
        L=list(islice(enumerate(Str),2,None))
    List=[t[1] for t in L]
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
        for line in sam_line.values():
            flag = flagBinary(line[0])
            if int(flag[-3]) == 1:
                unmapped_count += 1
                unmapped_fasta.write(str(line)[1:-1]) #retiré "toStringOutput"
        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
        return unmapped_count

def partiallyMapped2(sam_line):
    partially_mapped_count = 0
    with open ("../Results/only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("../Results/summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_line.values():
            flag = flagBinary(line[0]) # We compute the same 
            if int(flag[-2]) == 1: 
                if line[4] != "100M":
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(str(line)[1:-1]) #retiré "toStringOutput"
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

def toto(toto):
    for t in toto.values():
        if len(t[0]) != 3:
            print(t[0])

#### Summarise the results ####

#def Summary(fileName):
    
   

#### Main function ####

#def main(argv):
    

############### LAUNCH THE SCRIPT ###############

#if __name__ == "__main__":
#    main(sys.argv[1:])

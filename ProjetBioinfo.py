__authors__ = ("Mathilde Chatain", "XXX")
__contact__ = ("mathilde.chatain@etu.umontpellier.fr","XXX@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "18.11.2025"
__licence__ ="jjj."

## 1/ Check, 
#si repertoire est fichier, si fichier existe et si il est non vide
#verifier la ligne de commande

import os;
import sys;

def check_file(fichier):
    if not os.path.exists(fichier):
        print(f"Erreur : le chemin '{fichier}' n'existe pas.")
        return False

    if not os.path.isfile(fichier):
        print(f"Erreur : le chemin '{fichier}' n'est pas un fichier.")
        return False
    taille = os.path.getsize(fichier)
    if taille == 0:
        print(f"Erreur : le fichier '{fichier}' est vide.")
        return False

    print(f"OK : le fichier '{fichier}' existe et est non vide (taille : {taille} octets).")
    return True




## 2/ Read, 
#pouvoir lire fichier (readlines en une seule fois, readline en lisant ligne par ligne) puis split pour recuperer les infos qui nous interesse


## 3/ Store, stockage, reflechir structure donnee ( liste, dictionnaire, dictionnaire de liste, dictionnaire de dictionnaire de liste), quelle est la cle la plus pertinente?

## 4/ Analyse avec les quatres etapes - nombre reads mappés / combien ? -comment ces reads sont mappes ( nb de read pour chaque flag(comment la paire est alignee)) et CIGAR(comment est aligne le read)

def main(argv):
    check_file(argv)

main(sys.argv[1])
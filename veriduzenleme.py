import pandas as pd 
def veriyial():
    aminoasitler = ["H","R","K","I","F","L","W","A","M","P","C","N","V","G","S","Q","Y","D","E","T"]
    polarite = {aminoasitler[i]:[0,0,0,1,1,1,1,1,1,1,2,2,1,2,2,2,2,3,3,2][i] for i in range(len(aminoasitler))}
    hacim = {aminoasitler[i]:[153.2,173.4,168.6,166.7,189.9,166.7,227.8,88.6,162.9,112.7,108.5,114.1,140.0,60.1,89.0,143.8,193.6,111.1,138.4,116.1][i] for i in range(len(aminoasitler))}
    features = pd.read_csv("https://raw.githubusercontent.com/BiyoinformatikProje/Biyoinformatik-Projesi/main/features.csv")
    mutasyonlar = pd.read_csv("https://raw.githubusercontent.com/BiyoinformatikProje/Biyoinformatik-Projesi/main/RBD_all_mutations.csv")
    mkütle = {features["Amino acid code"][i]:features["Molecular Mass (Da)"][i] for i in range(len(features))}
    hidropati = {features["Amino acid code"][i]:features["Hydropathy index"][i] for i in range(len(features))}
    hidropatifarkı = []
    mkütlefarkı = []
    hacimfarkı = []
    polaritefarkı = []
    mutantlar = []
    wildtypelar = []
    for i in range(len(mutasyonlar)):
        wildtype = mutasyonlar["mutation"][i][0]
        mutant = mutasyonlar["mutation"][i][len(mutasyonlar["mutation"][i])-1]
        mutantlar.append(mutant)
        wildtypelar.append(wildtype)
        hidropatifarkı.append(hidropati[mutant] - hidropati[wildtype])
        mkütlefarkı.append(mkütle[mutant] - mkütle[wildtype])
        hacimfarkı.append(hacim[mutant] - hacim[wildtype])
        polaritefarkı.append(polarite[mutant] - polarite[wildtype])
    data = pd.DataFrame(data = {"wildtype":wildtypelar,"mutant":mutantlar,"hydropathy m-w":hidropatifarkı,"molecular mass m-w":mkütlefarkı,"polarity m-w":polaritefarkı,"volume m-w":hacimfarkı,"bind_avg":mutasyonlar["bind_avg"],"expr_avg":mutasyonlar["expr_avg"]})
    return data
def yeniozellikekle(data):
    A = "Positively Charged Side Chains"
    B = "Negatively Charged Side Chains"
    C = "Polar Uncharged Side Chains"
    D = "Special Case"
    E = "Hydrophobic Side Chain"
    aminoacids = ["R","H","K","D","E","S","T","N","Q","C","U","G","P","A","V","I","L","M","F","Y","W"]
    ring = [0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,2]
    dbond = [2,1,1,2,2,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1]
    state = [A,A,A,B,B,C,C,C,C,D,D,D,D,E,E,E,E,E,E,E,E]
    carbon   = [6,3,6,4,5,3,4,4,5,3,3,2,5,3,5,6,6,5,3,3,3]
    oxygen   = [2,2,2,4,4,3,2,3,3,2,2,2,2,2,2,2,2,2,2,2,2]
    hydrogen = [14,7,14,7,9,7,9,8,10,7,7,5,9,7,11,13,13,11,11,11,12]
    ring_d     = {}
    dbond_d    = {}
    state_d    = {}
    oxygen_d   = {}
    hydrogen_d = {}
    for i in range(len(aminoacids)):
        acid = aminoacids[i]
        dbond_d[acid] = dbond[i]
        state_d[acid] = state[i]
        oxygen_d[acid]= oxygen[i]
        hydrogen_d[acid] = hydrogen[i]
        ring_d[acid] = ring[i]
    dbonddiff = [dbond_d[data.iloc[i]["mutant"]] - dbond_d[data.iloc[i]["wildtype"]] for i in range(len(data))]
    oxygendiff = [oxygen_d[data.iloc[i]["mutant"]] - oxygen_d[data.iloc[i]["wildtype"]] for i in range(len(data))]
    hydrogendiff = [hydrogen_d[data.iloc[i]["mutant"]] - hydrogen_d[data.iloc[i]["wildtype"]] for i in range(len(data))]
    ringdiff = [ring_d[data.iloc[i]["mutant"]] - ring_d[data.iloc[i]["wildtype"]] for i in range(len(data))]
    df = pd.DataFrame(data = {"Double Bond m-w":dbonddiff,"Oxygen m-w":oxygendiff,"Hydrogen m-w":hydrogendiff,"Ring m-w":ringdiff})
    data["Double Bond m-w"] = dbonddiff 
    data["Oxygen m-w"] = oxygendiff 
    data["Hydrogen m-w"] = hydrogendiff 
    data["Ring m-w"] = ringdiff
    return data

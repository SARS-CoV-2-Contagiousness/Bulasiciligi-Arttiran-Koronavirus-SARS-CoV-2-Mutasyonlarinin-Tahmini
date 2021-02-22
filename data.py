import pandas as pd 
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
data.to_csv("data.csv",index = False)
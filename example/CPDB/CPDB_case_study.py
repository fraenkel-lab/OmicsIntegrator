import os

def conf_prep(mu,beta,D,w):
    file = open("conf.txt","w")
    file.writelines("w = %d\nb = %d\nD = %d\nmu = %f" % (w,beta,D,mu))
    file.close()

# Update this with the path to msgsteiner
msgpath = "./msgsteiner"

mu_range = [0,0.009,0.015]
beta_range = [1,5,10]
w_range = [1,2,3]
prize_file = "pyr_rho_mrna_noUBC.terminal"
edge_file = "../../data/iref_mitab_miscore_2013_08_12_interactome.txt"
conf_file = "conf.txt"
output_path = "../../results"
D = 5

# Create output directory if needed
if not os.path.exists(output_path):
    os.makedirs(output_path)

for mu in mu_range:
    for beta in beta_range:
        for w in w_range:
            conf_prep(mu,beta,D,w)
            out_label = "CPDB_w%f_beta%d_D%d_mu%f" %(w,beta,D,mu)
            os.system("python ../../scripts/forest.py --prize %s --edge %s --conf conf.txt --msgpath %s --outpath %s --outlabel %s" %(prize_file,edge_file,msgpath,output_path,out_label))

import os

def frange(start,end,step):
    return map(lambda x: x*step, range(int(start*1./step),int(end*1./step)))

def conf_prep(mu,beta,D,w):
    file = open("conf.txt","w")
    file.writelines("w = %d\nb = %d\nD = %d\nmu = %f" % (w,beta,D,mu))
    file.close()

mu_range = frange(0.002,0.004,0.002)
beta_range = frange(150,160,10)
w_range = frange(2,3,1)
prize_file = "gbm_prize.txt"
edge_file = "../../data/iref_mitab_miscore_2013_08_12_interactome.txt"
conf_file = "conf.txt"
output_path = "../../results/WT/"

D = 10

for mu in mu_range:
    for beta in beta_range:
        for w in w_range:
            conf_prep(mu,beta,D,w)
            output_path = "../../results/WT/"
            out_label = "WT_w%f_beta%d_D%d_mu%f" %(w,beta,D,mu)
            os.system("python ../../scripts/forest.py --prize %s --edge %s --conf conf.txt --msgpath /nfs/apps/bin/msgsteiner9 --outpath %s --outlabel %s" %(prize_file,edge_file,output_path,out_label))
            output_path = "../../results/KO/"
            out_label = "KO_w%f_beta%d_D%d_mu%f" %(w,beta,D,mu)
            os.system("python ../../scripts/forest.py --prize %s --edge %s --conf conf.txt --msgpath /nfs/apps/bin/msgsteiner9 --knockout EGFR --outpath %s --outlabel %s" %(prize_file,edge_file,output_path,out_label))

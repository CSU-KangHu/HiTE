import os ,sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def slice(bedfile,name, score_file):

    if not os.path.exists(name):
        os.mkdir(name)

    coords = []

    with open(bedfile,'r') as infile:

        lines = infile.readlines()

        for l in lines:
            fields = l.split()
            info = (fields[1],fields[2],fields[4],fields[5])
            coords.append(info)
    
    scores =[]

    with open(score_file,'r') as infile:
        scores = infile.readlines()
        

    for el in coords:

        print(el)
        start = int(el[0])
        end = int(el[1])

        

        mid1 = int(el[2]) 
        mid2 = int(el[3]) 
        

        bounds =[start,mid1,mid2,end-1]
        print(bounds)
        
        slice =name+"/"+el[0]+"-"+el[1]+".png"

        curr = scores[start:end]

        x = [int(l.split(',')[0]) for l in curr]
        y = [int(l.split(',')[1]) for l in curr]

        

        d = {"idx":x, "score":y}
        

        df = pd.DataFrame(d,index = d["idx"])

        graph = df.plot(x = "idx", y = "score")

        #f,ax = plt.subplots(1)

        graph.plot(bounds,[0.0]*len(bounds),"rv")

        fig = graph.get_figure()
        #plt.plot( bounds, marker = "v", markerfacecolor ="r")
        fig.savefig(slice)


if __name__ == "__main__":

    chrom = sys.argv[1]
    output_dir = sys.argv[2]
    score_file = sys.argv[3]



    slice(chrom,output_dir,score_file)

        
        



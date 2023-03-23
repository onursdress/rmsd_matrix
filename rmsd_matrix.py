from pymol import cmd, stored
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt

# Credits to PyMOL wiki for this function
def pick_Ca(userSelection):
    # this array will be used to hold the coordinates.  It
    # has access to PyMOL objects and, we have access to it.
    stored.alphaCarbons = []

    # let's just get the alpha carbons, so make the
    # selection just for them
    userSelection = userSelection + " & n. ca"
    return userSelection

# function calculates the rmsd of two structure objects with 'cmd.align()'
# 'align' is interchangable with other alignment algorithms, e.g. 'super'
def pairRMSD(seq1, seq2):

    # picking all the Ca atoms from mobile/target object 
    mobile = pick_Ca(seq1)
    target = pick_Ca(seq2)
    
    # temporary storage of one alignment
    aln = cmd.align(mobile, target)
    
    # this function returns the first output from cmd.align() which is the rmsd
    return aln[0]
    
# this function is called by PyMOL command
# this function stores all pairRMSD() with object names and rmsd or
# only rmsd
# directs the rmsd matrix to the function heatmap_rmsd()
def rmsd_matrix():
    # list with object names and rmsds 
    rmsd_o1_o2 = []
    # list with rmsds 
    rmsd_matrix = []

    allobjs = cmd.get_object_list('(all)')
    for o1 in allobjs:
        # temporary lists
        rmsd_all = []
        rmsd_only = []

        for o2 in allobjs:
            # if the alignment is between the same object, insert a '0' 
            if o1==o2:
                rmsd_all.extend([0, o1, o2])
                rmsd_only.extend([0])
            else: 
                rmsd_all.extend([pairRMSD(o1, o2), o1, o2])
                rmsd_only.extend([pairRMSD(o1, o2)])
        
        rmsd_o1_o2.append(rmsd_all)
        rmsd_matrix.append(rmsd_only)
            

    arr=np.array(rmsd_matrix)
    minimum=arr.min()
    maximum=arr.max()
    print(arr)
    
    # prints in the same directory the rmsd matrix as an excel file [optional]
    #df = pd.DataFrame(arr, index=allobjs, columns=allobjs)
    #print(df)
    #df.to_excel(***insert directory***)
    
    heatmap_rmsd(arr, allobjs, minimum, maximum)

# function to plot the rmsd matrix as a heatmap
def heatmap_rmsd(array, labels, minimum, maximum):
    fig, ax = plt.subplots()
    heat = ax.imshow(array, cmap="gist_earth", vmin=minimum, vmax=maximum)
    
    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_yticks(np.arange(len(labels)), labels=labels)
       
    plt.setp(ax.get_xticklabels(), rotation=90)
    max_round = round(maximum)
    plt.colorbar(heat, boundaries=[0, max_round*.125, max_round*.25,max_round*.375, max_round*.5,max_round*.625, max_round*.75, max_round*.875, max_round], drawedges=True)
    
    #change the title if you want, depending on the alignment algorithm
    ax.set_title("'align'-algorithm generated\nroot-mean-square deviations of C" + chr(945) + " atoms")
    fig.tight_layout()
    plt.show()
  
cmd.extend("rmsd_matrix", rmsd_matrix)
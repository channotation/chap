import numpy as np
import pandas
import seaborn as sns
import math
import json
import sys
import os
import matplotlib.pyplot as plt

def distance_orthogonal_projection(point, c_x, c_y):
    """
    Calculate the orthogonal projection distance from a point to a straight line defined by a contour.

    Parameters:
    - point (numpy.ndarray): The coordinates of the point (x, y) to calculate the distance from.
    - c_x (numpy.ndarray): The x-coordinates of the contour line defining the straight line.
    - c_y (numpy.ndarray): The y-coordinates of the contour line defining the straight line.

    Returns:
    float: The orthogonal projection distance from the point to the straight line defined by the contour.

    The function uses the orthogonal projection formula to find the closest point on the contour line
    to the given point and calculates the distance between them.

    Reference:
    - https://de.wikipedia.org/wiki/Orthogonalprojektion

    Example:
    >>> point = np.array([3, 4])
    >>> c_x = np.array([1, 2, 3, 4, 5])
    >>> c_y = np.array([0, 1, 3, 1, 0])
    >>> distance = distance_orthogonal_projection(point, c_x, c_y)
    >>> print(distance)
    0.4472135954999579
    """
    ### first and last point to define straight line ###
    ind1 = 0
    ind2 = -1
    ### use closest point on contour line with respect to x ###
    ind1 = np.argmin(abs(c_x-point[0]))
    ind2 = ind1 - 1
    #print('c_x', c_x[ind1], point[0])
    
    r0 = np.array([c_x[ind1], c_y[ind1]])
    u = np.array([c_x[ind2]-c_x[ind1], c_y[ind2]-c_y[ind1]])
    lambd = np.dot((point-r0), u)/np.dot(u,u)
    projection = r0 + lambd * u
    dist = np.linalg.norm(point - projection)
    return dist

def heuristic_score(grid, chap_output, path, title='GlyR', fonts=12, fac=1, f_size=14):
    """
    Generate a heuristic prediction score for pore-lining residues in a biological channel.
    Based on Rao et al., 2019.
    https://doi.org/10.1073%2Fpnas.1902702116

    Parameters:
    - grid (str): Path to a JSON file containing grid data with hydrophobicity, radius, and energy information.
        The grid can be found at: https://github.com/channotation/chap/blob/master/scripts/heuristic/heuristic_grid.json
    - chap_output (str): Path to a JSON file containing output data from CHAP.
    - path (str): Path to the directory where the resulting plot image will be saved.
    - title (str, optional): Title for the plot. Default is 'GlyR'.
    - fonts (int, optional): Font size for the plot. Default is 12.
    - fac (int, optional): Scaling factor for the channel coordinate. Default is 1.

    Returns:
    None

    This function reads grid and channel output data from JSON files, performs a heuristic prediction for structure,
    and generates a plot showing hydrophobicity-radius scatter and contour lines. It also calculates a heuristic score
    based on the orthogonal projection distance from pore-lining residues to the contour lines.

    The resulting plot is saved as a PNG file in the specified path.

    Example:
    >>> grid_file = "grid_data.json"
    >>> chap_output_file = "channel_output.json"
    >>> result_path = "./results/"
    >>> heuristic_score(grid_file, chap_output_file, result_path, title='MyChannel', fonts=14, fac=0.5)
    ### example ###
    grid = 'heuristic_grid.json'
    chap_output = 'Glycine.json'
    path = '/biggin/b198/orie4254/Documents/asymmetrical_Pores/GlyR/GlyR_heuristic_analysis/'

    heuristic_score(grid=path+grid, chap_output=path+chap_output, path=path, 
                    title='GlyR Glycine model (pdb)')
    """
    f_grid = open(grid,'r')
    data_grid = json.load(f_grid)
    #f_grid.close()

    f_chap = open(chap_output,'r')
    data_chap = json.load(f_chap)
    #f_chap.close()
    
    ### reformat grid ###
    hydrophobicity = []
    radius = []
    energy = []
    for i in range(len(data_grid)):
        hydrophobicity.append(data_grid[i]['hydrophobicity'])
        radius.append(data_grid[i]['radius'])
        energy.append(data_grid[i]['energy'])
    d ={'hydrophobicity': hydrophobicity, 'radius': radius, 'energy': energy}
    df_grid = pandas.DataFrame(data=d)
    
    ### Prepare data frame for pore-lining residues ###
    pore_facing = np.array(data_chap["residueSummary"]["poreFacing"]["mean"]) > 0.5
    narrow = np.array(data_chap['residueSummary']['poreRadius']['mean'] ) < 0.7
    print('pore_facing',np.unique(pore_facing), 'narrow', np.unique(narrow))

    name = []
    hydrophobicity = []
    radius = []
    channel_coord = []
    for residue in np.array(data_chap['residueSummary']['id'])[narrow * pore_facing]:
        nearest = np.argmin(abs(
            data_chap['pathwayProfile']['s']
             - np.array(data_chap['residueSummary']['s']['mean'])[ 
                 np.array(data_chap['residueSummary']['id']) == residue]))
        n = np.array(data_chap['residueSummary']['name'])[
            np.array(data_chap['residueSummary']['id']) == residue][0]
        h = np.array(data_chap['pathwayProfile']['pfHydrophobicityMean'])[nearest]
        r = np.array(data_chap['pathwayProfile']['radiusMean'])[nearest]
        s = np.array(data_chap['pathwayProfile']['s'])[nearest]
        #print(residue, nearest, n, h, r)
        name.append(n)
        hydrophobicity.append(h)
        radius.append(r)
        channel_coord.append(s)
    d ={'name': name, 'hydrophobicity': hydrophobicity, 
        'radius': radius, 'channel_coord': channel_coord}
    df = pandas.DataFrame(data=d)
    
    #### Heuristic prediction for structure ####
    pred = []
    for i in range(len(df)):
        #np.array(data_chap['
        p = np.array(df_grid['energy'])[
            np.argmin(
                pow(np.array(df_grid['hydrophobicity']) - df['hydrophobicity'].iloc[i], 2)
                + pow(np.array(df_grid['radius']) - df['radius'].iloc[i], 2)
            )
        ]
        #print(p)
        pred.append(p)
    df['pred'] = pred 
    
    #### Plot ###
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(7, 4),
                            gridspec_kw={'width_ratios': [2, 1]})
    axs[0].scatter(df['hydrophobicity'], df['radius'], marker='x')
    indices = df['pred'] > 2.6
    axs[0].scatter(df['hydrophobicity'][indices], df['radius'][indices], marker='x',  color='red')

    ### contour plot ###
    x = np.array(df_grid['hydrophobicity'])
    y = np.array(df_grid['radius'] )
    z = np.array(df_grid['energy'] )
    print(len(x), len(y),len(z))
    num_x = np.size(np.unique(x))
    num_y = np.size(np.unique(y))
    print('num_x, num_y', num_x, num_x)
    X = x.reshape(num_x, num_y)
    Y = y.reshape(num_x, num_y)
    Z = z.reshape(num_x, num_y)
    CS = axs[0].contour(X, Y, Z, 
                colors='black',
                levels = [2.6]           )
    #plt.clabel(CS, colors = 'k', fontsize=14)

    ### heuristic score ###
    n = np.array( df['name'][indices])
    hydro = np.array(df['hydrophobicity'][indices])
    rad = np.array( df['radius'][indices])
    s = np.array(df['channel_coord'][indices])
    N = len(rad)
    c_x, c_y = CS.collections[0].get_paths()[0].vertices.T
    #plt.scatter(c_x, c_y, marker='x')
    #plt.plot([c_x[0], c_x[-1]], [c_y[0], c_y[-1]], '-o', color='red')
    score = 0
    for i in range(N):
        dist = distance_orthogonal_projection(point=np.array([hydro[i], rad[i]]), 
                                        c_x=c_x, c_y=c_y)
        print(n[i],hydro[i], rad[i], s[i] , dist)
        score = score + dist
    print('score', score)

    #for i, txt in enumerate(n):
    #    ax.annotate(txt, (hydro[i], rad[i]))
    xlim = axs[0].get_xlim()
    ylim = axs[0].get_ylim()
    #plt.subplots_adjust(right=0.35)
    string = title+'\nHeuristic prediction score = ' +str(round(score,2)) + '\n\n'
    if len(n)!=0:
        string = string + 'Residues with energies\nlarger than 1 RT:\n\n'
    for i, txt in enumerate(n):
        if len(n)<11:
            string = string + txt + ' r={:.2f}'.format(rad[i]) +'  s={:.2f}'.format(fac*s[i])+' nm\n'
        else:
            if i%2:#+str(round(rad[i],2)) str(fac*round(s[i],2))
                string = string + txt + ' r={:.2f}'.format(rad[i]) +' s={:.2f}'.format(fac*s[i])+'\n'
            else:
                string = string + txt + ' r={:.2f}'.format(rad[i]) +' s={:.2f}'.format(fac*s[i])+'   '
    #ax.text(xlim[1],ylim[0],  neighbour_string, fontsize=8) #verticalalignment='top'
    axs[1].set_ylim([0,1])
    if len(n) == 0:
        print('no residues below line')
        axs[1].text(0,0.5,  string, fontsize=fonts, verticalalignment='center')
    elif len(n)<11:
        axs[1].text(0,1,  string, fontsize=fonts, verticalalignment='top') #1.05*xlim[1],0.7*ylim[0]
    else:
        axs[1].text(0,1,  string, fontsize=fonts, verticalalignment='top') #1.05*xlim[1],1*ylim[0]
    axs[1].set_axis_off()

    #ax.legend(loc="best")
    axs[0].set_ylabel("Pore radius (nm)", fontsize=f_size)
    axs[0].set_xlabel(r"Hydrophobicity", fontsize=f_size) #($\AA$)
    axs[0].tick_params(axis='both', which='major', labelsize=f_size)
    fig.tight_layout()
    fig.savefig(path + title + ".png", bbox_inches='tight')
    plt.show()

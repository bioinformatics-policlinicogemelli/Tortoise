import os
PATH_ALL_GENES="./../../NEW_DATA_NETWORK/Output/LUNG/Lung_Gene_Particular/EGFR_all/Arricchimento_egfr_del_19/"
#PATH_GENES_UNIQUE= "../../../DATA/POLMONE_NO_VUS_2/OUTPUT/Arricchimento/"

def plot_term_go(file,term_go,size_1=30,size_2=15):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})

    data=pd.read_csv(file,header=0)
    data=data[data["Adjusted.P.value"].astype(float)<0.05]
    data=data.sort_values(by="Adjusted.P.value")
    data=data[0:20]
    data=data[::-1]

    list_adjusted_term=[]
    min_adjp=0.05
    max_adjp=0
    for _r in data.iterrows():
        term=str(_r[1][0])
        adjusted=float(_r[1][3])
        gene=str(_r[1][-1]).split(";")
        list_adjusted_term.append([term,adjusted,len(gene)])
        if adjusted < min_adjp:
            min_adjp=adjusted
        
        if adjusted>max_adjp:
            max_adjp=adjusted
    #print(list_adjusted_term)
    colors = plt.cm.cool(np.linspace(0, 1, len(list_adjusted_term)))
    # mappa il valore del p-value su una scala di colori
    color_values = [np.interp(p[1], [min_adjp, max_adjp], [0, len(colors)-1]) for p in list_adjusted_term]
    # seleziona il colore corrispondente per ogni barra
    bar_colors = [colors[int(c)] for c in color_values]
    fig, ax = plt.subplots(figsize=(size_1,size_2))
    fig.subplots_adjust(left=0.3, right=1)
    
 
    ax.barh([_t[0].split("(GO:")[0].replace(",","\n") for _t in list_adjusted_term],[_t[-1] for _t in list_adjusted_term],color=bar_colors)
   
    ax.set_xlabel("Gene Count")
    #ax.set_title(f"Cluster {cluster}", fontsize=18)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(vmin=min_adjp, vmax=max_adjp), cmap=plt.cm.cool), ax=ax)
    # Definisce i valori estremi della scala del gradiente di colore
    cbar.set_ticks([min_adjp,max_adjp])
    # Aggiunge una label alla barra dei colori
    cbar.ax.set_ylabel('Adjusted p value', rotation=-90, va='bottom',fontsize=30)
    # Aggiunge una legenda accanto al grafico
    #fig.subplots_adjust(right=0.8)
    #cax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #cax.imshow(np.array([[np.min(adjusted)], [np.max(adjusted)]]), cmap=plt.cm.cool)
    #cax.set_yticks([np.min(adjusted), np.max(adjusted)])
    #cax.set_yticklabels(['{:.2f}'.format(np.min(adjusted)), '{:.2f}'.format(np.max(adjusted))])
    #cax.yaxis.tick_right()
    plt.savefig(f"{PATH_ALL_GENES}Images/{term_go}.png")



def plot_term_kegg(file,term_go,size_1=30,size_2=15):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})

    data=pd.read_csv(file,header=0)
    print(data.columns)
    data=data[data["KEGG_2021_Human.Adjusted.P.value"].astype(float)<0.05]
    data=data.sort_values(by="KEGG_2021_Human.Adjusted.P.value")
    data=data[0:20]
    data=data[::-1]

    list_adjusted_term=[]
    min_adjp=0.05
    max_adjp=0
    for _r in data.iterrows():
        term=str(_r[1][0])
        adjusted=float(_r[1][3])
        gene=str(_r[1][-1]).split(";")
        list_adjusted_term.append([term,adjusted,len(gene)])
        if adjusted < min_adjp:
            min_adjp=adjusted
        
        if adjusted>max_adjp:
            max_adjp=adjusted
    #print(list_adjusted_term)
    colors = plt.cm.cool(np.linspace(0, 1, len(list_adjusted_term)))
    # mappa il valore del p-value su una scala di colori
    color_values = [np.interp(p[1], [min_adjp, max_adjp], [0, len(colors)-1]) for p in list_adjusted_term]
    # seleziona il colore corrispondente per ogni barra
    bar_colors = [colors[int(c)] for c in color_values]
    fig, ax = plt.subplots(figsize=(size_1,size_2))
    fig.subplots_adjust(left=0.3, right=1)
    
 
    ax.barh([_t[0].replace(",","\n") for _t in list_adjusted_term],[_t[-1] for _t in list_adjusted_term],color=bar_colors)
   
    ax.set_xlabel("Gene Count")
    #ax.set_title(f"Cluster {cluster}", fontsize=18)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(vmin=min_adjp, vmax=max_adjp), cmap=plt.cm.cool), ax=ax)
    # Definisce i valori estremi della scala del gradiente di colore
    cbar.set_ticks([min_adjp,max_adjp])
    # Aggiunge una label alla barra dei colori
    cbar.ax.set_ylabel('Adjusted p value', rotation=-90, va='bottom',fontsize=30)
    # Aggiunge una legenda accanto al grafico
    #fig.subplots_adjust(right=0.8)
    #cax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #cax.imshow(np.array([[np.min(adjusted)], [np.max(adjusted)]]), cmap=plt.cm.cool)
    #cax.set_yticks([np.min(adjusted), np.max(adjusted)])
    #cax.set_yticklabels(['{:.2f}'.format(np.min(adjusted)), '{:.2f}'.format(np.max(adjusted))])
    #cax.yaxis.tick_right()
    plt.savefig(f"{PATH_ALL_GENES}Images/{term_go}.png")


def plot_term_pheno(file,term_go,size_1=30,size_2=15):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})

    data=pd.read_csv(file,header=0)
    print(data.columns)
    data=data[data["PhenGenI_Association_2021.Adjusted.P.value"].astype(float)<0.05]
    data=data.sort_values(by="PhenGenI_Association_2021.Adjusted.P.value")
    data=data[0:20]
    data=data[::-1]

    list_adjusted_term=[]
    min_adjp=0.05
    max_adjp=0
    for _r in data.iterrows():
        term=str(_r[1][0])
        adjusted=float(_r[1][3])
        gene=str(_r[1][-1]).split(";")
        list_adjusted_term.append([term,adjusted,len(gene)])
        if adjusted < min_adjp:
            min_adjp=adjusted
        
        if adjusted>max_adjp:
            max_adjp=adjusted
    #print(list_adjusted_term)
    colors = plt.cm.cool(np.linspace(0, 1, len(list_adjusted_term)))
    # mappa il valore del p-value su una scala di colori
    color_values = [np.interp(p[1], [min_adjp, max_adjp], [0, len(colors)-1]) for p in list_adjusted_term]
    # seleziona il colore corrispondente per ogni barra
    bar_colors = [colors[int(c)] for c in color_values]
    fig, ax = plt.subplots(figsize=(size_1,size_2))
    fig.subplots_adjust(left=0.3, right=1)
    
 
    ax.barh([_t[0].replace(",","\n") for _t in list_adjusted_term],[_t[-1] for _t in list_adjusted_term],color=bar_colors)
   
    ax.set_xlabel("Gene Count")
    #ax.set_title(f"Cluster {cluster}", fontsize=18)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(vmin=min_adjp, vmax=max_adjp), cmap=plt.cm.cool), ax=ax)
    # Definisce i valori estremi della scala del gradiente di colore
    cbar.set_ticks([min_adjp,max_adjp])
    # Aggiunge una label alla barra dei colori
    cbar.ax.set_ylabel('Adjusted p value', rotation=-90, va='bottom',fontsize=30)
    # Aggiunge una legenda accanto al grafico
    #fig.subplots_adjust(right=0.8)
    #cax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #cax.imshow(np.array([[np.min(adjusted)], [np.max(adjusted)]]), cmap=plt.cm.cool)
    #cax.set_yticks([np.min(adjusted), np.max(adjusted)])
    #cax.set_yticklabels(['{:.2f}'.format(np.min(adjusted)), '{:.2f}'.format(np.max(adjusted))])
    #cax.yaxis.tick_right()
    plt.savefig(f"{PATH_ALL_GENES}Images/{term_go}.png")


def plot_term_wiki(file,term_go,size_1=30,size_2=15):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    matplotlib.rcParams.update({'font.size': 18})

    data=pd.read_csv(file,header=0)
    print(data.columns)
    data=data[data["WikiPathway_2023_Human.Adjusted.P.value"].astype(float)<0.05]
    data=data.sort_values(by="WikiPathway_2023_Human.Adjusted.P.value")
    data=data[0:20]
    data=data[::-1]

    list_adjusted_term=[]
    min_adjp=0.05
    max_adjp=0
    for _r in data.iterrows():
        term=str(_r[1][0])
        adjusted=float(_r[1][3])
        gene=str(_r[1][-1]).split(";")
        list_adjusted_term.append([term,adjusted,len(gene)])
        if adjusted < min_adjp:
            min_adjp=adjusted
        
        if adjusted>max_adjp:
            max_adjp=adjusted
    #print(list_adjusted_term)
    colors = plt.cm.cool(np.linspace(0, 1, len(list_adjusted_term)))
    # mappa il valore del p-value su una scala di colori
    color_values = [np.interp(p[1], [min_adjp, max_adjp], [0, len(colors)-1]) for p in list_adjusted_term]
    # seleziona il colore corrispondente per ogni barra
    bar_colors = [colors[int(c)] for c in color_values]
    fig, ax = plt.subplots(figsize=(size_1,size_2))
    fig.subplots_adjust(left=0.3, right=1)
    
 
    ax.barh([_t[0].replace(",","\n") for _t in list_adjusted_term],[_t[-1] for _t in list_adjusted_term],color=bar_colors)
   
    ax.set_xlabel("Gene Count")
    #ax.set_title(f"Cluster {cluster}", fontsize=18)
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(vmin=min_adjp, vmax=max_adjp), cmap=plt.cm.cool), ax=ax)
    # Definisce i valori estremi della scala del gradiente di colore
    cbar.set_ticks([min_adjp,max_adjp])
    # Aggiunge una label alla barra dei colori
    cbar.ax.set_ylabel('Adjusted p value', rotation=-90, va='bottom',fontsize=30)
    # Aggiunge una legenda accanto al grafico
    #fig.subplots_adjust(right=0.8)
    #cax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #cax.imshow(np.array([[np.min(adjusted)], [np.max(adjusted)]]), cmap=plt.cm.cool)
    #cax.set_yticks([np.min(adjusted), np.max(adjusted)])
    #cax.set_yticklabels(['{:.2f}'.format(np.min(adjusted)), '{:.2f}'.format(np.max(adjusted))])
    #cax.yaxis.tick_right()
    plt.savefig(f"{PATH_ALL_GENES}Images/{term_go}.png")






for f in os.listdir(PATH_ALL_GENES):
    if os.path.isfile(os.path.join(PATH_ALL_GENES,f)):
        type_go=f.split("_")[0]
        #=f.split("_")[1].split(".")[0]
        if f.startswith("kegg"):
            plot_term_kegg(f"{PATH_ALL_GENES}{f}",type_go)
        elif f.startswith("phen"):
            plot_term_pheno(f"{PATH_ALL_GENES}{f}",type_go)
        
        elif f.startswith("wiki"):
            plot_term_wiki(f"{PATH_ALL_GENES}{f}",type_go)
        
        else:
            plot_term_go(f"{PATH_ALL_GENES}{f}",type_go)


import pandas as pd
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import seaborn as sns
import random as random
import random
import datetime
import sys
from multiprocessing import Pool
from itertools import combinations
import matplotlib.pyplot as plt
from comut import comut
from comut import fileparsers
import os
import colorsys




# IMPORTIAMO I DATI
# - MUTAZIONALI
# - CLINICI DEL PAZIENTE
# - CLINICI DEL SAMPLE
def read_file(path_mutational, path_clinical_sample="", path_clinical_patient=""):
    data_mutational=pd.read_csv(f"{path_mutational}",sep="\t",low_memory=False)
    
    if path_clinical_sample!="":
        data_sample=pd.read_csv(f"{path_clinical_sample}",sep="\t")
    else:
        data_sample=[]

    if path_clinical_patient!="":
        data_patient=pd.read_csv(f"{path_clinical_patient}",sep="\t")
    else:
        data_patient=[]

    return data_mutational,data_sample,data_patient

#FIXME
# FUNZIONE PER PLOTTARE LE INFORMAZIONI CHE VOGLIAMO A PARTIRE DAI DATI CLINICI DEL PAZIENTE E/O SAMPLE
def summury_dataset(dataset,column_name,path_saved):
    import os

    path_folder=f"{path_saved}/Clinical_Figure/"
    if not os.path.exists(path_folder):
        os.makedirs(path_folder)
    if column_name not in list(dataset.columns):
        print(f"Colonna {column_name} non presente nel dataset")
    else:
        labels=list(dataset[column_name].unique())
        ax=plt.subplot()
        
        #ax.pie(list(dataset[column_name].value_counts(dropna=False)),labels=labels,autopct='%1.1f%%',pctdistance=1.1, labeldistance=.6,textprops={'size': 'smaller'})
        ax.pie(list(dataset[column_name].value_counts(dropna=False)),labels=labels,autopct='%1.1f%%')
        plt.title(f"Descriptive data for {column_name}")
        plt.savefig(f"{path_folder}/{column_name}_descriptive.png")
        plt.show()


#funzione per aggiungere una colonna che abbia nome del gene e sostituzione amminoacidica
def adding_category_mutation(data_mutational,gene_name,hgsvp_short,variant_classification,hgvsc):
    print(hgsvp_short,variant_classification,hgvsc)
    if hgsvp_short=="None" and hgvsc!="":
        nuovi_nomi={}
        data_mutational.fillna({f'{gene_name}': 'N/D', f'{variant_classification}': 'N/D',f'{hgvsc}': 'N/D'}, inplace=True)
        gruppi_mutazioni = data_mutational.groupby([f'{gene_name}',f'{variant_classification}',f'{hgvsc}'])
        for nome_gruppo, gruppo in gruppi_mutazioni:
            if nome_gruppo[2]=="N/D":
                nuovi_nomi[nome_gruppo]= f'Mut_{nome_gruppo[0]}_{nome_gruppo[1]}'
            else:
                nuovi_nomi[nome_gruppo]= f'Mut_{nome_gruppo[0]}_{nome_gruppo[2]}'
        data_mutational['nome_mutazione'] = data_mutational.apply(lambda row: nuovi_nomi[(row[f'{gene_name}'], row[f'{variant_classification}'],row[f'{hgvsc}'])], axis=1)
        return data_mutational

 
    elif hgvsc=="None" and hgsvp_short!= "":
        print("ciao")
        nuovi_nomi={}
        data_mutational.fillna({f'{gene_name}': 'N/D', f'{hgsvp_short}': 'N/D', f'{variant_classification}': 'N/D'}, inplace=True)
        gruppi_mutazioni = data_mutational.groupby([f'{gene_name}', f'{hgsvp_short}', f'{variant_classification}'])
        for nome_gruppo, gruppo in gruppi_mutazioni:
            if nome_gruppo[1]=="N/D":
                nuovi_nomi[nome_gruppo]= f'Mut_{nome_gruppo[0]}_{nome_gruppo[2]}'
            else:
                nuovi_nomi[nome_gruppo]= f'Mut_{nome_gruppo[0]}_{nome_gruppo[1]}'
        data_mutational['nome_mutazione'] = data_mutational.apply(lambda row: nuovi_nomi[(row[f'{gene_name}'], row[f'{hgsvp_short}'], row[f'{variant_classification}'])], axis=1)
        return data_mutational
    
    elif (hgsvp_short=="None") and (hgvsc=="None"):
        print("Hgvsp_Short e Hgvsc non present")
        

    #return data_mutational
    #data_mutational.fillna({'Hugo_Symbol': 'N/D', 'HGVSp_Short': 'N/D', 'Variant_Classification': 'N/D','HGVSc': 'N/D'}, inplace=True)
    #data_mutational.fillna({f'{gene_name}': 'N/D', f'{hgsvp_short}': 'N/D', f'{variant_classification}': 'N/D',f'{hgvsc}': 'N/D'}, inplace=True)
    #gruppi_mutazioni = data_mutational.groupby([f'{gene_name}', f'{hgsvp_short}', f'{variant_classification}',f'{hgvsc}'])
    #print(gruppi_mutazioni)
    #nuovi_nomi={}
    
    #for nome_gruppo, gruppo in gruppi_mutazioni:
       # try:
            #nuovi_nomi[nome_gruppo]=f'Mut_{nome_gruppo[0]}'

          #  if nome_gruppo[1]=="N/D" and nome_gruppo[2]=="N/D":
              # nuovi_nomi[nome_gruppo]= f'Mut_{nome_gruppo[0]}_{nome_gruppo[3]}'

          #  if nome_gruppo[1]=="N/D" and nome_gruppo[2]!="N/D":
             # nuovi_nomi[nome_gruppo]= f'Mut_{nome_gruppo[0]}_{nome_gruppo[2]}'
#
           # else:
                #print(nome_gruppo)
              #  nuovi_nomi[nome_gruppo]= f'Mut_{nome_gruppo[0]}_{nome_gruppo[1]}'
       # except:
          #  print(f"Not present {hgsvp_short} or {hgvsc} in your dataset")

  
    #data_mutational['nome_mutazione'] = data_mutational.apply(lambda row: nuovi_nomi[(row[f'{gene_name}'], row[f'{hgsvp_short}'], row[f'{variant_classification}'],row[f'{hgvsc}'])], axis=1)

    #return data_mutational

#funzione per cacolare, se assente, la colonna della VAF
def calculated_vaf(riga):
    return (riga['t_alt_count'])/(riga['t_alt_count'] + riga["t_ref_count"]) 



#FIXME
#funzione per contare le mutazioni del gene target di interesse e scrivere tali info in un file csv
def check_target_mutation_all_tissue(data_mutational, gene_target,sample_id=[],primary=False):
    gene_mutation={}
    list_gene=[]
    for row in data_mutational.iterrows():
        sample=row[1]["Tumor_Sample_Barcode"]
        if primary:
            if str(row[1]["Hugo_Symbol"])==gene_target and sample in sample_id:
                gene=str(row[1]["Hugo_Symbol"])+"##"+str(row[1]["Chromosome"])+"##"+str(row[1]["Start_Position"])+"##"+str(row[1]["End_Position"])+"##"+str(row[1]["HGVSp_Short"])+"##"+str(row[1]["Exon_Number"])
                list_gene.append(gene)
        else:
            if str(row[1]["Hugo_Symbol"])==gene_target:
                gene=str(row[1]["Hugo_Symbol"])+"##"+str(row[1]["Chromosome"])+"##"+str(row[1]["Start_Position"])+"##"+str(row[1]["End_Position"])+"##"+str(row[1]["HGVSp_Short"])+"##"+str(row[1]["Exon_Number"])
                list_gene.append(gene)


    for element in list_gene:
        numerosity=list_gene.count(element)
        if element not in gene_mutation.keys():
            gene_mutation[element]=0
        gene_mutation[element]=numerosity

    gene_mutation=dict(sorted(gene_mutation.items(), key=lambda x: x[1], reverse=True))
    
    name_file=""
    if primary:
        name_file=f"./Target_Gene_Infos/mutation_{gene_target}_primary.csv"
    else:
        name_file=f"./Target_Gene_Infos/mutation_{gene_target}.csv"

    with open(name_file,"w") as f:
        f.write("Gene\tChromosome\tStart\tEnd\tSost_amm\tExon\tCount\tPercentage\n")
        for gene,count in gene_mutation.items():
            total=len(list_gene)
            percent=round(((count*100)/total),4)
            gene=gene.split("##")
            for g in gene:
                f.write(g+"\t")
            f.write(str(count)+"\t"+str(percent)+"\n")

#FIXME
#plot count geni target
def plot_target_mutation(filename):
    data_egfr_primary = pd.read_csv(filename,sep="\t")
    gene=filename.split("_")[1]
    #plottiamo le prime 13 mutazioni la cui comparsa va da 205 a 10
    #data_egfr_sub=data_egfr.head(13)
    plt.figure(figsize=(15, 10)) 
    sns.barplot(data_egfr_primary,x="Count",y="Sost_amm",hue="Sost_amm")
    plt.title(f"{gene} Mutation Primary Tissue")
    plt.yticks(fontsize=9)
    plt.xticks(fontsize=9)
    plt.savefig(f"./Target_Gene_Infos/plot_{gene}_mutation")


#funzione per creare mappa paziente e mutazione con categorizzazione
def create_maps(data_mutational,data_config,column_mutation,gene_interest=""):
    map_variants={}
    map_patients={}
    map_variant_count={}
    map_tp53={}
    map_consequence={}
    for _row in data_mutational.iterrows():
    #GET INFOS FROM DATA
        _paz=str(_row[1]["Tumor_Sample_Barcode"])
        _gene=str(_row[1][data_config["Mutation"]["column_gene_name"]])
        _chrom = str(_row[1]["Chromosome"])
        _chrom_start = str(_row[1]["Start_Position"])
        _chrom_stop = str(_row[1]["End_Position"])
        _cons=str(_row[1]["Variant_Classification"])
        _variant_type=str(_row[1]["Variant_Type"])
        
        #se la colonna nome mutazione non è presente nel dataset, allora verrà applicata la categorizzazione di default
        if column_mutation=="":
            #data_mutational=adding_category_mutation(data_mutational)
            _category=str(_row[1]["nome_mutazione"])
        #altrimenti verrà preso come nome mutazione, quello presente nella colonna indicata dall'utente
        else:
            _category=str(_row[1][column_mutation])
        
        #se questa colonna è vuota tornerà una stringa vuota
        if "HGVSp_Short" in data_mutational.columns:
            _sost_amm=str(_row[1]["HGVSp_Short"])
        else:
            _sost_amm=""
        #se questa colonna è vuota tornerà una stringa vuota (dire all'utente che la colonna corrispondente alla VAF dovrà essere chiamata t_AF)
        if "t_AF" in data_mutational.columns:
             _vaf=float(_row[1]["t_AF"])
        else:
            _vaf=""
        
        #assegnazione del nome della variante
        _var=_category
       
        #creazione dizionario delle varianti
        if _var not in map_variants.keys():
            map_variants[_var]={}
            map_variants[_var]["patients"]=set()
            map_variants[_var]["sost_amm"]=_sost_amm
            map_variants[_var]["gene"]=_gene
            map_variants[_var]["cons"]=_cons
            map_variants[_var]["vaf"]=_vaf
            
        map_variants[_var]["patients"].add(_paz)

        #creazione dizionario dei pazienti
        if _paz not in map_patients.keys():
            map_patients[_paz]={}
            map_patients[_paz]["variants"]=set()
        map_patients[_paz]["variants"].add(_var)

        #creazione dizionario delle consequence
        if _gene not in map_consequence.keys():
            map_consequence[_gene]=[]
        map_consequence[_gene].append(_cons)

    return map_patients,map_variants, map_consequence

#plot delle mutazioni dei geni di interesse

def plot_mutation_gene(map_consequence,list_gene_interest,path_save):
    gene_mutation_list = [(gene, mutation) for gene, mutations in map_consequence.items() for mutation in mutations if gene in list_gene_interest]
    df_mut = pd.DataFrame(gene_mutation_list, columns=['gene', 'mutazione'])
    pivot_df = df_mut.pivot_table(index='gene', columns='mutazione', aggfunc='size', fill_value=0)
    pivot_df.plot(kind='bar', stacked=True)
    plt.xlabel('Gene')
    plt.ylabel('Count of Mutations')
    plt.title('Stacked Barplot of Mutation Types per Gene')
    plt.legend(title='Mutation Types', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f"{path_save}/class_mutation_gene_interest.png")

#Creazione del Grafo
def graph_creation(map_patients,map_variants):
    edges=[]
    for _k_variant, _v_variant in map_variants.items():
        for _k_patient, _v_patient in map_patients.items():
            if _k_patient in [_x for _x in _v_variant["patients"]]:
                edges.append((_k_variant, _k_patient))

    graph=ig.Graph()
    graph.add_vertices(list(map_variants.keys()),attributes={"vertex_type":"VARIANT","color":"blue","gene":[variant.split("_")[1]for variant  in map_variants.keys()],"sost_amm":[f'{value["gene"]}_{value["sost_amm"]}' for key,value  in map_variants.items()]})
    graph.add_vertices(list(map_patients.keys()),attributes={"vertex_type":"PATIENT","color":"red"})

    graph.add_edges(edges)
    return graph




#count dei geni presenti in relazione alle singole mutazioni
def count_gene(graph):
    gene_total_count={}
    for vertex in graph.vs:
        if vertex["vertex_type"]=="VARIANT":
            if vertex["gene"] not in gene_total_count.keys():
                gene_total_count[vertex["gene"]]=0
            gene_total_count[vertex["gene"]]+=1
    gene_total_count=dict(sorted(gene_total_count.items(), key=lambda kv: kv[1],reverse=True))
    return gene_total_count



# SELEZIONE DEL SEED CHE Dà VALORE DI MODULARITà PIù ALTA A SEGUITO DEL LEIDEN ALGORITHM
def selected_seed(GRAPH):
    best_seed=0
    if False :#sys.platform.startswith('win') or sys.platform.startswith("linux"):
        def process_data(_seed):
            random.seed(_seed)
            _dendro_2=GRAPH.community_leiden(objective_function="modularity")
            modularity=_dendro_2.modularity
            return modularity

        data = list(range(10_000))
        #print(__name__)
        with Pool() as p:
            mod_results = p.map(process_data, data)

        best_seed = mod_results.index(max(mod_results))

    else:
        def process_data(_seed):
            random.seed(_seed)
            _dendro_2=GRAPH.community_leiden(objective_function="modularity")
            modularity=_dendro_2.modularity
            return modularity

        data = list(range(10_000))
        mod_results = []
        for s in data:
            mod_results.append(process_data(s))
        best_seed = mod_results.index(max(mod_results))


    best_seed = mod_results.index(max(mod_results))
    return best_seed


# SELEZIONE DEL SEED E LANCIO DELL'ALGORITMO DI CLUSTERIZZAZIONE
def leiden_clustering(graph, best_seed):
    random.seed(best_seed)
    dendro=graph.community_leiden(objective_function="modularity")
    print("numero di clusters:", len(list(dendro)),"Modularità:",dendro.modularity)
    return dendro


def adding_graph_color(graph,dendro):
    num_clusters = len(dendro)
    # Genera colori in base al numero di cluster
    colors = []
    for i in range(num_clusters):
        hue = i / num_clusters  # variare la tonalità in base al numero di cluster
        rgb = colorsys.hsv_to_rgb(hue, 1, 1)  # converti da spazio colore HSV a RGB
        colors.append('#%02x%02x%02x' % tuple(int(c * 255) for c in rgb))  # formato esadecimale RGB
    # Assegna i colori alle singole comunità nel grafo
    graph.vs["color"] = [colors[cluster] for cluster in dendro.membership]

    return graph

def plot_graph(graph,path_save,gene):
    ig.plot(graph,f"{path_save}/plot_{gene}.pdf",
        **{
        "vertex_color": graph.vs["color"],
        "bbox": (2000,1000),
        "edge_curved": 0.1,
        "vertex_label":  graph.vs["name"],
        "vertex_color": graph.vs["color"],
        "vertex_label_size": 1.5,
        "vertex_size":20,
        "edge_color":"grey",
        "vertex_shape":["triangle" if v["vertex_type"]=="PATIENT" else "circle" for v in graph.vs ],
        #"edge_color": [get_color(x, _max_conn, _min_conn) for x in _new_g.es["weight"]],
        "edge_widht":0.3,
        "vertex_frame_width":0.05,
        "layout":"fr"
        })

def write_graph_to_cytoscape(graph,path_save):
    graph.write_graphml(f"{path_save}/grafo_cytoscape.graphml")



# CREAZIONE DI UN DIZIONARIO MAP_CLUSTER 
def map_cluster_creation(graph,dendro):
    map_cluster={}
    for cluster_index in range(len(dendro)):
        map_cluster[cluster_index]=[]
        for  element in dendro[cluster_index]:
            _vertex = graph.vs()[element]
            map_cluster[cluster_index].append((_vertex["name"],_vertex["vertex_type"]))
            
    return map_cluster

#funzione per aggiungere il cluster alla mappa dei pazienti e delle varianti
def adding_cluster_to_map(map_cluster,map_patients,map_variants):
    for cluster, infos in map_cluster.items():
        for info in infos:
            if info[1]== "PATIENT":
                map_patients[info[0]]["cluster"]=cluster
            elif info[1]=="VARIANT":
                map_variants[info[0]]["cluster"]=cluster
                
    return map_patients,map_variants

#funzione per scrivere i centroidi di ciascun cluster in un file
def centroids_cluster(dendro,path_save):
    with open (f"{path_save}/Centroidi_Mutazioni.csv","w") as f:
        for _i in range(len(dendro)):
            sub_graph = dendro.subgraph(_i)
            max_value = 0
            _list_centroids = []
            for _v in sub_graph.vs():
                print(_v["name"])
                if _v["vertex_type"] != "VARIANT":
                    continue
                temp_val = sub_graph.neighborhood_size(_v)
                if temp_val > max_value:
                    max_value = temp_val
                    _list_centroids = [(_v["name"],max_value-1)]
                elif temp_val == max_value:
                    _list_centroids.append((_v["name"],max_value-1))
            #print(path_save)
            f.write(f"CLUSTER {_i} CENTROIDS: {_list_centroids}\n")

#funzione per scrivere e trovare le coppie centroide-elemento per ciascun cluster
def couple_centroid_element(dendro,map_cluster,path_save):
    all_pairs=[]
    for _i in range(len(dendro)):
        sub_graph = dendro.subgraph(_i)
        max_value = 0
        _list_centroids = []
        #individuazione centroide per ogni cluster
        for _v in sub_graph.vs():
            temp_val = sub_graph.neighborhood_size(_v)
            if temp_val > max_value:
                max_value = temp_val
                _list_centroids = [(_v["name"],max_value-1)]
            elif temp_val == max_value:
                _list_centroids.append((_v["name"],max_value-1))
        
        _list_pairs=[]
        #creazione coppie centroide-elemento per ogni cluster
        for element in map_cluster[_i]:
            if element[0] not in [ x[0] for x in _list_centroids]:
                for c in _list_centroids:
                    _list_pairs.append((c[0],element[0]))
        all_pairs+=_list_pairs
    with open(f"{path_save}/coppie_centroide.csv","w") as f:
        f.write("Centroide\tElemento\n")
        for elements in all_pairs:
            f.write(elements[0]+"\t"+elements[1]+"\n")

    return all_pairs



# AGGIUNTA DELLE INFORMAZIONI CLINICHE DI INTERESSE ALLA MAPPA DEI PAZIENTI + AGGIUNTA DEL CLUSTER DI APPARTENENZA ALLA MAPPA DEI PAZIENTI E DELLE MITAZIONI
def enriched_sample_data(data_sample,MAP_CLUSTER,MAP_PATIENTS,sample_name,patient_name):
    #generalizzazione --> CLINICAL SAMPLE
    list_clinical_parameters=[col for col in data_sample.columns]
    for row in data_sample.iterrows():
        _sample=row[1][sample_name]
        _paz=row[1][patient_name]
        for parameter in list_clinical_parameters:
            variable=row[1][parameter]
            if _sample in MAP_PATIENTS.keys():
                MAP_PATIENTS[_sample][parameter]=variable
    return MAP_PATIENTS

def enriched_patient_data(data_patient,MAP_PATIENTS,patient_name):
    #generalizzazione --> CLINICAL PATIENT
    list_clinical_patient = [col for col in data_patient.columns if col != patient_name]

    for _row in data_patient.iterrows():
            _paz=_row[1][patient_name]
            for parameters in list_clinical_patient:

                variable=_row[1][parameters]
                for _sample,value in MAP_PATIENTS.items():
                    if value[patient_name]==_paz:
                        MAP_PATIENTS[_sample][parameters]=variable
    return MAP_PATIENTS


#funzione per assegnare alle mutazioni dei geni target, il cluster di appartenenza
def cluster_target_file(data_target,file_name,MAP_VARIANTS):
    temp=[]
    for row in data_target.iterrows():
        #mutation=row[1]["Gene"]+"_"+str(row[1]["Chromosome"])+"_"+str(row[1]["Start"])+"_"+str(row[1]["End"])
        mutation=row[1]["Gene"]+"_"+row[1]["Sost_amm"]
        cluster=MAP_VARIANTS[mutation]["cluster"]
        temp.append(cluster)

    data_target["Cluster"]=temp
    data_target.to_csv(file_name,index=False,sep="\t")
    return data_target

#funzione per assegnatre l'attributo oncotreecode al grafo
def attributed_oncotreecode_graph(GRAPH,MAP_PATIENTS):
    for v in GRAPH.vs():
        if v["vertex_type"]=="PATIENT":
            v["tumor"]=MAP_PATIENTS[v["name"]]["ONCOTREE_CODE"]


#FIXME da migliorare la parte dell'algoritmo
def percentage_tumor_cluster(data_sample,MAP_PATIENTS):
    from columnar import columnar
    TUMORS=list(data_sample["ONCOTREE_CODE"].unique())
    if len(TUMORS)>1:
        _map_cluster_tumor={}
        for _k, _v in MAP_PATIENTS.items():
            _cluster_id=_v["cluster"]
            _tumor=_v["ONCOTREE_CODE"]
            if _cluster_id not in _map_cluster_tumor.keys():
                _map_cluster_tumor[_cluster_id]={}
            if _tumor not in _map_cluster_tumor[_cluster_id].keys():
                _map_cluster_tumor[_cluster_id][_tumor]=0
            _map_cluster_tumor[_cluster_id][_tumor]+=1
        _map_cluster_tumor=dict(sorted(_map_cluster_tumor.items(),key=lambda x:x[0]))
        ## PRINT TABLE
        _data_table=[]
        for _k,_v in _map_cluster_tumor.items():
            _temp=[]
            _temp.append(_k)
            for _t in TUMORS:
                try:
                    _temp.append(_map_cluster_tumor[_k][_t])
                except:
                     _temp.append(0)
            _data_table.append(_temp)
        #table = columnar(_data_table, headers=['CLUSTER']+TUMORS,min_column_width=5,row_sep="_")
        #print(table)
        
        _data_table_percent=[]
        for _t in _data_table:
            _temp=[]
            _temp.append(_t[0])
            for _i in range(1,len(_t)):
                _percent=round((_t[_i]/sum(_t[1:]))*100,2)
                _temp.append(_percent)
            #percent_patient.append(_temp)
            _data_table_percent.append(_temp)
                
        table_2= columnar(_data_table_percent,headers=['CLUSTER']+TUMORS,row_sep="_")
        #print(table_2)
        
    else:
        print("Single Tumor")
    
    return table_2




# CREAZIONE DEL FILE "CLUSTER_CLINICAL_DATA" IN CUI INSERIAMO IL CLUSTER DI APPARTENENZA DI OGI SAMPLES + INFORMAZIONI CLINICHE
def creation_cluster_clinical_data(MAP_PATIENTS,path_saved):
    header_clinical_data=sorted(list(set([k2 for k,v in MAP_PATIENTS.items() for k2 in v.keys() if k2!="variants"])))
    with open (f"{path_saved}/cluster_clinical_data.csv","w") as f:
        f.write('\t'.join(header_clinical_data)+"\n")
        for k,v in MAP_PATIENTS.items():
            temp=[]
            for data in header_clinical_data:
                temp.append(str(v[data]))
            f.write("\t".join(temp)+"\n")
    df = pd.read_csv(f"{path_saved}/cluster_clinical_data.csv",sep="\t")
    return df


# creazione di un file con il numero di pazienti e varianti per ogni cluster  + CREAZIONE DI UN ARRAY CON I CLUSTER CON > 1 PAZIENTE + UN ARRAY CON I CLUSTER CON ==1 PAZIENTE
def cluster_division(MAP_CLUSTER):
    cluster_one_patient=[]
    cluster_more_patient=[]
    with open ("numerosity_cluster.csv","w") as f:
        f.write("Cluster\tPatient\tVariant\n")
        for cluster,value in MAP_CLUSTER.items():
            count_variant=0
            count_patient=0
            for v in value:
                if v[1]=="VARIANT":
                    count_variant+=1
                else:
                    count_patient+=1
            f.write(str(cluster)+"\t"+str(count_patient)+"\t"+str(count_variant)+"\n")
            if count_patient==1:
                cluster_one_patient.append(cluster)
            else:
                cluster_more_patient.append(cluster)

    return cluster_more_patient,cluster_one_patient






#defizione del numero di connessioni per ogni variante all'interno del cluster     
def variant_conncection_patient(_dendro_2):
    variant_patient_connection_count={}

    for _i in range(len(_dendro_2)):
        variant_patient_connection_count[_i]=[]
        sub_graph = _dendro_2.subgraph(_i)
        for _v in sub_graph.vs():
            if _v["vertex_type"] != "VARIANT":
                continue
            temp_val = sub_graph.neighborhood_size(_v)
            variant_patient_connection_count[_i].append((_v["name"],(temp_val-1)))
            
    variant_patient_connection_count={key:sorted(value, key=lambda x:x[1], reverse=True) for key, value in variant_patient_connection_count.items()}

    return variant_patient_connection_count


#defizione del numero di connessioni per ogni paziente all'interno del cluster
def patient_connection_variant(_dendro_2):
    patient_variant_connection_count={}
    for _i in range(len(_dendro_2)):
        patient_variant_connection_count[_i]=[]
        sub_graph = _dendro_2.subgraph(_i)
        for _v in sub_graph.vs():
            if _v["vertex_type"] != "PATIENT":
                continue
            temp_val = sub_graph.neighborhood_size(_v)
            patient_variant_connection_count[_i].append((_v["name"],(temp_val-1)))

    patient_variant_connection_count={key:sorted(value, key=lambda x:x[1], reverse=True) for key, value in patient_variant_connection_count.items()}     

    return patient_variant_connection_count



# CREAZIONE DI UN FILE "CONNECTION_VARIANT" IN CUI SONO INDICATE IL NUMERO DI VARIANTI COMUNI TRA I VARI PAZIENTI DI UN CLUSTER
def file_connection_variant(MAP_CLUSTER,MAP_PATIENTS,path_saved):
    with open(f"{path_saved}/connection_variant.csv","w") as f:
        f.write("Cluster\tPatient_1\tPatient_2\tVariant\tNumber_variant\n")
        for cluster,values in MAP_CLUSTER.items():
            #f.write(str(cluster)+"\t")
            patients=[]
            for value in values:
                if value[1]=="PATIENT":
                    patients.append(value[0])
            for patient in list(combinations(patients,2)):
                paz_1=patient[0]
                paz_2=patient[1]
                variant_p1=MAP_PATIENTS[paz_1]["variants"]
                variant_p2=MAP_PATIENTS[paz_2]["variants"]
                variant_common=list(variant_p1.intersection(variant_p2))
                f.write(str(cluster)+"\t"+str(paz_1)+"\t"+str(paz_2)+"\t"+",".join(variant_common)+"\t"+str(len(variant_common))+"\n")
                




# CREAZIONE DI UN FILE "CONCCECTION_PATIENT" IN CUI SONO INDICATI IL NUMERO DI PAZIENTI COMUNI TRA LE VARIE VARIANTI DI UN CLUSTER
def file_connection_patient(MAP_CLUSTER,MAP_VARIANTS,path_saved):
    with open(f"{path_saved}/connection_patient.csv","w") as f:
        f.write("Cluster\tVariant_1\tVariant_2\tPatient\tNumber_patient\n")
        for cluster,values in MAP_CLUSTER.items():
            #f.write(str(cluster)+"\t")
            variants=[]
            for value in values:
                if value[1]=="VARIANT":
                    variants.append(value[0])
            for variant in list(combinations(variants,2)):
                var_1=variant[0]
                var_2=variant[1]
                patient_v1=MAP_VARIANTS[var_1]["patients"]
                patient_v2=MAP_VARIANTS[var_2]["patients"]
                patient_common=list(patient_v1.intersection(patient_v2))
                f.write(str(cluster)+"\t"+str(var_1)+"\t"+str(var_2)+"\t"+",".join(patient_common)+"\t"+str(len(patient_common))+"\n")



#funzione che scrive un file che riassume le info presenti in ogni cluster      
def numerosity_info(path_save,map_cluster):
    with open (f"{path_save}/numerosity_cluster.csv","w") as f:
        f.write("Cluster\tPatient\tVariant\tGene\n")
        for cluster,value in map_cluster.items():
            count_variant=0
            count_patient=0
            gene_count=set()
            for v in value:
                if v[1]=="VARIANT":
                    count_variant+=1
                    gene=v[0].split("_")[1]
                    gene_count.add(gene)
                else:
                    count_patient+=1
            f.write(str(cluster)+"\t"+str(count_patient)+"\t"+str(count_variant)+"\t"+str(len(gene_count))+"\n")

 #riassunto delle informazioni (mutazioni e pazienti per ciascun cluster)
def summary_info(path_save,map_cluster,map_patients,patient_name):
    with open(f"{path_save}/summury_file.csv","w") as f:
        f.write("Cluster\tNode_Name\tType\tGene\tPatient_id\n")
        for cluster,infos in map_cluster.items():
            for info in infos:
                if info[1]=="VARIANT":
                    f.write(str(cluster)+"\t"+str(info[0])+"\t"+str(info[1])+"\t"+str(info[0].split("_")[1])+"\t"+"None"+"\n")
                if info[1]=="PATIENT":
                    try:
                        f.write(str(cluster)+"\t"+str(info[0])+"\t"+str(info[1])+"\t"+"None"+"\t"+str(map_patients[info[0]][patient_name])+"\n")  
                    except:
                        f.write(str(cluster)+"\t"+str(info[0])+"\t"+str(info[1])+"\t"+"None"+"\t"+str("nan")+"\n")  


#aggiunta size ai vertici mutazione 
def add_size_node(GRAPH,variant_patient_connection_count):        
    for cluster,elements in variant_patient_connection_count.items():
        for e in elements:
            variant=e[0]
            presence=e[1]
            _vx=GRAPH.vs.find(name=variant)
            _vx["size"]=presence
    return GRAPH


# COUNT DEI GENI PRESENTI IN OGNI CLUSTER + VALORE PERCENTUALE DI APPARTENENZA DI OGNI GENE A OGNI CLUSTER, CON SALVATAGGIO DELLE INFO IN UN FILE "DISTRIBUTION_GENE_CLUSTER"
def count_gene_abs_percent(map_cluster,gene_total_count,path_save):
    MAP_CLUSTER_GENE={}
    for cluster,values in map_cluster.items():
        _temp=[]
        if cluster not in MAP_CLUSTER_GENE.keys():
            MAP_CLUSTER_GENE[cluster]={}
        for value in values:
            if value[1]=="VARIANT":
                _gene=value[0].split("_")[1]
                _temp.append(_gene)
        for _g in _temp:
            count=_temp.count(_g)
            MAP_CLUSTER_GENE[cluster][_g]=count


    MAP_CLUSTER_GENE_PERCEN={}
    for cluster, genes in MAP_CLUSTER_GENE.items():
        if cluster not in MAP_CLUSTER_GENE_PERCEN.keys():
            MAP_CLUSTER_GENE_PERCEN[cluster]={}
        for gene, count in genes.items():
            if gene not in MAP_CLUSTER_GENE_PERCEN[cluster].keys():
                #print(gene)
                #print(gene_total_count[gene])
                MAP_CLUSTER_GENE_PERCEN[cluster][gene]=(MAP_CLUSTER_GENE[cluster][gene]*100)/gene_total_count[gene]

    MAP_CLUSTER_GENE_PERCEN={k: dict(sorted(v.items(), key=lambda x: x[1],reverse=True)) for k, v in MAP_CLUSTER_GENE_PERCEN.items()}

    with open(f"{path_save}/distribution_gene_cluster.csv","w") as f:
        f.write("Cluster\tGene\tPercentage\n")
        for cluster,genes in MAP_CLUSTER_GENE_PERCEN.items():
            for gene,count in genes.items():
                f.write(str(cluster)+"\t"+str(gene)+"\t"+str(count)+"\n")
    
    return MAP_CLUSTER_GENE,MAP_CLUSTER_GENE_PERCEN

# creazione di un file per ogni cluster, contenente i geni presenti
def genes_single_cluster(map_cluster,path_save):
    if not os.path.exists(f"{path_save}/Gene"):
        os.makedirs(f"{path_save}/Gene")

    for cluster, infos in map_cluster.items():
        with open(f"{path_save}/Gene/genes_cluster_{cluster}.csv","w") as f:
            set_gene=set()
            for info in infos:
                
                if info[1]=="VARIANT":
                    gene=info[0].split("_")[1]
                    set_gene.add(gene)
                #print(set_gene)
            for genes in set_gene:
                f.write(genes+"\n")


# AGGIUNTA DELL'ATTRIBUTO CLUSTER AI NODI DEL GRAFO
def cluster_noded_attributes(GRAPH,MAP_PATIENTS,MAP_VARIANTS):
    for v in GRAPH.vs():
        if v["vertex_type"]=="PATIENT":
            v["cluster"]=MAP_PATIENTS[v["name"]]["cluster"]
        else:
            v["cluster"]=MAP_VARIANTS[v["name"]]["cluster"]
    return GRAPH


# aggiunta dell'attributo "sost_amm" ai nodi del grafo
##TODO

def sost_amm_attributes(MAP_VARIANTS,GRAPH):
    for mutation,attributes in MAP_VARIANTS.items():
        gene=mutation.split("_")[0]
        for v in GRAPH.vs:
            if v["name"]==mutation:
                v["sost_amm"]=str(gene)+"_"+str(attributes["sost_amm"])

    


# CREAZIONE, PER OGNI CLUSTER, DI UN GRAFO I CUI NODI SONO LE VARIANTI E IN CUI:
# - LA GRANDEZZA DEI NODI è = AL NUMERO DI PAZIENTI DEL CLUSTER CHE HANNO QUELLA MUTAZIONE
# - 2 MUTAZIONI SONO CONNESSE SE CO-MUTATE (QUINDI PRESENTI IN ALMENO 2 PAZIENTI DEL CLUSTER)
# - GRANDEZZA DEI NODI è PARI AL NUMERO DI PAZIENTI CHE HANNO QUELLE DUE MUTAZIONI NEL CLUSTER

#connection_df = pd.read_csv("../Data/FILTRED/connection_patient.csv", sep="\t")
def get_color_comutated_cluster(value, _max, _min, _scaling=1):
    _range=(_max + 1 - _min) * _scaling
    _range_sector = _range / 4

    _index_color = (value+1)// _range_sector

    if _index_color==0:
        #RED
        return (1.0, 0.0, 0.0, 0.5)
    elif _index_color == 1:
        #ORANGE
        return (1.0, 0.5, 0.0, 0.5)
    elif _index_color == 2:
        #YELLOW
        return (1.0, 1.0, 0.0, 0.5)
    else:
        #GREEN
        return (0.0, 1.0, 0.0, 0.5)



def plot_distance_comutated_cluster_variants(_dendro_2,cluster_index,connection_df,path_saved):
    import os
    if not os.path.exists(f"{path_saved}/plot_mutation_connection"):
        os.makedirs(f"{path_saved}/plot_mutation_connection")

    sub_graph = _dendro_2.subgraph(cluster_index)
    list_variants = list(sub_graph.vs.select(vertex_type="VARIANT"))
    _new_g = ig.Graph()
    for var in list_variants:
        _new_g.add_vertex(name=var["name"], **{"size": min(100, var["size"]*5),"sost_amm":var["sost_amm"]})
    _temp_edge=set()

    _df_cluster = connection_df[connection_df["Cluster"]==cluster_index]
    _max_conn = max(_df_cluster["Number_patient"].values)
    _min_conn = min(_df_cluster["Number_patient"].values)

    for var in list_variants:
        _var_name = var["name"]
        #_var_name = var["sost_amm"]
        _filtered_edge = _df_cluster[(_df_cluster["Variant_1"]==_var_name) | (_df_cluster["Variant_2"]==_var_name)]

        for index, row in _filtered_edge.iterrows():
            _var1 = row["Variant_1"]
            _var2 = row["Variant_2"]
            _w = row["Number_patient"]
            if _w > 0:
                if _var1 > _var2:
                    _temp_edge.add((_var1, _var2, _w))
                else:
                    _temp_edge.add((_var2, _var1, _w))

    for e in _temp_edge:
        _new_g.add_edge(e[0], e[1], weight = e[2])
    
    weight_arch=[]
    
    if  len(_new_g.es)==0:
        weight_arch=[]
    else:
        weight_arch=_new_g.es["weight"]

    #fig, ax = plt.subplots(dpi=500)
    ig.plot(_new_g,
        #target=ax,
        f"{path_saved}/plot_mutation_connection/cluster_{cluster_index}.pdf",
        **{
        #"edge_width": _new_g.es["weight"],
        "edge_width":weight_arch,
        "vertex_color": "cyan",
        "bbox": (1290,820),
        "edge_curved": 0.1,
        #"vertex_label": [x.replace("_", "\n") for x in _new_g.vs["name"]],
        "vertex_label": [x.replace("_","\n") for x in _new_g.vs["sost_amm"]],
        "vertex_label_size": 2.5,
        #"edge_color": [get_color(x, _max_conn, _min_conn) for x in _new_g.es["weight"]],
        "edge_color": [get_color_comutated_cluster(x, _max_conn, _min_conn) for x in weight_arch],
        "edge_label": weight_arch,
        "edge_label_size": [1.5 + (n * 2) for n in weight_arch],
        "vertex_frame_width":0.05,
        "background":(0.3, 0.4, 0.5, 1)
        }
    )

'''
if sys.platform.startswith("win") or sys.platform.startswith("linux"):
    with Pool() as p:
        p.map(plot_distance, cluster_more_patient)
else:
    for i in cluster_more_patient:
        plot_distance(i)'''


## patient graph
#connection_df = pd.read_csv("connection_variant.csv", sep="\t")

def plot_distance_comutated_cluster_patients(_dendro_2,cluster_index,connection_df,path_saved):
    import os
    if not os.path.exists(f"{path_saved}/plot_patient_connection"):
        os.makedirs(f"{path_saved}/plot_patient_connection")

    sub_graph = _dendro_2.subgraph(cluster_index)
    list_patients = list(sub_graph.vs.select(vertex_type="PATIENT"))
    _new_g = ig.Graph()
    for paz in list_patients:
        _new_g.add_vertex(name=paz["name"])
    _temp_edge=set()

    _df_cluster = connection_df[connection_df["Cluster"]==cluster_index]
    _max_conn = max(_df_cluster["Number_variant"].values)
    _min_conn = min(_df_cluster["Number_variant"].values)

    for paz in list_patients:
        _paz_name = paz["name"]
        #_var_name = var["sost_amm"]
        _filtered_edge = _df_cluster[(_df_cluster["Patient_1"]==_paz_name) | (_df_cluster["Patient_2"]==_paz_name)]

        for index, row in _filtered_edge.iterrows():
            _var1 = row["Patient_1"]
            _var2 = row["Patient_2"]
            _w = row["Number_variant"]
            if _w > 0:
                if _var1 > _var2:
                    _temp_edge.add((_var1, _var2, _w))
                else:
                    _temp_edge.add((_var2, _var1, _w))

    for e in _temp_edge:
        _new_g.add_edge(e[0], e[1], weight = e[2])
    
    weight_arch=[]
    
    if  len(_new_g.es)==0:
        weight_arch=[]
    else:
        weight_arch=_new_g.es["weight"]

    #fig, ax = plt.subplots(dpi=500)
    ig.plot(_new_g,
        #target=ax,
        f"{path_saved}/plot_patient_connection/cluster_{cluster_index}.pdf",
        **{
        #"edge_width": _new_g.es["weight"],
        "edge_width":weight_arch,
        "vertex_color": "cyan",
        "bbox": (1280,720),
        "edge_curved": 0.1,
        #"vertex_label": [x.replace("_", "\n") for x in _new_g.vs["name"]],
        "vertex_label":  _new_g.vs["name"],
        "vertex_label_size": 4,
        #"edge_color": [get_color(x, _max_conn, _min_conn) for x in _new_g.es["weight"]],
        "edge_color": [get_color_comutated_cluster(x, _max_conn, _min_conn) for x in weight_arch],
        "edge_label": weight_arch,
        "edge_label_size": [5 + (n * 2) for n in weight_arch],
        "vertex_frame_width":0.05,
        "background":(0.3, 0.4, 0.5, 1)
        }
    )



# ONCOPLOT CON INFO QUALI:
def oncoplot(MAP_PATIENTS,MAP_VARIANTS,cluster_more_patient,path_saved):
    import os


    DF_IDENTITY=pd.DataFrame(columns=["sample", "vars_same_cluster", "vars_outside_cluster", "tot_vars", "cluster_patient"])
    DF_VARS=pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    DF_COUNTER_VARS=pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    DF_SIDE_BAR=pd.DataFrame(columns=["category", "data", "cluster_patient"])
    DF_TMB=pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    DF_MSI=pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    DF_MMR=pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    DF_CRS=pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    DF_SMOKER_STATUS = pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    #DF_SEX= pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])
    DF_ONCO_TREE=pd.DataFrame(columns=["sample", "category", "value", "cluster_patient"])

    temp_info_genes={}

    for patient,infos in MAP_PATIENTS.items():
        variant_inside=0
        variant_outside=0
        cluster_paz=infos["cluster"]
        variants=list(infos["variants"])
        smoker_status = infos["SMOKER"]
        #sex= infos["SEX"]
        onco_tree=infos["ONCOTREE_CODE"]
        #sample_class = infos["SAMPLE_CLASS"]
        tmb=infos["TMB.CLASS"]
        msi=infos["MSI.CLASS"]
        mmr=infos["MMR"]
        crs=infos["CRS"]
        for variant in variants:
            cons=MAP_VARIANTS[variant]["cons"]
            gene=variant.split("_")[0]
            cluster_variant=MAP_VARIANTS[variant]["cluster"]
            #VARS
            _new_record = pd.DataFrame([{'sample': patient, 'category':gene,'value': cons, 'cluster_patient': cluster_paz}])
            DF_VARS = pd.concat([DF_VARS, _new_record], ignore_index=True)

            #COUNTER_GENES
            temp_info_genes[cluster_paz]=temp_info_genes.get(cluster_paz, {})
            temp_info_genes[cluster_paz][gene]=temp_info_genes[cluster_paz].get(gene, 0)
            temp_info_genes[cluster_paz][gene]+=1



            if cluster_paz==cluster_variant:
                variant_inside+=1
            else:
                variant_outside+=1
        
        #IDENTITY
        tot_variant=variant_inside+variant_outside
        #tot_variant=len(infos["variants"])
        percent_inside=(variant_inside*100)/tot_variant
        percent_outside=(variant_outside*100)/tot_variant
        _new_record = pd.DataFrame([{'sample': patient, 'vars_same_cluster': percent_inside, 'vars_outside_cluster': percent_outside, "tot_vars": tot_variant,'cluster_patient': cluster_paz}])
        DF_IDENTITY = pd.concat([DF_IDENTITY, _new_record], ignore_index=True)
        
        #COUNTER_VARS
        _new_record = pd.DataFrame([{'sample': patient, 'category': "Variants", 'value': len(infos["variants"]), 'cluster_patient': cluster_paz}])
        DF_COUNTER_VARS = pd.concat([DF_COUNTER_VARS, _new_record], ignore_index=True)

        #TMB VALUE
        _new_record= pd.DataFrame([{'sample':patient,'category':"TMB.CLASS",'value':tmb,'cluster_patient':cluster_paz}])
        DF_TMB=pd.concat([DF_TMB,_new_record],ignore_index=True)


        #MSI VALUE
        _new_record= pd.DataFrame([{'sample':patient,'category':"MSI.CLASS",'value':msi,'cluster_patient':cluster_paz}])
        DF_MSI=pd.concat([DF_MSI,_new_record],ignore_index=True)
        DF_MSI.fillna("nan",inplace=True)

        #MMR VALUE
        _new_record= pd.DataFrame([{'sample':patient,'category':"MMR",'value':mmr,'cluster_patient':cluster_paz}])
        DF_MMR=pd.concat([DF_MMR,_new_record],ignore_index=True)
        DF_MMR.fillna("nan",inplace=True)

        #CRS VALUE
        _new_record= pd.DataFrame([{'sample':patient,'category':"CRS",'value':crs,'cluster_patient':cluster_paz}])
        DF_CRS=pd.concat([DF_CRS,_new_record],ignore_index=True)
        DF_CRS.fillna("nan",inplace=True)



        #SMOKING STATUS
        _new_record = pd.DataFrame([{'sample':patient,'category':"SMOKER",'value':smoker_status,'cluster_patient':cluster_paz}])
        DF_SMOKER_STATUS =pd.concat([DF_SMOKER_STATUS, _new_record],ignore_index=True)
        DF_SMOKER_STATUS.fillna("nan",inplace=True)

        #sex
        #_new_record = pd.DataFrame([{'sample':patient,'category':"SEX",'value':sex,'cluster_patient':cluster_paz}])
        #DF_SEX =pd.concat([DF_SEX, _new_record],ignore_index=True)
        #DF_SEX.fillna("nan",inplace=True)

        #ONCO TREE CODE
        _new_record = pd.DataFrame([{'sample':patient,'category':"ONCOTREE_CODE",'value':onco_tree,'cluster_patient':cluster_paz}])
        DF_ONCO_TREE =pd.concat([DF_ONCO_TREE, _new_record],ignore_index=True)
        DF_ONCO_TREE.fillna("nan",inplace=True)

    #media di adesione tra le mutazioni dei pazienti di un cluster e le mutazioni presenti nello stesso cluster
    media_cluster=DF_IDENTITY.groupby("cluster_patient")["vars_same_cluster"].mean()
    #df.groupby('cluster_patient')['vars_inside'].mean()

    if os.path.exists(f"{path_saved}/oncoplot"):
        os.makedirs(f"{path_saved}/oncoplot")

    with open(f"{path_saved}/oncoplot/mean_patient_variant_adhesion.csv","w") as f:
        f.write("Cluster\tMean_Adesion\n")
        for i in cluster_more_patient:
            f.write(str(i)+"\t"+str(media_cluster[i])+"\n")

    for _cluster_id, _v in temp_info_genes.items():
        for _gene, _value in _v.items():
            #SIDE_BAR
            _new_record = pd.DataFrame([{'category': _gene, 'data': _value, 'cluster_patient': _cluster_id}])
            DF_SIDE_BAR = pd.concat([DF_SIDE_BAR, _new_record], ignore_index=True)

        #generazione di 20 colori diversi da dare alle consequences
        consequence=DF_VARS["value"].unique()
        colors = [
            "orange", "orchid", "yellowgreen", "gold", "forestgreen", "lightseagreen",
            "aqua", "dodgerblue", "navy", "slategray", "pink", "olive",
            "deepskyblue", "lavender", "darkcyan", "darkgoldenrod", "seagreen", "darkorange",
            "fuchsia", "springgreen", "deeppink", "coral", "tan", "peru","mediumpurple","crimson","powderblue","steelblue","lightpink","crimson","c","skyblue",
            "azure","mediumslateblue","wheat","ivory","black","gold","violet","palevioletred","indianred","dimgray","grey","snow","tomato","limegreen"
        ]
        mut_mapping={}
        for termine, colore in zip(consequence, colors):
            mut_mapping[termine]=colore

        mut_mapping['Absent']={'facecolor': 'grey','alpha':0.2}


    def plot_comut(_CLUSTER_ID):
        bar_mapping = {'vars_same_cluster': 'turquoise', 'vars_outside_cluster': 'royalblue'}
        tmb_mapping={
            "High":"navy",
            "Medium":"steelblue",
            "Low":"deeppink",
            "NI":"grey"
            
        }
        os_mapping={
            "1:DECEASED":"darkviolet",
            "0:LIVING":"crimson",
        }
        smoker_mapping={
            "Yes":"dodgerblue",
            "No": "cadetblue",
            "nan":"grey"

        }
        crs_mapping={
            "3":"pink",
            "1 - 2":"blue",
            "nan":"grey",
            "NI": "dimgray"
        },
        mmr_mapping={
            "Proficient":"mediumpurple",
            "nan":"grey",
            "Deficient":"fuchsia",
            "NI":"dimgray"
        },
        msi_mapping={
            "Stable":"forestgreen",
            "Instable":"dodgerblue",
            "NI":"grey"

        },


        oncotree_code_mapping ={
            "UTERUS":"pink",
            "OVARY":"magenta",
            "NSCLC":"yellow",
            "BOWEL":"chocolate"
        }

        toy_comut = comut.CoMut()

        toy_comut.add_categorical_data(DF_VARS[DF_VARS["cluster_patient"]==_CLUSTER_ID], name = 'Consequence', mapping=mut_mapping, category_order=list(DF_SIDE_BAR[DF_SIDE_BAR["cluster_patient"]==_CLUSTER_ID].sort_values(by=["data"], ascending=True)["category"]))
        toy_comut.add_bar_data(DF_IDENTITY[DF_IDENTITY["cluster_patient"]==_CLUSTER_ID][["sample", "vars_same_cluster", "vars_outside_cluster"]], mapping=bar_mapping, stacked = True, name="variants in cluster", ylabel="identity")
        #toy_comut.add_continuous_data(DF_COUNTER_VARS[DF_COUNTER_VARS["cluster_patient"]==_CLUSTER_ID], name = 'Variants',mapping = 'Purples')
        toy_comut.add_side_bar_data(DF_SIDE_BAR[DF_SIDE_BAR["cluster_patient"]==_CLUSTER_ID], paired_name = 'Consequence',name="bar_vars")
        toy_comut.add_categorical_data(DF_TMB[DF_TMB["cluster_patient"]==_CLUSTER_ID],name="TMB",mapping=tmb_mapping)

        #toy_comut.add_categorical_data(DF_MSI[DF_MSI["cluster_patient"]==_CLUSTER_ID],name="MSI",mapping=msi_mapping)
        #toy_comut.add_categorical_data(DF_MMR[DF_MMR["cluster_patient"]==_CLUSTER_ID],name="MMR",mapping=mmr_mapping)
        #toy_comut.add_categorical_data(DF_CRS[DF_CRS["cluster_patient"]==_CLUSTER_ID],name="CRS",mapping=crs_mapping)

        toy_comut.add_categorical_data(DF_SMOKER_STATUS[DF_SMOKER_STATUS["cluster_patient"]==_CLUSTER_ID],name="SMOKER_STATUS",mapping=smoker_mapping)
        #toy_comut.add_categorical_data(DF_SEX[DF_SEX["cluster_patient"]==_CLUSTER_ID],name="SEX",mapping=sex_mapping)
        toy_comut.add_categorical_data(DF_ONCO_TREE[DF_ONCO_TREE["cluster_patient"]==_CLUSTER_ID],name="ONCOTREE_CODE",mapping=oncotree_code_mapping)
        toy_comut.plot_comut(figsize = (16, 15), x_padding = 0.1, y_padding = 0.1, hspace=0.01, wspace=0.01)
        toy_comut.add_unified_legend(axis_name="bar_vars")
        #toy_comut.add_unified_legend(axis_name="")
        toy_comut.figure.savefig(f'{path_saved}/oncoplot/comut_.{_CLUSTER_ID}.pdf', bbox_inches = 'tight', dpi = 300)

    if sys.platform.startswith('win') or sys.platform.startswith("linux"):
        with Pool() as p:
            mod_results = p.map(plot_comut, cluster_more_patient)
    else:
        for _i in cluster_more_patient:
            plot_comut(_i)



#prendi paziente 0 del cluster A
#calcolo indice di somiglianza con tutti i glia altri pazienti contenuti nel cluster B (tranne se stesso)
#i valori vengono aggiunti ad un array temporaneo
#passo al paziente successivo, una volta terminato il cliclo avrò un array con i valori di somiglianza
#return

def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(set(list1)) + len(set(list2))) - intersection
    return float(intersection) / union

def somiglianze_cluster(cluster_a, cluster_b,MAP_CLUSTER,MAP_PATIENTS):
    list_patients_cluster_a = [x for x in MAP_CLUSTER[cluster_a] if x[1] == "PATIENT"]
    list_patients_cluster_b = [x for x in MAP_CLUSTER[cluster_b] if x[1] == "PATIENT"]
    values = []

    #IN
    if cluster_a == cluster_b:
        for paz in list(combinations(list_patients_cluster_a, 2)):
            variant_paz_a = MAP_PATIENTS[paz[0][0]]["variants"]
            variant_paz_b = MAP_PATIENTS[paz[1][0]]["variants"]
            values.append(jaccard_similarity(variant_paz_a, variant_paz_b))
    #OUT
    else:
        for paz_a in list_patients_cluster_a:
            for paz_b in list_patients_cluster_b:
                variant_paz_a = MAP_PATIENTS[paz_a[0]]["variants"]
                variant_paz_b = MAP_PATIENTS[paz_b[0]]["variants"]
                values.append(jaccard_similarity(variant_paz_a, variant_paz_b))

    return values


def box_plot_similitudine_one_cluster(cluster_more_patient,path_saved):
    import tap

    cluster_base = 0
    values = somiglianze_cluster(cluster_base, cluster_base)
    values_out = []
    for i in cluster_more_patient:
        if i == cluster_base:
            continue
        values_out+=somiglianze_cluster(cluster_base,i)
    df_temp = pd.DataFrame({
        "cluster_base": [cluster_base]*len(values) + [cluster_base]*len(values_out),
        "category": [f"C_{cluster_base}_IN"]*len(values) + [f"C_{cluster_base}_OUT"]*len(values_out),
        "values": values + values_out
    })
    print(values)
    tap.plot_stats(df_temp, "category", "values",filename=f"{path_saved}/adhesion_cluster_{cluster_base}.png")


def box_plot_similitudine_all_clusters(cluster_more_patient,):
    import tap
    for index_cluster in cluster_more_patient:
        #IN
        df_values = pd.DataFrame(columns=["Cluster", "Category", "Value"])
        values_in = somiglianze_cluster(index_cluster, index_cluster)
        for value in values_in:
            df_values.loc[len(df_values.index)] = [index_cluster, f"C_{index_cluster}_IN", value]
        #OUT
        for index_cluster_out in cluster_more_patient:
            if index_cluster == index_cluster_out:
                continue
            values_out = somiglianze_cluster(index_cluster, index_cluster_out)
            for value in values_out:
                df_values.loc[len(df_values.index)] = [index_cluster, f"C_{index_cluster}_OUT", value]

        tap.plot_stats(df_values, "Category", "Value")



'''def main():# NOME DELLO STUDIO
    import os 
    name_study="./NEW_DATA_NETWORK/Data/LUNG"
    path_results_first="./NEW_DATA_NETWORK/Data/LUNG/output_egfr_categories"

    path_results_second="./NEW_DATA_NETWORK/Data/LUNG/output_egfr_categories/output_comutation"

    if not os.path.exists(path_results_first):
        os.makedirs(path_results_first)

    if not os.path.exists(path_results_second):
        os.makedirs(path_results_second)


    
    path_mutational=f"{name_study}/data_mutations_extended.txt"
    path_sample=f"{name_study}/data_clinical_sample.txt"
    path_clinical=f"{name_study}/data_clinical_patient.txt"

    data_mutational,data_sample,data_patient=read_file(path_mutational,path_sample,path_clinical)

    #summury_dataset(data_sample,"ONCOTREE_CODE",name_study)
    summury_dataset(data_patient,"GENDER",path_results_first)
    #summury_dataset(data_sample,"SAMPLE_CLASS",name_study)
    summury_dataset(data_patient,"SMOKER",path_results_first)
    summury_dataset(data_patient,"SMOKING.STATUS",path_results_first)
    summury_dataset(data_patient,"PD.L1.CATEGORIES",path_results_first)


    EGFR=True
    TP53=False

    map_patients,map_variants=created_maps(data_mutational,EGFR,TP53)
    graph = graph_creation(map_patients,map_variants)
    best_seed= selected_seed(graph)
    clustering_graph_first = leiden_clustering(graph,best_seed)

    #MAP CLUSTER
    map_cluster_graph_first=map_cluster(graph,clustering_graph_first)

    #add data sample
    enriched_sample_data(data_sample,map_cluster_graph_first,map_patients,map_variants)

    #add patient data
    enriched_patient_data(data_patient,map_patients)

    #enriched sample data
    creation_cluster_clinical_data(map_patients,path_results_first)

    #division of patient, variant and cluster
    cluster_more_patient,cluster_one_patient= cluster_division(map_cluster_graph_first)

    #
    variant_patient_connection_count=variant_conncection_patient(clustering_graph_first)

    patient_variant_connection_count=patient_connection_variant(clustering_graph_first)

    file_connection_variant(map_cluster_graph_first,map_patients,path_results_first)

    file_connection_patient(map_cluster_graph_first,map_variants,path_results_first)

    #cREAZIONE DI UN FILE "VARIANT_CLUSTER" CHE RIASSUMA CLUSTER+MUTAZIONE+GENE
    file_summary_variant_gene(path_results_first,map_cluster_graph_first,map_variants)

    #add size nodi primo grafo
    add_size_node(graph,variant_patient_connection_count)

    #conteggio numero totale di mutazioni per gene
    gene_total_count=count_gene(graph)
    #conteggio valore assoluto e percentuale di presenza di ciascun gene nei diversi cluster
    map_cluster_gene_count,map_cluster_gene_percent=gene_cluster_total_and_percent(map_cluster_graph_first,gene_total_count)
    #file salvataggio % dei geni
    file_gene_percentage(map_cluster_gene_percent,path_results_first)
    #creazione di un file con la lista dei geni per ogni cluster
    genes_single_cluster(map_cluster_gene_count,path_results_first)

    #aggiunta del cluster i nodi
    cluster_noded_attributes(graph,map_patients,map_variants)
    #aggiunta di sost_amm ai nodi
    sost_amm_attributes(map_variants,graph)

    #plot delle comutazioni a partire dal grafo iniziale
    connection_mutation_df = pd.read_csv(f"{path_results_first}/connection_patient.csv", sep="\t")
    for cluster_index in cluster_more_patient:
        plot_distance_comutated_cluster_variants(clustering_graph_first,cluster_index,connection_mutation_df,path_results_first)
    if sys.platform.startswith("win") or sys.platform.startswith("linux"):
        with Pool() as p:
            p.map(plot_distance_comutated_cluster_variants, cluster_more_patient)
    else:
        for i in cluster_more_patient:
            plot_distance_comutated_cluster_variants(i)




    
    #********COMUTATED ALL MUTATIONS************+
    graph_all_comutations,lista_edge_only_varianti=graph_comutated_variants(map_variants)
    clustering_all_comutated = leiden_comutated_variants(graph_all_comutations)
    file_comutated_variant_cytoscape(lista_edge_only_varianti,2,map_variants,path_results_second)
    #MAP CLUSTER MUTATION
    map_cluster_graph_mutation = map_cluster(graph_all_comutations,clustering_all_comutated)


    


if __name__=="__main__":
    main()'''


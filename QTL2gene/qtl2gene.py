from SPARQLWrapper import SPARQLWrapper, JSON

import pandas as pd        
import pickle, hashlib   
import scipy.stats as stats
from statutils.multi_comparison import p_adjust

class SEARCH:

    def __init__(self, url_pbg):
        #define cache directory
        self.cache = "cache/"
        #define url
        self.url_pbg = url_pbg
        #define sparql
        self.sparql_pbg = SPARQLWrapper(self.url_pbg)
        self.sparql_pbg.setReturnFormat(JSON)        
        
    def cache_name(self, method, parameters) :
        key = method+"_"+hashlib.md5(pickle.dumps(parameters)).hexdigest()
        return(key)
    
    
    def get_location(self, id):
        filename = self.cache + self.cache_name("get_location", id)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/gene_location.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_pbg.setQuery(query % id)
            # JSON example
            response = self.sparql_pbg.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    result.append([
                    item["gene_id"]["value"],
                    item["chromosome"]["value"],
                    item["begin_ref"]["value"],
                    item["begin_pos"]["value"],
                    item["end_ref"]["value"],
                    item["end_pos"]["value"]])
                df = pd.DataFrame(result)  
                df.columns = ["gene_id", "chromosome", "begin_ref", "begin_pos", "end_ref", "end_pos" ]
                df = df.set_index("gene_id")
                df["begin_pos"] = pd.to_numeric(df["begin_pos"])
                df["end_pos"] = pd.to_numeric(df["end_pos"])
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()   
        
    def compute_interval(self, g1, g2):  
        locations = pd.concat([self.get_location(g1), self.get_location(g2)])
        display(locations[["location"]])
        if(len(locations.index)!=2) :
            print("unexpected number of rows in locations:",len(locations.index))
        elif(locations.iloc[0]['end_pos']>locations.iloc[1]['begin_pos']) & (g1!=g2) :
            print("unexpected order",locations.index[0],"and",locations.index[1])
        else :
            result = []
            if locations.iloc[0]["end_pos"]>locations.iloc[0]["begin_pos"] :
              result.append(["begin", locations.iloc[0]["end_ref"], locations.iloc[0]["end_pos"]])
            else :
              result.append(["begin", locations.iloc[0]["begin_ref"], locations.iloc[0]["begin_pos"]])
            if locations.iloc[1]["begin_pos"]<locations.iloc[1]["end_pos"] :
              result.append(["end", locations.iloc[1]["begin_ref"], locations.iloc[1]["begin_pos"]])
            else :
              result.append(["end", locations.iloc[1]["end_ref"], locations.iloc[1]["end_pos"]])
            df = pd.DataFrame(result)
            df.columns = ["type", "ref", "pos" ]
            df = df.set_index("type")
            return df    
        
    def make_interval(self, ref, start, end):
        result = []
        result.append(["begin",ref,start])
        result.append(["end",ref,end])
        df = pd.DataFrame(result)
        df.columns = ["type", "ref", "pos" ]
        df = df.set_index("type")
        return df    
        
    def interval_genes(self, interval):
        filename = self.cache + self.cache_name("interval_genes", interval)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/interval_genes.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_pbg.setQuery(query % {"beginRef" : interval.loc["begin"]["ref"], "beginPos" : interval.loc["begin"]["pos"], "endRef" : interval.loc["end"]["ref"], "endPos" : interval.loc["end"]["pos"]})
            # JSON example
            response = self.sparql_pbg.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = []
                    row.append(item["gene_id"]["value"])
                    row.append(item["chromosome"]["value"])
                    row.append(item["begin_pos"]["value"])
                    row.append(item["end_pos"]["value"])
                    result.append(row)
                df = pd.DataFrame(result)  
                df.columns = ["gene_id", "chromosome", "start", "end"]
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()  
            
    def go_genes(self, graphEnsembl, graphUniprot, go):
        filename = self.cache + self.cache_name("go_genes", [graphEnsembl, graphUniprot, go])
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/go_genes.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_pbg.setQuery(query % {"graphEnsembl" : graphEnsembl, "graphUniprot" : graphUniprot, "go" : go})
            # JSON example
            response = self.sparql_pbg.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = []
                    row.append(item["gene_count"]["value"])
                    row.append(item["gene_with_go_count"]["value"])
                    result.append(row)
                df = pd.DataFrame(result)  
                df.columns = ["gene_count", "gene_with_go_count"]
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()              
            
    def gene_goterms(self, id):  
        filename = self.cache + self.cache_name("gene_goterms", id)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/gene_goterm.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_pbg.setQuery(query % id)
            # JSON example
            response = self.sparql_pbg.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = []
                    row.append(item["gene_id"]["value"])
                    row.append(item["go_id"]["value"])
                    row.append(item["go_term"]["value"])
                    row.append(item["go_cat"]["value"])
                    row.append(item["graph_ensembl"]["value"])
                    row.append(item["graph_uniprot"]["value"])
                    result.append(row)
                df = pd.DataFrame(result)  
                df.columns = ["gene_id", "go_id", "go_term", "go_cat", "graph_ensembl", "graph_uniprot"]
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()  
            
    def genes_goterms(self, ids):
        list = []
        for id in ids:
            list.append(self.gene_goterms(id))
        return pd.concat(list).reset_index(drop=True)
    
    def get_go_numbers(self, goterms, genes):
        #construct the number of genes with/without goterm
        graphs = list(goterms.groupby(["graph_ensembl","graph_uniprot"]).indices.keys())
        golist = goterms["go_id"].unique()
        #construct df
        df = pd.DataFrame(goterms.groupby("go_id").size(), columns=["interval_genes_annotated"])
        #add gene numbers
        df["interval_genes_not_annotated"] = len(genes.index) - df["interval_genes_annotated"]
        df["outside_genes_annotated"] = 0
        df["outside_genes_not_annotated"] = 0
        df["total_genes"] = 0
        for go in golist:
            for graph in graphs:    
              result = self.go_genes(graph[0], graph[1], go)
              df.loc[go,"outside_genes_annotated"] = df.loc[go,"outside_genes_annotated"] + int(result.loc[0,"gene_with_go_count"])  
              df.loc[go,"total_genes"] = df.loc[go,"total_genes"] + int(result.loc[0,"gene_count"])  
            df.loc[go,"outside_genes_annotated"] = df.loc[go,"outside_genes_annotated"] - df.loc[go,"interval_genes_annotated"]
            df.loc[go,"outside_genes_not_annotated"] = df.loc[go,"total_genes"] - df.loc[go,"outside_genes_annotated"] - df.loc[go,"interval_genes_annotated"] - df.loc[go,"interval_genes_not_annotated"]
        #do fisher tests
        for go in golist:
            m = [[df.loc[go,"interval_genes_annotated"], df.loc[go,"outside_genes_annotated"]],[df.loc[go,"interval_genes_not_annotated"], df.loc[go,"outside_genes_not_annotated"]]]
            df.loc[go,"p_less"] = stats.fisher_exact(m, alternative="less")[1]
            df.loc[go,"p_greater"] = stats.fisher_exact(m, alternative="greater")[1]
        df["p_adjusted"] = p_adjust(df["p_greater"], method="BH")    
        return df
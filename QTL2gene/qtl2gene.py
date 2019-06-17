from SPARQLWrapper import SPARQLWrapper, JSON

import pandas as pd    
    
import pickle, hashlib   

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
                    item["location"]["value"],
                    item["begin_ref"]["value"],
                    item["begin_pos"]["value"],
                    item["end_ref"]["value"],
                    item["end_pos"]["value"]])
                df = pd.DataFrame(result)  
                df.columns = ["gene_id", "location", "begin_ref", "begin_pos", "end_ref", "end_pos" ]
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
                    row.append(item["location"]["value"])
                    result.append(row)
                df = pd.DataFrame(result)  
                df.columns = ["gene_id", "location"]
                df = df.set_index("gene_id")
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()      
        

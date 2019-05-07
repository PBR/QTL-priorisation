from SPARQLWrapper import SPARQLWrapper, JSON

import pandas as pd
        
class SEARCH:

    def __init__(self, url_pbg, url_oma, url_uniprot):
        #define url
        self.url_pbg = url_pbg
        self.url_oma = url_oma
        self.url_uniprot = url_uniprot
        #define sparql
        self.sparql_pbg = SPARQLWrapper(self.url_pbg)
        self.sparql_pbg.setReturnFormat(JSON)
        self.sparql_oma = SPARQLWrapper(self.url_oma)
        self.sparql_oma.setReturnFormat(JSON)
        self.sparql_uniprot = SPARQLWrapper(self.url_uniprot)
        self.sparql_uniprot.setReturnFormat(JSON)
        
    def url(self):
        print(self.url)
        
    def get_location(self, id):
        file = open("qtlsearch/gene_location.sparql", "r") 
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
        file = open("qtlsearch/interval_genes.sparql", "r") 
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
            return df 
        else:
            return pd.DataFrame()      
        
        
    def get_parent_groups(self, protein):
        file = open("qtlsearch/parent_groups.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql_oma.setQuery(query % {"protein" : "\""+protein+"\""})
        # JSON example
        response = self.sparql_oma.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                result.append([
                item["level"]["value"],
                item["group"]["value"],
                item["type"]["value"],
                item["protein"]["value"]])
            df = pd.DataFrame(result)  
            df.columns = ["level", "group", "type", "protein" ]
            df.drop_duplicates(subset ="group", keep = "first", inplace = True) 
            df = df.set_index("group")
            df["level"] = pd.to_numeric(df["level"])
            return df 
        else:
            return pd.DataFrame()  
        
    #get only group paths that end with a uniprot annotated protein
    def get_child_groups(self, parent):
        file = open("qtlsearch/child_groups.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql_oma.setQuery(query % {"parent" : "<"+parent+">"})
        # JSON example
        response = self.sparql_oma.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                row = [
                  item["group"]["value"],
                  item["type"]["value"]  
                ]
                if "parent" in item.keys() :
                    row.append(item["parent"]["value"])
                else:
                    row.append(None)
                if "parent_type" in item.keys() :
                    row.append(item["parent_type"]["value"])
                else:
                    row.append(None) 
                if "label" in item.keys() :
                    row.append(item["label"]["value"])
                else:
                    row.append(None)
                if "parent_label" in item.keys():
                    row.append(item["parent_label"]["value"])
                else:
                    row.append(None)    
                result.append(row)                
            df = pd.DataFrame(result)  
            df.columns = ["group", "type", "parent", "parent_type", "label", "parent_label" ]
            df.drop_duplicates(subset ="group", keep = "first", inplace = True) 
            df = df.set_index("group")
            return df 
        else:
            return pd.DataFrame()    
        
    #get uniprot annotated proteins and their group
    def get_child_proteins_uniprot(self, parent):
        file = open("qtlsearch/child_proteins_uniprot.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql_oma.setQuery(query % {"parent" : "<"+parent+">"})
        # JSON example
        response = self.sparql_oma.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                row = [
                  item["group"]["value"],                  
                  item["protein"]["value"]
                ] 
                if "protein_uniprot" in item.keys() :
                    row.append(item["protein_uniprot"]["value"])
                else:
                    row.append(None)
                if "group_label" in item.keys() :
                    row.append(item["group_label"]["value"])
                else:
                    row.append(None)
                result.append(row)                
            df = pd.DataFrame(result)  
            df.columns = ["group", "protein", "uniprot", "label" ]
            df.drop_duplicates(subset ="uniprot", keep = "first", inplace = True) 
            df = df.set_index("uniprot")
            return df 
        else:
            return pd.DataFrame()    
        
    #get proteins and their group
    def get_child_proteins(self, parent):
        file = open("qtlsearch/child_proteins.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql_oma.setQuery(query % {"parent" : "<"+parent+">"})
        # JSON example
        response = self.sparql_oma.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                row = [
                  item["group"]["value"],                  
                  item["protein"]["value"]
                ] 
                if "group_label" in item.keys() :
                    row.append(item["group_label"]["value"])
                else:
                    row.append(None)
                result.append(row)                
            df = pd.DataFrame(result)  
            df.columns = ["group", "protein", "label" ]
            df.drop_duplicates(subset ="protein", keep = "first", inplace = True) 
            df = df.set_index("protein")
            return df 
        else:
            return pd.DataFrame()    
        
        
    #get child annotations
    def get_child_annotations(self, go_annotation):
        file = open("qtlsearch/child_annotations.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql_pbg.setQuery(query % {"annotation" : go_annotation})
        # JSON example
        response = self.sparql_pbg.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                row = [
                  item["go_annotation"]["value"],
                  item["label"]["value"]
                ]                
                result.append(row)                
            df = pd.DataFrame(result)  
            df.columns = ["go_annotation", "label" ]
            df.drop_duplicates(subset ="go_annotation", keep = "first", inplace = True) 
            df = df.set_index("go_annotation")
            return df 
        else:
            return pd.DataFrame()         
        
    #check uniprot for proteins and GO annotations
    def check_uniprot_annotations(self, uniprot_proteins, go_annotations):
        file = open("qtlsearch/check_uniprot_annotations.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql_uniprot.setQuery(query % {"proteins" : "<"+">,<".join(uniprot_proteins)+">", "annotations": "<"+">,<".join(go_annotations)+">"})
        # JSON example
        response = self.sparql_uniprot.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                row = [
                  item["uniprot"]["value"],
                  item["reviewed"]["value"]
                ]                
                result.append(row)                
            df = pd.DataFrame(result)  
            df.columns = ["uniprot", "reviewed" ]
            #df["reviewed"] = df["reviewed"].astype("bool")
            df.drop_duplicates(subset ="uniprot", keep = "first", inplace = True) 
            df = df.set_index("uniprot")
            return df 
        else:
            return pd.DataFrame()
        
    
    

    
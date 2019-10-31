import yaml

# This class implements a very simple ASCII parameter file database.
# Parameters are accessed as if they were parts of a class/object hierarchy.
   
class DB:
    
    def __init__(self):
        self.yaml_doc = None
    
    def get_entry_from_yaml_doc(self, entry_keys):       
        # split the key path into individual params
        entry_keys = entry_keys.split('/') 
        entry = self.yaml_doc
        for key in entry_keys:
            #print('key=',end='')
            #print(key)
            entry = entry[key]
        
        return entry
    
    def load_file(self,file):
        fh = open(file, 'r')
        self.yaml_doc = yaml.load(fh) 
        fh.close()

    def param_defined(self,param_path):
        param_path = param_path.split('/')         
        entry = self.yaml_doc
        for key in param_path:            
            try:
                entry = entry[key]
            except KeyError:
                return False
        
        return True
                
    def DB_GET_STRING(self,param_path):
        # split the key path into individual params
        
        db_string = self.get_entry_from_yaml_doc(param_path)        
        return db_string


    def DB_GET_DOUBLE(self,param_path):
        db_double = self.get_entry_from_yaml_doc(param_path)        
        return db_double


    def DB_GET_INTEGER(self,param_path):
        db_int = self.get_entry_from_yaml_doc(param_path)        
        return db_int


    def DB_GET_BOOL(self,param_path):
        db_bool = self.get_entry_from_yaml_doc(param_path)        
        return db_bool

    def DB_GET_STRING_LIST(self,param_path):
        db_string = self.get_entry_from_yaml_doc(param_path)   
        print(db_string)
        db_string = db_string.split(',')
        db_string_list = []
        for s in db_string:
            db_string_list.append(s.strip())
        
        return db_string_list



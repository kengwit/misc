import yaml
from elements.element_factory import ElementFactory

def strip_line(line):
    return [x.strip() for x in line.split(',')]

def MeshReadYAML(mec,filename):            
    fh = open(filename, 'r')
    doc = yaml.load(fh) 
    fh.close()
    
    mec.ndim = doc['mesh']['ndim']
    nodes = doc['node']['data']
    nnodes = len(nodes)
    mec.nnodes = nnodes
    for i in range(0,mec.nnodes):
        data = nodes[i]
        node_id = data[0]
        coords  = data[1:]
        mec.X_coords[node_id] = coords
        mec.x_coords[node_id] = coords
        mec.U_tau[node_id] = [0]*len(coords)
        mec.U_t[node_id] = [0]*len(coords)
        mec.dU[node_id] = [0]*len(coords)
        
    fe_data = doc['element']['data']
    fe_type = doc['element']['type'].upper()
    ndof_per_node = doc['element']['ndof_per_node']
    nelems = len(fe_data)
    mec.nelems = nelems
    for i in range(0,nelems):
        fe_id   = int(fe_data[i][0])
        fe_conn = fe_data[i][1:]
        
        el      = ElementFactory(fe_type)           
        el.id   = fe_id
        el.conn = fe_conn
        
        mec.ElementMap[fe_id] = el       
        
    mec.ndofs   = nnodes*ndof_per_node
    #mesh.U_tau   = np.zeros((nnodes,ndof_per_node))
    #mesh.U_t     = np.zeros((nnodes,ndof_per_node))
    #mesh.dU      = np.zeros((nnodes,ndof_per_node))
    #mesh.fext    = np.zeros((nnodes,ndof_per_node))
    #mesh.fint    = np.zeros((nnodes,ndof_per_node))
    #mesh.frct    = np.zeros((nnodes,ndof_per_node))
     
    


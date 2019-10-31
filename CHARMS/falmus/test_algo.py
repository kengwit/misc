from algo_elasticity import ALGO_ELASTICITY
from db import DB
from mgr import MGR



jobname = 'Mandible'

db = DB()
db.load_file(jobname + '.yaml')

mgr = MGR(db)

a = ALGO_ELASTICITY(jobname, mgr);
a.setup()


#for gcell in a.mgr().mesh_mgr().gmesh('m1')._gsubmeshes[0]._gcell_groups[0]._gcells:
#    print('GCELL ID = ' + str(gcell.id()) + ', GCELL TYPE = ' + gcell.type_name() + ', conn = ',end='')
#    conn = gcell._conn
#    #print(conn.fen(0))
#    for i in range(0, conn.nfens()):
#        print(str(conn.fen(i).id()) + ' ',end='')
#        
#    print('')


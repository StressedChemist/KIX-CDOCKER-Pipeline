import os
if os.getenv('CHARMM_LIB_DIR') == None:
    # set path to CHARMM_LIB_DIR
    pyCHARMM_LIB = 'CHARMM_LIB_DIR=/home/boloyede/charmm/pycharmm/lib'
    os.environ['CHARMM_LIB_DIR'] = pyCHARMM_LIB
import pycharmm
import pycharmm.lib as lib
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.settings as settings
from pycharmm.cdocker import Rigid_CDOCKER
home = os.getcwd()

settings.set_bomb_level(-1)
owl = settings.set_warn_level(-5)
read.rtf('"../toppar/top_all36_prot.rtf"')
read.rtf('"../toppar/top_all36_cgenff.rtf"', append = True)
read.prm('"../toppar/par_all36m_prot.prm"', flex = True)
read.prm('"../toppar/par_all36_cgenff.prm"', append = True, flex = True)
settings.set_bomb_level(0)
lingo.charmm_script("stream ligandrtf")
settings.set_warn_level(owl)

## File name and pathway
ligPDB = "mymol.pdb"
confDir = "conformer/"
receptorPDB = "2agh_prob.pdb"
receptorPSF = "2agh_prob.psf"


## Rigid CDOCKER standard docking protocol
clusterResult, dockResult = Rigid_CDOCKER(xcen =  -10.29894514001519, ycen = 1.0328956229548085, zcen = 4.694294212106652,
                                        flag_suppress_print = False,
                                        maxlen = 15, ligPDB = ligPDB, receptorPDB = receptorPDB,
                                        probeFile = "../toppar/fftdock_c36prot_cgenff_probes.txt",
                                        receptorPSF = receptorPSF, confDir = confDir, flag_grid = False,
                                        flag_fast_placement = False, flag_delete_conformer = False, numPlace = 5)


print(clusterResult)
print(dockResult)

if not os.path.exists('conformer'): os.mkdir('conformer')
for r in range(1,21):
    #os.system(f'obrotamer "output.pdb" > conformer/{r}.pdb')
    os.system(f'obrotamer "mymol.pdb" > conformer/{r}.pdb')

stop

import py3Dmol
view = py3Dmol.view(width=500, height=400)
view.setBackgroundColor('grey')
for i,struct in enumerate(['2agh_prob.pdb','2agh_prob.psf','dockresult/cluster/top_1.pdb',
                          'dockresult/top_ener/top_1.pdb']):
    view.addModel(open(struct,'r').read(),'pdb')
    if i > 1:
        view.setStyle({'model':i},{'stick':{'colorscheme':'greenCarbon'}})
    elif i == 0:
        view.setStyle({'model':i},{'stick':{'colorscheme':'yellowCarbon' }})
        view.addSurface(py3Dmol.VDW,{'opacity': 0.7},{'model':i})
    elif i == 1:
        view.setStyle({'model':i},{'cartoon':{'color': 'lightblue'}})
        view.addSurface(py3Dmol.SES,{'opacity': 0.8,'color': 'white'},{'model':i})
view.zoomTo()
view.show()
#view

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
home = os.getcwd()

settings.set_bomb_level(-1)
owl = settings.set_warn_level(-5)
read.rtf('"../toppar/top_all36_prot.rtf"')
read.rtf('"../toppar/top_all36_cgenff.rtf"', append = True)
read.prm('"../toppar/par_all36m_prot.prm"', flex = True)
read.prm('"../toppar/par_all36_cgenff.prm"', append = True, flex = True)
settings.set_bomb_level(0)
lingo.charmm_script("stream mymol.str")
#lingo.charmm_script("stream nolp2.str")
settings.set_warn_level(owl)

from pycharmm.cdocker import Rigid_CDOCKER
#help(Rigid_CDOCKER) #used for instructions

read.psf_card("yw_protein.psf", append = True)

exit()

if not os.path.exists('conformer'): os.mkdir('conformer')
for r in range(1,21):
    os.system(f'obrotamer "output.pdb" > conformer/{r}.pdb')
    #os.system(f'obrotamer "nolp.pdb" > conformer/{r}.pdb')

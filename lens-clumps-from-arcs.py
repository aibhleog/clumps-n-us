'''

The below only works if the arcs.dat file is provided.

'''
from fitting_ifu_spectra import *


# reading in data
# arcs = pd.read_csv('lensing/new-model/model files/arcs.dat',
# arcs = pd.read_csv('lensing/SGAS1110/pedram_source2_arcs1_SGAS1110.dat',
arcs = pd.read_csv('lensing/SGAS1110/pedram_source2ONLY_arcs1_SGAS1110.dat',
                   sep=r'\s+',comment='#',
                   index_col=None,names=['sID','ra','dec','s1','s2','s3','z','s4'])

# -- not necessary since adding the "comment" kwarg above
# t = arcs.sID.values
# mask = [tt[0]=='#' for tt in t]
# arcs = arcs.loc[np.logical_not(mask)].copy()



reghead = '''# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5'''


# making region file
with open('lensing/SGAS1110/source2-regions.reg','w') as f:
    print(reghead,file=f)
    for i in arcs.index.values:
        # if arcs.loc[i,'z'] == 5.043:
        # filler = f'circle({arcs.loc[i,'ra']},{arcs.loc[i,'dec']},{arcs.loc[i,'s1']}") # text='+'{'+str(arcs.loc[i,'sID'])+'}'
        filler = f'circle({arcs.loc[i,'ra']},{arcs.loc[i,'dec']},{arcs.loc[i,'s1']}") # text='+'{'+str(arcs.loc[i,'sID'])+'}'
        print(filler,file=f)
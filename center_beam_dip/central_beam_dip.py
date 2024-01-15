import os,mrcfile
from utils import*                      ;imp.reload(dsp)
from EDutils import utilities as ut     ;imp.reload(ut)
from EDutils import pets as pt          ;imp.reload(pt)
from EDutils import dials_utils as dials;imp.reload(dials)
# from EDutils import display as EDdisp   ;imp.reload(EDdisp)
from blochwave import bloch_pp as bl    ;imp.reload(bl)
from subprocess import check_output
from scipy.signal import find_peaks

plt.close('all')
# edly_path='/home/tarik/Documents/git/ccp4/src/edly'
mol_path='../biotin_test'

#options : Solve(S) p(plot rocking) n(new)
opts = 'a'
name = 'central_beam_fine'
thicks = [20,100,300,500,1000]  #Angstrom
Smax   = 0.01
Nmax   = 13
n_scans  = 3
n_frames = 140
a=np.deg2rad(140)
# cif_file='dat/1ejg.pdb'
cif_file='biotin.cif'


### paths
rocks_path=os.path.join(mol_path,'rocks')
figs_path= os.path.join(rocks_path,'figs')
rock_path=lambda i: os.path.join(rocks_path,'%s_%d' %(name,i))
rock_file=lambda i: os.path.join(rock_path(i),'rock_.pkl')
fig_file =lambda i: os.path.join(figs_path,'%s_%d.png' %(name,i))
if not os.path.exists(figs_path):
    check_output('mkdir %s' %figs_path,shell=True)



if 'n' in opts:check_output('rm -rf %s/*' %rocks_path,shell=True)
for i in range(n_scans):
    ## Simulate
    if 'S' in opts or not os.path.exists(rock_path(i)):
        u_0 = np.random.rand(3)
        u_1 = np.cross(u_0,[0,1,0])
        u_end = np.cos(a)*u_0+np.sin(a)*u_1
        u_end/=np.linalg.norm(u_end)
        uvw   = ut.get_uvw_cont(u_0,u_end,n_frames,show=0)#;plt.show()
        rock  = bl.Bloch_cont(path=rock_path(i),uvw=uvw,tag='',
            Sargs = {'cif_file':os.path.join(rocks_path,cif_file),'thick':thicks[0],
                'Smax':Smax,'Nmax':Nmax,'solve':1,'keV':200},
            )
        rock.set_beams_vs_thickness(thicks=thicks)
        rock.save()

    ### plot
    if 'p' in opts:
        rock=pd.read_pickle(rock_file(i))
        rock.change_path(rock_path(i))
        rock.plot_rocking(refl=[str((0,0,0))],zs=thicks,x='Frame',
            cmap='viridis',opt='sc',name=fig_file(i))

    if 'a' in opts:
        rock=pd.read_pickle(rock_file(i))
        rock.change_path(rock_path(i))
        full_rocks=rock.get_rocking(refl=[str((0,0,0))],zs=500)
        y=full_rocks[1][str((0,0,0))].T[0]

        idx = find_peaks(1-y,height=0.2)[0]
        b=rock.load(0)
        Arec=np.linalg.inv(b.lat_vec.T)
        uvw=rock.uvw[idx]
        hkl=Arec.dot(uvw.T)
        print(hkl)
# print(rock_file(0))
# rock=ut.load_pkl(rock_file(0))
# rock.rock_frames=[1,n_frames]
# rock.path=rocks_path
# rock.save()

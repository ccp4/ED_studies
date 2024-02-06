from utils import*;imp.reload(dsp)
import gemmi
from EDutils import utilities as ut 

def get_hkls(all_beams,N,h,index=False):
    hkls = [hs.index[:N].values.tolist() for hs in all_beams]
    _hkls=[]
    for hs in hkls : 
        if h not in hs: hs[-1]=h
        if not index:
            hs=np.array([eval(h) for h in hs])
        _hkls.append(hs)
    return _hkls
    
def get_uvw(u,npts=20,osc = 0.2):
    ez = [0,0,1]
    rot_axis = np.cross(ez,u)

    u1 = ut.rotation_matrix(rot_axis,-osc,deg=True).dot(u)
    u2 = ut.rotation_matrix(rot_axis, osc,deg=True).dot(u)
    uvw=ut.get_uvw_cont(u1,u2,nframes=npts)
    return uvw
    
def gemmi_sf(cif_file:str='',dmin=2):
    st = gemmi.read_structure(cif_file)
    if len(st):st=st[0]
    dc = gemmi.DensityCalculatorX()
    dc.d_min = dmin
    dc.addends.subtract_z()
    dc.set_grid_cell_and_spacegroup(st)
    dc.set_refmac_compatible_blur(st)
    dc.put_model_density_on_grid(st)
    grid = gemmi.transform_map_to_f_phi(dc.grid); print('sf : ',grid.shape)
    Nmax = (np.array(grid.shape)//2-1).min()
    Fhkl = np.zeros((2*Nmax+1,)*3,dtype=complex)
    hklF = np.meshgrid(*[np.arange(-Nmax,Nmax+1)]*3)
    hkls = np.array([i.flatten() for i in hklF]).T
    idxs = np.array([i.flatten() for i in np.meshgrid(*[np.arange(2*Nmax+1)]*3)]).T
    #TODO : multiple size grid
    # nu2,nv2,nw2 = np.array(grid.shape)//2-1
    # Fhkl = np.zeros((2*nu2+1,2*nv2+1,2*nw2+1),dtype=complex)
    # hkls=np.array([i.flatten() for i in np.meshgrid(range(-nu2,nu2+1),range(-nv2,nv2+1),range(-nw2,nw2+1))]).T
    # idxs=np.array([i.flatten() for i in np.meshgrid(range(2*nu2+1),range(2*nv2+1),range(2*nw2+1))]).T
    # hkls = pd.DataFrame()
    print(colors.blue+'...filling...'+colors.black)
    for idx,hkl in zip(idxs,hkls):
        Fhkl[tuple(idx)] = dc.mott_bethe_factor(hkl) * grid.get_value(*hkl)
    return hklF,Fhkl
    
def import_fcf(file_path = 'shelx_thick_10A.fcf'):

    # Initialize an empty list to store your data
    data = []
    # Open the file and read it line by line
    with open(file_path, 'r') as file:
        # A flag to mark when the relevant data block starts
        data_block_start = False
        for line in file:
            # Check for the start of the data block
            if line.strip() == '_refln_index_h':
                data_block_start = True
                continue
            
            # If in the data block and the line starts with an underscore, it's a column header, so skip it
            if data_block_start and line.strip().startswith('_'):
                continue
    
            # If in the data block and the line doesn't start with an underscore, it's a data row
            if data_block_start and not line.strip().startswith('_'):
                # Check for the end of the data block (empty line or new section starting)
                if line.strip() == '' or line.strip().startswith('_'):
                    data_block_start = False
                    continue
                # Split the line by whitespace and add the resulting list to the data
                data.append(line.split())
    
    # Define the column names
    columns = ['h', 'k', 'l', 'Fc-squared', 'Fo-squared', 'sigma(Fo-squared)', 'status flag']
    
    # # Create a DataFrame
    df = pd.DataFrame(data, columns=columns)
    
    # Convert numerical columns to appropriate types
    df[['h', 'k', 'l']] = df[['h', 'k', 'l']].astype(int)
    df[['Fc-squared', 'Fo-squared', 'sigma(Fo-squared)']] = df[['Fc-squared', 'Fo-squared', 'sigma(Fo-squared)']].astype(float)
    
    df.index=[str((h,k,l)) for h,k,l in df[['h','k','l']].values]
    df['hkl']=df.index    
    return df

def plot_rocking():
    h = hkls[0]
    rock.do('set_thickness',verbose=False,thick=100,v=0)
    Sw,frames = rock.beams.loc[h,['Sw','Frame']]
    I = [rock.load(f).df_G.loc[h,'I'] for f in frames ]
    
    # x='Sw'
    x='Frame'
    plts = [rock.beams.loc[h,x],I,'b-o']
    xlab = {'Frame':'frame','Sw':r'$S_w(\AA^{-1})$','theta':r'$\theta(deg)$'}[x]
    dsp.stddisp(plts,labs=[xlab,'$I$'],figsize=(12,7))    

def plot_integrated(self,refl,new:bool=False,cm='Spectral',kin=False,zs=[],
    **kwargs):
    """plot the integrated rocking curves for selected beams as function of z

    Parameters
    ----------
    refl
        beams to consider
    kin
        True to overlay kinematic integration
    zs 
    new
        update integration
    """
    self._integrate_rocking(refl=refl,new=new)
    nbs = len(refl)
    z = self.load(0).z

    Iz = np.array([self.Iz_dyn[h] for i,h in enumerate(refl)])

    cs,legElt = dsp.getCs(cm,nbs),{}
    plts = [[z,Iz[i],cs[i],'%s' %h] for i,h in enumerate(refl)]
    if kin:
        Iz = np.array([self.Iz_kin[h] for i,h in enumerate(refl)])
        plts += [[z,Iz[i],[cs[i],'--'],''] for i,h in enumerate(refl)]
        legElt={'dyn':'k-','kin':'k--'}
    if any(zs):
        idx = [np.argmin(abs(z-z0)) for z0 in zs]
        plts+=[[zs,Iz[i][idx],[cs[i],'o'],''] for i,h in enumerate(refl)]
        
    fig,ax=dsp.stddisp(plts,labs=[r'$z(\AA)$','$I_{int}$'],legElt=legElt,**kwargs)    
    return fig,ax

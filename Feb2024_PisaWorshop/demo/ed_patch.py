from utils import displayStandards as dsp
import numpy as np

def plot_integrated(self,refl,cm='Spectral',kin=False,zs=[],
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
        
    fig,ax=dsp.stddisp(plts,labs=[r'$z(\AA)$','$I_{int}$'],legElt=legElt,**kwargs); 
    return fig,ax


def plot_Idyn_vs_Ikin(hklf,rock,Swm=0.0095):
    iZs,nzs=rock._get_iZs(iZs=None,zs=[t])
    iZ=iZs[0]#;print(t,iZ)
    
    df['Idyns'] = np.array([rock.Iz_dyn[h][iZ] for h in hklf])
    df['Ikins'] = np.array([rock.Iz_kin[h][iZ] for h in hklf])
    # df['Ikins'] = np.abs([integrate.trapz(Ikin(sw,t,xi_g),sw) for xi_g,sw in zip(df.xi_g,rock.beams.loc[hklf,'Sw']) ])
    df['rel']=np.abs(df.Idyns-df.Ikins)/df.Ikins
    
    formats={'xi_g':'{:>7.0f}','Idyns':'{:>6.2e}','Ikins':'{:>6.2e}'}
    print(df.sort_values('Ikins',ascending=False)[['xi_g','Idyns','Ikins','rel']][:10].to_string(formatters={k: v.format for k, v in formats.items()}))
    Rfactor=np.abs(df.Idyns-df.Ikins).sum()/df.Ikins.sum()*100
    print('Rfactor=%.1f' %Rfactor)
    # print(len(df.loc[df.rel<0.1].sort_values('rel')),len(hklf))

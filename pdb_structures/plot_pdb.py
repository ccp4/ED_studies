from utils import *;imp.reload(dsp)
import make_pdb_pkl as pdbf;imp.reload(pdbf)

df = pd.read_pickle(pdbf.pkl_file)
exp = df.exp.unique().astype(str)
# exp_diff=[s for s in exp if 'DIFFRACTION' in s]
exp_diff=[s for s in exp if 'ELECTRON CRYSTALLOGRAPHY' in s]
df = df.loc[df.exp.isin(exp_diff)]
df = df.loc[df.V>10]
df = df.loc[df.res<3.0]
# df = df.loc[df.sp=='P 21 21 21']
# df = df.loc[df.sp=='P 1']
df['std'] = list(map(lambda x:x.std()/x.mean(),
    [np.array(eval(s)[:3]) for s in df.params],
    ))

cols=pdbf.col_names
c_col='res'
# c_col='std'

# fig,ax=plt.subplots(layout="constrained")
fig,ax=dsp.create_fig()
cs = ax.scatter(np.log10(df.V.astype(int)),np.log10(df.nb_atoms.astype(int)),
    s=30,c=df[c_col],cmap='RdYlGn_r',vmin=0.5,vmax=4)
fig,ax = dsp.stddisp(fig=fig,ax=ax,#scat=[np.log10(df.V.astype(int)),np.log10(df.nb_atoms.astype(int)),df.res],
    labs=['$log_{10} %s$' %cols['V'],'$log_{10}$(%s)' %cols['nb_atoms']])#,imOpt='c')
# fig.colorbar(cs)
(xm,ym),(xM,yM)=cs.get_datalim(ax.transData).get_points()
to_axes=lambda x,y:[(x-xm)/(xM-xm),(y-ym)/(yM-ym)]
annot = ax.annotate("",xy=(0,0),xycoords='data',xytext=(100,100),textcoords="offset points",
                    ha='left',va='bottom',
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"),
                    fontsize=13,
                    in_layout=True)
annot.get_bbox_patch().set_facecolor('b')#cmap(norm(c[ind["ind"][0]])))
annot.get_bbox_patch().set_alpha(0.1)

annot.set_visible(False)

xylims=[]

def update_annot(ind):
    pos = cs.get_offsets()[ind["ind"][0]]#;print(pos,to_axes(pos[0],pos[1]))
    annot.xy = pos
    dat=df.iloc[ind["ind"][0]]
    a,b,c,alpha,beta,gamma = eval(dat.params)
    text = r"""{name} : {sp}
    res={res:.1f} A
    nb_atoms={nb_atoms},
    a={a:4.1f}A, α={alpha:4.1f},
    b={b:4.1f}A, β={beta: 4.1f},
    c={c:4.1f}A, γ={gamma:4.1f}""".format(
        name=dat.name,sp=dat.sp,res=dat.res,
        a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,
        nb_atoms=dat.nb_atoms)
    annot.set_text(text)
    # annot.set_usetex(True) #too slow
    (_,_),(w,h)=annot.get_bbox_patch().get_bbox().get_points()
    x0,y0=to_axes(pos[0],pos[1])
    X,Y=20,20
    if x0>0.5:
        X=-w
    if y0>0.5:
        Y=-h
    annot.xyann=(X,Y)

def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = cs.contains(event)
        if cont:
            update_annot(ind)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)

plt.show()

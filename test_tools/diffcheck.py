#!/usr/bin/python3
import visualize as vs
pd_df, nxy_df = vs.geom_info('../hyper_res/')
newr1s = vs.parse_results_all('../hyper_res')
refr1s = vs.parse_results_all('../ref_res')
d_df , rd_df = vs.diff(newr1s[10],refr1s[10])

part_dfs = [ pd_df.loc[ pd_df['P'] == i ] for i in range(int(pd_df.values.min()),int(pd_df.values.max()+1))]
parts = [ part_df.join(nxy_df)[['x','y']].values for part_df in part_dfs ]


def show(col, pm = False) :
    vs.plt.figure(1)
    vs.plot_stuff(nxy_df,d_df,col)
    vs.plt.title("Delta, " + col)
    if pm:
        for part in parts:
            vs.plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
    vs.plt.figure(2)
    
    vs.plot_stuff(nxy_df,rd_df,col)
    vs.plt.title("Relative Delta, " + col)
    if pm: 
        for part in parts:
            vs.plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
    vs.plt.show()


show('ND')
show('U')
show('V')
show('RHO')
show('PS')
show('TEMP')

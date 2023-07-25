import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
import numpy as np
import scipy.constants as pc
import astropy.constants as ac
import os

def h_norm(m_a, m_d, f, r):
    '''
    Compute h_norm for a given binary system. This computation assumes that the binary system is detached.
    m_a: mass of accretor, in M_sun
    m_d: mass of donor, in M_sun
    f: orbital frequency
    r: distance to system
    '''
    M_tot = m_a + m_d
    q = m_d / m_a
    Q = q / ((1 + q) ** 2)
    M_ch = M_tot * Q ** (3/5)

    J_orb = (pc.G ** 2 * M_tot ** 5 / (np.pi * f)) ** (1 / 3) * Q
    h_norm = 4*pc.G**3/(r*pc.c**4)*M_ch**5/J_orb**2
    return h_norm

def rh_norm(m_a_in, m_d_in, f, to_kg=True):
    '''
    Compute rh_norm for a given binary system, the normalized gravitational strain of the system.
    This computation assumes that the binary system can be treated as two point masses.
    m_a: mass of accretor, in M_sun
    m_d: mass of donor, in M_sun
    f: orbital frequency
    '''
    if to_kg:
        m_a = m_a_in*ac.M_sun.value
        m_d = m_d_in*ac.M_sun.value
    else:
        m_a = m_a_in
        m_d = m_d_in

    M_tot = m_a + m_d
    q = m_d / m_a
    Q = q / ((1 + q) ** 2)
    M_ch = (m_a*m_d)**(3/5)/(m_a+m_d)**(1/5)

    J_orb = (pc.G ** 2 * M_tot ** 5 / (np.pi * f)) ** (1 / 3) * Q
    h_norm = 4*pc.G**3/(pc.c**4)*M_ch**5/J_orb**2
    return h_norm

def f_gw(P_orb, to_sec=True):
    '''
    Gravitational wave frequency from the orbital period. If to_sec, input is in days
    and output is hz. Otherwise, a conversion needs to be made after the function is run.
    '''
    if to_sec:
        f_orb = 1/(P_orb*60*60*24)
    else:
        f_orb = 1/P_orb

    return 2*f_orb

def r_RL(m_a,m_d,a):
    '''
    Computes the Roche Lobe radius from the component masses and sep
    '''
    q = m_d / m_a
    R_l = a*(0.49*q**2/3)/(0.6*q**(2/3)+np.log(1+q**(1/3)))
    return R_l

def is_contact_binary(m_a,m_d,a):
    '''
    Determine if the binary system is a contact binary undergoing mass transfer.
    a: orbital radius of binary system (in solar radii)
    m_a: mass of accretor, in M_sun
    m_d: mass of donor, in M_sun
    '''
    q = m_d / m_a

    R_d = 0.0114*((m_d/1.44)**(-2/3)-(m_d/1.44)**(2/3))**(1/2)*(1+3.5*(m_d/0.00057)**(-2/3)+(m_d/0.00057)**(-1))**(-2/3) #solar radii
    R_l = a*(0.49*q**2/3)/(0.6*q**(2/3)+np.log(1+q**(1/3)))

    return R_d >= R_l

def dwd_bin_type(kstar_1,kstar_2):
    if kstar_1 < kstar_2:
        klg = kstar_2
        ksm = kstar_1
    else:
        klg = kstar_1
        ksm = kstar_2

    if ksm < 10 or klg > 14:
        return 'Not a binary'

    t = {10:'He',11:'Co',12:'O/Ne',13:'NS',14:'BH'}
    #100 for He-He, 110 for He-Co, 120 for He-O/N
    return f'{t[ksm]}-{t[klg]}'

def df_bin_type(df):
    return (df[df['WDbintype'] == 'He-He'], df[df['WDbintype'] == 'He-Co'], df[df['WDbintype'] == 'Co-Co'],
                    df[df['WDbintype'] == 'He-O/Ne'], df[df['WDbintype'] == 'Co-O/Ne'], df[df['WDbintype'] == 'O/Ne-O/Ne'])


def process_taumaline(f_name,nrows=None):
    #columns to use
    usecols = ['#bin_num','mass_1','mass_2','kstar_1','kstar_2','t_birth','sep_final','porb_final']
    df = pd.read_csv(f_name,sep=',',usecols=usecols,nrows=nrows)

    df.insert(2,'DWD_bin_type', \
                [dwd_bin_type(df.at[i,'kstar_1'],df.at[i,'kstar_2']) for i in range(df.shape[0])])
    df.drop(labels=['kstar_1','kstar_2'],axis=1,inplace=True)

    df.insert(3,'Is_contact_binary',is_contact_binary(df['mass_1'], df['mass_2'],1/2*df['sep_final']))
    #df.drop(labels=['sep_final'],axis=1,inplace=True)

    #calculate gravitational wave frequency from planetary orbital period
    df.insert(4, 'f_gw', f_gw(df['porb_final']))
    #df.drop(labels=['porb_final'],axis=1,inplace=True) #we don't need planetary orbital period anymore
    #calculate rh_norm
    df.insert(5,'rh_norm',rh_norm(df['mass_1'], df['mass_2'], df['f_gw']))
    #df.drop(labels=['mass_1','mass_2'], axis=1, inplace=True) #we do not need mass anymore


    return df

def dwd_type_cb_subplots(df,outfile=None):
    '''
    Take df and make subplots by binary system type and if the system is in mass transfer
    '''
    dwd_type = ['He-He','He-Co','Co-Co','He-O/Ne','Co-O/Ne']
    is_cb = [False, True]
    cb_label = ['Detached', 'Contact']

    l_f_max = np.max(np.log10(np.array(np.array(df['f_gw']))))
    l_rh_max = np.max(np.log10(np.array(np.array(df['rh_norm']))))
    l_f_min = np.min(np.log10(np.array(np.array(df['f_gw']))))
    l_rh_min = np.min(np.log10(np.array(np.array(df['rh_norm']))))

    t_birth_max = df['t_birth'].max()
    t_birth_min = df['t_birth'].min()

    fig, axes = plt.subplots(len(dwd_type),len(is_cb),sharex=True,sharey=True,figsize=(8,6))
    norm = mpl.colors.Normalize(vmin=t_birth_min,vmax=t_birth_max)
    cmap = mpl.colormaps['cividis']

    for i in range(len(dwd_type)):
        for j in range(len(is_cb)):
            df_ij = df[df['DWD_bin_type'] == dwd_type[i]][df['Is_contact_binary'] == is_cb[j]]
            axes[i][j].set(xlim=(l_f_min,l_f_max),ylim=(l_rh_min,l_rh_max))
            axes[i][j].scatter(np.log10(np.array(df_ij['f_gw'])),
                                np.log10(np.array(df_ij['rh_norm'])),
                                c=np.array(df_ij['t_birth']),
                                norm=norm,
                                cmap=cmap,
                                label=f'{dwd_type[i]}',
                                marker='.',lw=0,s=1,alpha=0.5)

            if j == 0:
                axes[i][j].set_ylabel(r'$\log_{10}{(rh_{norm})}$')
            if j == (len(is_cb) - 1):
                axes[i][j].yaxis.set_label_position('right')
                axes[i][j].set_ylabel(f'{dwd_type[i]}')
            #axes[i][j].legend()

    axes[0][0].set_title('Detached Binary')
    axes[0][1].set_title('Contact Binary')

    axes[-1][0].set_xlabel(r'$\log_{10}{(f_{gw} (hz))}$')
    axes[-1][1].set_xlabel(r'$\log_{10}{(f_{gw} (hz))}$')
    fig.suptitle('Taumaline Galaxy Kt Plots')


    fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes)
    if outfile:
        plt.savefig(outfile)
    plt.show()

def m_1_m_2_from(K,q):
    '''
    Calculate mass of primary and secondary from mass parameter K and mass ratio q
    as defined in KT. Note that 0<K<1 (for White Dwarfs) and 0<q<1
    '''
    m_1 = 1.44*K/(2**(1/5))*((1+q)/q**3)**(1/5) #M_ch*K/(2^1/5)*((1+q)/q^3)^1/5
                                                #where M_ch is the Chandrasekar mass
                                                #Note that m_1 is taken to be primary by COSMIC
                                                #and is M_d in KT
    m_2 = q*m_1 #Mass of secondary related by mass ratio to mass of primary
    return (m_1,m_2)

def bound_A(log10_f):
    return 0.731 + 2/3 * log10_f

def bound_B(log10_f):
    return 0.703 + 0.637*log10_f - 0.017*log10_f**2 + 0.298*log10_f**3 + 0.061*log10_f**4

def bound_C(log10_f):
    return 0.761 + 1.005*log10_f + 0.700*log10_f**2 +0.700*log10_f**3 + 0.214*log10_f**4 + 0.023*log10_f**5

def bound_D(log10_f):
    return 2.141+1.686*log10_f - 0.141*log10_f**2 + 0.007*log10_f**3

def bound_E(log10_f):
    return -1.381 - 2.108*log10_f - 1.394*log10_f**2 - 0.167*log10_f**3

def mt_lim_from_M_q(M,q):
    '''
    Determine mass transfer limit for system with total mass M=m_a+m_d and mass ratio q=m_d/m_a
    For q=2/3 and all M, essentially equal to bound_D
    '''
    m_d = M/(1+1/q)
    m_a = m_d/q
    Q = q/((1+q)**2)
    M_ch = 1.44 #Chandrasekar mass in M_sun
    M_p = 0.00057
    M_sun = ac.M_sun.value
    R_sun = ac.R_sun.value


    a = 0.0114 * ((m_d/M_ch) ** (-2/3) - (m_d/M_ch) ** (2/3)) ** (1/2) * (1+3.5*(m_d/M_p)**(-2/3)+(m_d/M_p)**(-1))**(-2/3) * \
                (0.6*q**(2/3)+np.log(1+q**(1/3))) / (0.49*q**(2/3))

    Omega = np.sqrt(pc.G*M*M_sun/((a*R_sun)**3)) #last term to convert from /year to /s
    f_gw = Omega/np.pi
    J_orb = (m_a*M_sun*m_d*M_sun*(a*R_sun)**2)/(M*M_sun) * Omega
    rh_norm = 4*pc.G**3/pc.c**4*((M*M_sun)**5*Q**3)/J_orb**2

    return (np.log10(f_gw),np.log10(rh_norm))

def mt_lim_from_m_a_q(m_a,q):
    '''
    Determine mass transfer limit for system with m_a=M_ch and mass ratio q=m_d/m_a
    '''
    m_d = q*m_a
    M = m_a + m_d
    Q = q/((1+q)**2)
    M_ch = 1.44 #Chandrasekar mass in M_sun
    M_p = 0.00057 #
    M_sun = ac.M_sun.value
    R_sun = ac.R_sun.value


    a = 0.0114 * ((m_d/M_ch) ** (-2/3) - (m_d/M_ch) ** (2/3)) ** (1/2) * (1+3.5*(m_d/M_p)**(-2/3)+(m_d/M_p)**(-1))**(-2/3) * \
                (0.6*q**(2/3)+np.log(1+q**(1/3))) / (0.49*q**(2/3))

    Omega = np.sqrt(pc.G*M*M_sun/((a*R_sun)**3)) #last term to convert from /year to /s
    f_gw = Omega/np.pi
    J_orb = (m_a*M_sun*m_d*M_sun*(a*R_sun)**2)/(M*M_sun) * Omega
    rh_norm = 4*pc.G**3/pc.c**4*((M*M_sun)**5*Q**3)/J_orb**2

    return (np.log10(f_gw),np.log10(rh_norm))

def plot_bin_trajectory(bpp, bcm, fig, ax, system_label=None, marker='x', \
                            end_inspiral_label=None,cmap=None,norm=None,color=None):
    M_sun = ac.M_sun.value
    mass_1 = bcm['mass_1']
    mass_2 = bcm['mass_2']

    bcm['f_gw'] = f_gw(bcm['porb'])
    bcm['rh_norm'] = rh_norm(mass_1, mass_2, bcm['f_gw'])

    #plot trajectory
    x = np.log10(np.array(bcm['f_gw']))
    y = np.log10(np.array(bcm['rh_norm']))
    if cmap:
        ax.scatter(x,y,c=bcm['tphys'],marker=marker,label=system_label,cmap=cmap,norm=norm)
    else:
        ax.scatter(x,y,marker=marker,label=system_label,color=color,lw=0.5,s=0.5)


    #plot transition to mass transfer
    porb = bpp[bpp['evol_type']==3]['porb']
    f_mt = f_gw(porb)
    rh_norm_mt = rh_norm(mass_1.iloc[0],mass_2.iloc[0],f_mt)
    ax.scatter(np.log10(f_mt), np.log10(rh_norm_mt), label=end_inspiral_label, marker='*',color=color)

    #ax.legend()

    return (fig, ax)

def KT_from_hdf5_dir(d_name, out_dir='Graphs'):
    '''
    Create a KT diagram for each h5 file in a directory of h5 files
    '''
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    fig, ax = plt.subplots()
    colors = ['black','red','cadetblue','pink','deeppink','palevioletred','magenta','purple','darkviolet','blueviolet',
              'mediumslateblue','blue','navy','coral','aqua','bisque','khaki']

    q = np.linspace(0,1,100)
    M = np.linspace(0,2.88,100)
    log10_f_1 = np.linspace(-4,0.1,100)
    log10_f_2 = np.linspace(-4,-1.38,100)
    log10_f_3 = np.linspace(-4,-1.55,100)

    for f in os.listdir(d_name):
        fig_s,ax_s = plt.subplots()
        legend_elements = []
        title = f.strip('gx_dat_.h5')
        out_name = os.path.join(out_dir, title + '.jpg')

        f_name = os.path.join(d_name, f)
        print(f'[*] Reading in {f_name}')

        try:
            df = pd.read_hdf(f_name,key='LISA_population')
        except:
            continue

        df['f_gw'] = f_gw(df['porb_final'])
        df['rh_norm'] = rh_norm(df['mass_1'], df['mass_2'], df['f_gw'])
        evol_types = df['evol_type'].unique()
        print(f'[*] Found evol_types: {evol_types}')



        for t in evol_types:
            df_t = df[df['evol_type'] == t]
            ax.scatter(np.log10(df_t['f_gw']), np.log10(df_t['rh_norm']), color=colors[int(t)],marker='.',
                       lw=0.1,s=1)
            ax_s.scatter(np.log10(df_t['f_gw']), np.log10(df_t['rh_norm']), color=colors[int(t)],marker='.',
                       lw=0.5,s=1)
            legend_elements.append(Line2D([0], [0], color=colors[int(t)], marker='o', \
                                            label=f'Evol type: {t}', markersize=10))

        ax_s.plot(log10_f_1,bound_A(log10_f_1),color='tab:red')
        ax_s.plot(*mt_lim_from_M_q(M,1),color='tab:red')
        ax_s.plot(*mt_lim_from_m_a_q(1.44,q),color='tab:purple',lw=0.9)
        ax_s.plot(log10_f_2,bound_D(log10_f_2),color='tab:cyan',lw=0.9)
        ax_s.plot(log10_f_3,bound_E(log10_f_3),color='tab:green',lw=0.9)

        ax_s.legend(handles=legend_elements)
        ax_s.set_xlabel(r'$\log_{10}{(f_{gw})}$')
        ax_s.set_ylabel(r'$\log_{10}{(rh_{norm})}$')
        ax_s.set_title(f'{title}')
        ax_s.grid()
        ax_s.set_xlim(-5,-1)
        ax_s.set_ylim(-4,-1)
        fig_s.savefig(out_name)

        legend_elements = []


    ax.plot(log10_f_1,bound_A(log10_f_1),color='tab:red')
    ax.plot(*mt_lim_from_M_q(M,1),color='tab:red')
    ax.plot(*mt_lim_from_m_a_q(1.44,q),color='tab:purple',lw=0.9)
    ax.plot(log10_f_2,bound_D(log10_f_2),color='tab:cyan',lw=0.9)
    ax.plot(log10_f_3,bound_E(log10_f_3),color='tab:green',lw=0.9)
    #plt.legend()
    ax.set_xlabel(r'$\log_{10}{(f_{gw})}$')
    ax.set_ylabel(r'$\log_{10}{(rh_{norm})}$')
    ax.set_title('All Binary Systems')
    ax.grid()
    ax.set_xlim(-5,-1)
    ax.set_ylim(-4,-1)
    fig_s.savefig(os.path.join(out_dir,'All_Systems.pdf'))
    plt.show()

def evol_type_tphys(gal_pop, fix_pop, out_dir='Evol_Types'):
    title = gal_pop.strip('/Tau/hdf5/gx_dat_.h5')
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    out_2 = os.path.join(out_dir,'final_type_2')
    out_4 = os.path.join(out_dir,'final_type_4')
    try:
        os.mkdir(out_2)
    except FileExistsError:
        pass

    try:
        os.mkdir(out_4)
    except FileExistsError:
        pass


    #first, we want to get unique bin nums in the galaxy population
    print(f'reading in: {gal_pop}')
    df = pd.read_hdf(gal_pop,key='LISA_population')
    bn_2 = np.array(df[df['evol_type']==2]['bin_num'].unique()) #unique bin nums in df for evol type 2
    bn_4 = np.array(df[df['evol_type']==4]['bin_num'].unique()) #and the same for evol type 4
    del df #clear up this memory so it can be used when we load up the fix_pop

    #Get the bpp dataframe from the fix_pop, which includes evolutionairy history
    print(f'reading in: {fix_pop}')
    bpp = pd.read_hdf(fix_pop,key='bpp')
    bpp_2 = bpp[bpp['bin_num'].isin(bn_2)] #filter for stars with evol_type 2 in galaxy population
    bpp_4 = bpp[bpp['bin_num'].isin(bn_4)] #and same for evol_type 4 in galaxy population


    colors = ['white','red','green','hotpink','crimson','violet','darkorange','yellow','lightsteelblue','black', \
                        'silver', 'magenta','limegreen','cyan','blue','navy','deeppink'] #colors for each possible evol type
                                                                                            #(white place holder so access by evol type)
    #create plot for bpp_2
    evol_types = bpp_2['evol_type'].unique()
    fig_2, ax_2 = plt.subplots()

    for t in evol_types:
        ax_2.scatter(bpp_2[bpp_2['evol_type']==t]['tphys'], \
                     np.log10(bpp_2[bpp_2['evol_type']==t]['kstar_1']*bpp_2[bpp_2['evol_type']==t]['kstar_2']), \
                        color=colors[int(t)],marker='o',alpha=0.2,label=f'evol_type:{t}')

    ax_2.set_xlabel(r'Evolution Time (Myr)')
    ax_2.set_ylabel(r'$\log{kstar_1*kstar_2}$')
    ax_2.set_title(title)
    ax_2.legend()
    ax_2.grid()
    fig_2.savefig(os.path.join(out_2,title + '.jpg'))

    #create plot for bpp_4
    evol_types = bpp_4['evol_type'].unique()
    fig_4, ax_4 = plt.subplots()

    for t in evol_types:
        ax_4.scatter(bpp_4[bpp_4['evol_type']==t]['tphys'], \
                     np.log10(bpp_4[bpp_4['evol_type']==t]['kstar_1']*bpp_4[bpp_4['evol_type']==t]['kstar_2']), \
                        color=colors[int(t)],marker='o',alpha=0.2,label=f'evol_type:{t}')

    ax_4.set_xlabel(r'Evolution Time (Myr)')
    ax_4.set_ylabel(r'$\log{kstar_1*kstar_2}$')

    ax_4.set_title(title)
    ax_4.grid()
    ax_4.legend()
    fig_4.savefig(os.path.join(out_4,title + '.jpg'))


def get_final_fix(bin_nums,bpp):
    '''
    Create a dataframe from a bpp dataframe which includes the final state for each star, as well as the evol_type it was in
    when the simulation ended.

    bin_nums: list of bin_nums to include in the output array
    bpp: bpp dataframe
    '''
    fix_list = [pd.DataFrame(columns=bpp.columns)]

    for num in bin_nums:
        sys = bpp[bpp.bin_num == num]
        row = sys.iloc[[-1]]
        row.evol_type = float(sys.iloc[-2].evol_type)

        fix_list.append(row)


    return pd.concat(fix_list)

def evol_graph(kstar1, kstar2, final_evol_type, fix_pop, nstars, seed, fig=None, ax=None, out_dir=None):
    '''
    Create a KT diagram tracking the evolution of a sample of stars as they move throught their stellar evolution

    '''
    evol_colors = {1.0 : 'lightgrey', 2.0 : 'purple', 3.0 : 'cyan', 4.0 : 'deeppink', 5.0 : 'lightgreen',
               6.0 : 'darkorange', 7.0 : 'lime', 8.0 : 'blue', 9.0 : 'black', 10.0 : '#0f0f0f0f', 11.0 : 'red',
               12.0 : 'mediumaquamarine', 13.0 : 'aquamarine', 14.0 : 'navy', 15.0 : 'magenta', 16.0 : 'crimson'}

    fp = {'bulge' : 'fixed_fix_pop_bulge.csv',
          'thinDisk' : 'fixed_fix_pop_thinDisk.csv',
          'thickDisk' : 'fixed_fix_pop_thickDisk.csv'}
    b_name = ''
    if kstar1 == 12:
        b_name = 'dat_kstar1_12_kstar2_10_12.h5'
    else:
        b_name = f'dat_kstar1_{int(kstar1)}_kstar2_{int(kstar2)}.h5'
    #First, look in fixed pop to get some stars
    fix = pd.read_csv(fp[fix_pop])
    fix = fix.query(f'kstar_1 == {kstar1} & kstar_2 == {kstar2} & evol_type == {final_evol_type}')
    bin_nums = np.array(fix.sample(nstars,random_state=seed)['bin_num'])
    del fix

    #Now, lets get these stars stellar histories
    base = 'tau_hdf5'
    fixed = os.path.join('tau_hdf5','fixed_pop')
    subs = {'bulge' : 'fixed_bulge','thickDisk' : 'fixed_thickDisk','thinDisk' : 'fixed_thinDisk'}
    b_name = os.path.join(os.path.join(fixed,subs[fix_pop]),b_name)

    bpp = pd.read_hdf(b_name, key='bpp') #Full stellar history
    bpp = bpp[bpp['bin_num'].isin(bin_nums)]
    bpp['f_gw'] = f_gw(bpp['porb'])
    bpp['rh_norm'] = rh_norm(bpp['mass_1'],bpp['mass_2'],bpp['f_gw'])

    if ax == None:
        fig,ax = plt.subplots()
    ax.set_xlabel(r'$\log_{10}{(f_{gw})}$')
    ax.set_ylabel(r'$\log_{10}{(rh_{norm})}$')
    ax.set_title(f'Evolutionairy History For kstar1={kstar1}, kstar2={kstar2}')
    ax.axvline(x=-5,color='black',linestyle='dashed')
    ax.grid()

    cmap = mpl.colormaps['viridis']
    max_diff = np.diff(bpp.tphys).max()
    norm = mpl.colors.PowerNorm(1/3,0,max_diff)
    #fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap))

    for i in bin_nums:
        df = bpp[bpp['bin_num'] == i]
        dt = np.diff(df.tphys)

        log10_f_gw = np.log10(np.array(df['f_gw']))
        log10_rh_norm = np.log10(np.array(df['rh_norm']))
        for j in range(1,len(dt)+1):
            color = cmap(norm(dt[j-1]+0.1))
            ax.plot([log10_f_gw[j-1],log10_f_gw[j]],[log10_rh_norm[j-1],log10_rh_norm[j]],marker='.', \
                            color=color,lw=0.2,markersize=0)

        for key in evol_colors.keys():
            df_evol = df[df['evol_type']==key]
            ax.scatter(np.log10(df_evol['f_gw']),np.log10(df_evol['rh_norm']),marker='*',color=evol_colors[key])

    if out_dir:
        try:
            os.mkdir(out_dir)
        except FileExistsError:
            pass

        fig.savefig(os.path.join(out_dir, f'kstar1_{kstar1}_kstar2_{kstar2}_evol_type_{final_evol_type}_seed_{seed}.pdf'))

    return bpp, fig, ax

def plot_by_kstar(fig,ax,df):
    '''
    Create a KT diagram from a df with f_gw and rh_norm columns
    color code by kstar-kstar
    '''
    colors = ['tab:green','tab:red','tab:cyan','tab:purple','tab:pink','tab:orange']
    ax.scatter(
        np.log10(df[(df['kstar_1']==10) & (df['kstar_2']==10)]['f_gw']),
        np.log10(df[(df['kstar_1']==10) & (df['kstar_2']==10)]['rh_norm']),
        marker='.',color=colors[0],lw=1,s=1)
    ax.scatter(
        np.log10(df[(df['kstar_1']==11) & (df['kstar_2']==10)]['f_gw']),
        np.log10(df[(df['kstar_1']==11) & (df['kstar_2']==10)]['rh_norm']),
        marker='.',color=colors[1],lw=1,s=1)
    ax.scatter(
        np.log10(df[(df['kstar_1']==11) & (df['kstar_2']==11)]['f_gw']),
        np.log10(df[(df['kstar_1']==11) & (df['kstar_2']==11)]['rh_norm']),
        marker='.',color=colors[2],lw=1,s=1)
    ax.scatter(
        np.log10(df[(df['kstar_1']==12) & (df['kstar_2']==10)]['f_gw']),
        np.log10(df[(df['kstar_1']==12) & (df['kstar_2']==10)]['rh_norm']),
        marker='.',color=colors[3],lw=1,s=1)
    ax.scatter(
        np.log10(df[(df['kstar_1']==12) & (df['kstar_2']==11)]['f_gw']),
        np.log10(df[(df['kstar_1']==12) & (df['kstar_2']==11)]['rh_norm']),
        marker='.',color=colors[4],lw=1,s=1)
    ax.scatter(
        np.log10(df[(df['kstar_1']==12) & (df['kstar_2']==12)]['f_gw']),
        np.log10(df[(df['kstar_1']==12) & (df['kstar_2']==12)]['rh_norm']),
        marker='.',color=colors[5],lw=1,s=1)


    legend_elements = [Line2D([0],[0],color=colors[0],label=f'He-He',marker='o',lw=0),
                       Line2D([0],[0],color=colors[1],label=f'C-He',marker='o',lw=0),
                       Line2D([0],[0],color=colors[2],label=f'C-C',marker='o',lw=0),
                       Line2D([0],[0],color=colors[3],label=f'O/Ne-He',marker='o',lw=0),
                       Line2D([0],[0],color=colors[4],label=f'O/Ne-C',marker='o',lw=0),
                       Line2D([0],[0],color=colors[5],label=f'O/Ne-O/Ne',marker='o',lw=0)]

    ax.legend(handles=legend_elements,loc='lower right')

    return fig, ax


def plot_boundaries(fig,ax):
    '''
    Plot KT boundaries
    '''
    M = np.linspace(0,2.88,100)
    q = np.linspace(0,1,100)
    log10_f = np.linspace(-5,0,100)

    ax.plot(log10_f, bound_A(log10_f),color='black')
    ax.plot(*mt_lim_from_M_q(M,1),color='black')
    ax.plot(*mt_lim_from_m_a_q(1.44,q),color='black',ls='dashed')
    ax.fill_between(*mt_lim_from_M_q(M,1),-5,color='lightgrey')
    ax.fill_between(log10_f, 1, bound_A(log10_f), color='lightgrey')

    return fig, ax

def csv_from_fix_pops(base_dir,out_name='agate_csv.csv'):
    '''
    Take as input a directory which contains sub directories 'Bulge', 'thickDisk', and 'thinDisk'
    which each contain one or more fixed pop .h5 files. Take the bcm df from each of these files,
    concat them all together, and save them to a csv!
    '''
    sub_dirs = ['Bulge','thickDisk','thinDisk']
    df_list = []

    for dir in sub_dirs:
        d = os.path.join(base_dir,dir)
        for file in os.listdir(d):
            if file[-3:] == '.h5':
                #Add the 'bcm' file, which will contain the final state of all the stars
                #(We lose evol_type this way, but the alternative–using get_final_fix–would take very long)
                df_list.append(pd.read_hdf(os.path.join(d,file),'bcm'))

    out_df = pd.concat(df_list)
    out_df['f_gw'] = f_gw(out_df['porb'])
    out_df['rh_norm'] = rh_norm(out_df['mass_1'],out_df['mass_2'],out_df['f_gw'])
    out_df.to_csv(out_name)
    return

def snr(df, nc_file):
    '''
    Calculate the LISA signal to noise ratio for all systems in the given df
    from a LISA noise curve file
    df: dataframe to calculate the signal to noise ratio for.
    nc_file: path to noise curve file

    Returns: numpy array with signal to noise ratios for each system in df.
    '''
    dat = np.loadtxt(os.path.join('Agate','LISA2018_esaSTD.csv'),delimiter=',')
    f = dat[:,0]
    noise = dat[:,1]

    m_chirp = (np.array(df.mass_1)*np.array(df.mass_2))**(3/5)/(np.array(df.mass_1)+np.array(df.mass_2))**(1/5)
    kappa = (5 / (256 * np.pi ** (8/3)) / (pc.G/pc.c**3)**(5/3)) * (1 / ((m_chirp * ac.M_sun.value) ** (5/3)))
    t_obs = 1 * (60*60*24*365.25) #Observation time in years
    stat_freq = (8/3*(kappa/t_obs**2)**(3/11)) #Cutoff frequency to be considered a stationairy source

    is_monochrome = [df.f_gw < stat_freq]
    h_norm = df['rh_norm'] / (8*ac.kpc.value) #h_norm if source were at a distance of 8 parsecs from detector

    snr = np.empty(m_chirp.size)
    snr = h_norm * np.sqrt(t_obs) / np.interp(np.array(df.f_gw),f,noise) #np.interp interpolates the noisecurve at f_gw

    return snr 

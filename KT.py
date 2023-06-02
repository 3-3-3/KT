import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as pc

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

def rh_norm(m_a, m_d, f):
    '''
    Compute h_norm for a given binary system. This computation assumes that the binary system is detached.
    m_a: mass of accretor, in M_sun
    m_d: mass of donor, in M_sun
    f: orbital frequency
    '''
    M_tot = m_a + m_d
    q = m_d / m_a
    Q = q / ((1 + q) ** 2)
    M_ch = (m_a*m_d)**(3/5)/(m_a+m_d)**(1/5)

    J_orb = (pc.G ** 2 * M_tot ** 5 / (np.pi * f)) ** (1 / 3) * Q
    h_norm = 4*pc.G**3/(pc.c**4)*M_ch**5/J_orb**2
    return h_norm

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

def bin_type(kstar_1,kstar_2):
    t = {10:'He',11:'Co',12:'O/Ne',13:'NS',14:'BH'}
    return f'{t[kstar_1]}-{t[kstar_2]}'

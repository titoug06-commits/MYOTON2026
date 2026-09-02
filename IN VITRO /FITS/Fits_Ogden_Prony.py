"""
Script d'identification des paramètres hyperélastiques et viscoélastiques.
Compare les modèles Néo-Hookéen vs Ogden ainsi que les séries de Prony (n=1 vs n=2).

A exécuter dans le même répertoire que là ou se situent les données brutes (.tracking.csv)

Génère les graphiques PNG et un fichier de synthèse CSV.
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Configuration 
diametre_piston = 29.51
A0 = np.pi * (diametre_piston / 2)**2

materials = [
    {   'name': 'Ecoflex 00-10 (3.36 mm)',
        'L0': 3.36,
        'cyc_file': 'Ecoflex00-10_336_1.Boucle 1.tracking.csv',
        'rel_file': 'Ecoflex00-10_336_relax_test2.steps.tracking.csv',
        'rel_tstart': 1.0, 'rel_tpeak': 5.25, 'rel_tend': 305.85,
        'cyc_tstart': 640.0, 'cyc_tpeak': 709.0},

    {   'name': 'Ecoflex 00-10 (11.41 mm)',
        'L0': 11.41,
        'cyc_file': 'Ecoflex00-10_1141_4.Boucle 1.tracking.csv',
        'rel_file': 'Ecoflex00-10_1141_relax_test1.steps.tracking.csv',
        'rel_tstart': 0.6, 'rel_tpeak': 4.05, 'rel_tend': 304.60,
        'cyc_tstart': 475.0, 'cyc_tpeak': 519.15},

    {   'name': 'Ecoflex 00-50 (26.03 mm)',
        'L0': 26.03,
        'cyc_file': 'Ecoflex00-50_2603_1.Boucle 1.tracking.csv',
        'rel_file': 'Ecoflex00-50_2603_relax_test1.steps.tracking.csv',
        'rel_tstart': 0.15, 'rel_tpeak': 1.6, 'rel_tend': 302.10,
        'cyc_tstart': 280.0, 'cyc_tpeak': 288.35},

    {   'name': 'Dragon Skin 30 (1.56 mm)',
        'L0': 1.56,
        'cyc_file': 'Dragon30_156_1.Boucle 1.tracking.csv',
        'rel_file': 'Dragon30_156_relax_test3.steps.tracking.csv',
        'rel_tstart': 4.0, 'rel_tpeak': 12.15, 'rel_tend': 313.90,
        'cyc_tstart': 906.0, 'cyc_tpeak': 1003.35},

    {   'name': 'Dragon Skin 30 (2.70 mm)',
        'L0': 2.70,
        'cyc_file': 'Dragon30_270_2.Boucle 1.tracking.csv',
        'rel_file': 'Dragon30_270_relax_test1.steps.tracking.csv',
        'rel_tstart': 0.5, 'rel_tpeak': 7.6, 'rel_tend': 308.10,
        'cyc_tstart': 608.0, 'cyc_tpeak': 678.3}
]

# Fonctions utilitaires
def compute_r2(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1.0 - (ss_res / ss_tot)

def ogden_stress(lmbda, mu1, alpha1):
    return mu1 * (lmbda**(alpha1 - 1.0) - lmbda**(-0.5 * alpha1 - 1.0))

def force_boltzmann_n1(t_array, F_max, t_ramp, g1, tau1):
    g_inf = 1.0 - g1
    D = g_inf * t_ramp + g1 * tau1 * (1.0 - np.exp(-t_ramp / tau1))
    
    F_pred = np.zeros_like(t_array)
    mask_charge = t_array <= t_ramp
    t_c = t_array[mask_charge]
    F_pred[mask_charge] = (F_max / D) * (g_inf * t_c + g1 * tau1 * (1.0 - np.exp(-t_c / tau1)))
    
    mask_relax = t_array > t_ramp
    t_r = t_array[mask_relax] - t_ramp
    F_pred[mask_relax] = (F_max / D) * (g_inf * t_ramp + g1 * tau1 * (1.0 - np.exp(-t_ramp / tau1)) * np.exp(-t_r / tau1))
    
    return F_pred

def force_boltzmann_complete(t_array, F_max, t_ramp, g1, tau1, g2, tau2):
    g_inf = 1.0 - g1 - g2
    D = g_inf * t_ramp + g1 * tau1 * (1.0 - np.exp(-t_ramp / tau1)) + g2 * tau2 * (1.0 - np.exp(-t_ramp / tau2))
    
    F_pred = np.zeros_like(t_array)
    mask_charge = t_array <= t_ramp
    t_c = t_array[mask_charge]
    F_pred[mask_charge] = (F_max / D) * (g_inf * t_c + g1 * tau1 * (1.0 - np.exp(-t_c / tau1)) + g2 * tau2 * (1.0 - np.exp(-t_c / tau2)))
    
    mask_relax = t_array > t_ramp
    t_r = t_array[mask_relax] - t_ramp
    F_pred[mask_relax] = (F_max / D) * (g_inf * t_ramp + g1 * tau1 * (1.0 - np.exp(-t_ramp / tau1)) * np.exp(-t_r / tau1) + g2 * tau2 * (1.0 - np.exp(-t_ramp / tau2)) * np.exp(-t_r / tau2))
    
    return F_pred

# Boucle d'identification
results = []

for mat in materials:
    # 1. Chargement & Fit Hyperélastique
    cyc_df = pd.read_csv(mat['cyc_file'], sep=';', decimal=',')
    time_c = cyc_df['Temps total (s)'].values
    pos_c = cyc_df['Pos (µm)'].values / 1000.0
    force_c = cyc_df['Effort_5N (N)'].values

    mask_c = (time_c >= mat['cyc_tstart']) & (time_c <= mat['cyc_tpeak'])
    u_c = pos_c[mask_c] - pos_c[mask_c][0]
    F_c = force_c[mask_c] - force_c[mask_c][0]

    eps_c = u_c / mat['L0']
    lambda_c = 1.0 + eps_c 
    P_exp = F_c / A0

    # A. Fit Néo-Hookéen
    def obj_nh(params):
        mu1 = params[0]
        P_pred = ogden_stress(lambda_c, mu1, 2.0)
        return np.mean((P_exp - P_pred)**2)

    res_nh = minimize(obj_nh, [0.01], bounds=[(1e-5, 50.0)])
    mu1_nh = res_nh.x[0]
    P_pred_nh = ogden_stress(lambda_c, mu1_nh, 2.0)
    r2_nh = compute_r2(P_exp, P_pred_nh)

    # B. Fit Ogden
    def obj_ogden(params):
        mu1, alpha1 = params
        P_pred = ogden_stress(lambda_c, mu1, alpha1)
        return np.mean((P_exp - P_pred)**2)
        
    bounds_ogden_n2 = [(1e-5, 50.0), (0.1, 5.0)]
    res_ogden = minimize(obj_ogden, [0.01, 2.0], method='Nelder-Mead', bounds=bounds_ogden_n2,
                         options={'xatol': 1e-6, 'fatol': 1e-12, 'maxiter': 20000})
    mu1_ogden, alpha1_ogden = res_ogden.x
    P_pred_ogden = ogden_stress(lambda_c, mu1_ogden, alpha1_ogden)
    r2_ogden = compute_r2(P_exp, P_pred_ogden)

    # 2. Données de relaxation
    rel_df = pd.read_csv(mat['rel_file'], sep=';', decimal=',')
    time_r_ms = rel_df['Temps total (s)'].values * 1000.0
    force_r = np.abs(rel_df['Effort_5N (N)'].values)

    t_start_ms = mat['rel_tstart'] * 1000.0
    t_peak_ms = mat['rel_tpeak'] * 1000.0
    t_end_ms = mat['rel_tend'] * 1000.0
    t_ramp_ms = t_peak_ms - t_start_ms

    mask_r = (time_r_ms >= t_start_ms) & (time_r_ms <= t_end_ms)
    t_all_ms = time_r_ms[mask_r] - t_start_ms
    F_all = force_r[mask_r]
    F_max_exp = np.max(F_all)

    # 3. Fit Prony n=1
    def obj_prony_n1(params):
        g1, tau1_ms = params
        if (g1 <= 0.0) or (g1 >= 1.0) or (tau1_ms <= 0):
            return 1e6
        F_pred = force_boltzmann_n1(t_all_ms, F_max_exp, t_ramp_ms, g1, tau1_ms)
        return np.mean((F_all - F_pred)**2)

    p0_prony_n1 = [0.15, 2000.0]
    bounds_prony_n1 = [(0.0001, 0.4), (1.0, 1000000.0)]
    res_prony_n1 = minimize(obj_prony_n1, p0_prony_n1, method='Nelder-Mead', bounds=bounds_prony_n1,
                            options={'xatol': 1e-6, 'fatol': 1e-12, 'maxiter': 20000, 'maxfev': 20000})
    g1_n1, tau1_n1_ms = res_prony_n1.x
    g_inf_n1 = 1.0 - g1_n1
    F_pred_n1 = force_boltzmann_n1(t_all_ms, F_max_exp, t_ramp_ms, g1_n1, tau1_n1_ms)
    r2_prony_n1 = compute_r2(F_all, F_pred_n1)

    # 4. Fit Prony n=2
    def obj_prony_n2(params):
        g1, tau1_ms, g2, tau2_ms = params
        if (g1 + g2 >= 1.0) or (tau1_ms <= 0) or (tau2_ms <= 0) or (tau1_ms >= tau2_ms):
            return 1e6
        F_pred = force_boltzmann_complete(t_all_ms, F_max_exp, t_ramp_ms, g1, tau1_ms, g2, tau2_ms)
        return np.mean((F_all - F_pred)**2)

    p0_prony_n2 = [0.15, 2000.0, 0.15, 50000.0]
    bounds_prony_n2 = [(0.0001, 0.4), (1.0, 20000.0), (0.0001, 0.4), (10000.0, 1000000.0)]
    res_prony_n2 = minimize(obj_prony_n2, p0_prony_n2, method='Nelder-Mead', bounds=bounds_prony_n2,
                            options={'xatol': 1e-6, 'fatol': 1e-12, 'maxiter': 20000, 'maxfev': 20000})
    g1_n2, tau1_n2_ms, g2_n2, tau2_n2_ms = res_prony_n2.x  
    g_inf_n2 = 1.0 - g1_n2 - g2_n2
    F_pred_n2 = force_boltzmann_complete(t_all_ms, F_max_exp, t_ramp_ms, g1_n2, tau1_n2_ms, g2_n2, tau2_n2_ms)
    r2_prony_n2 = compute_r2(F_all, F_pred_n2)

    # Stockage pour le .csv
    results.append({
        'Material': mat['name'],
        'Thickness_mm': mat['L0'],
        'mu1_NH_MPa': mu1_nh,
        'alpha1_NH': 2.0,
        'R2_NH': r2_nh,
        'mu1_Ogden_MPa': mu1_ogden,
        'alpha1_Ogden': alpha1_ogden,
        'R2_Ogden': r2_ogden,
        'Delta_R2_Hyper': r2_ogden - r2_nh,
        'Prony_n1_g1': g1_n1,
        'Prony_n1_tau1_ms': tau1_n1_ms,
        'Prony_n1_g_inf': g_inf_n1,
        'Prony_n1_R2': r2_prony_n1,
        'Prony_n2_g1': g1_n2,
        'Prony_n2_tau1_ms': tau1_n2_ms,
        'Prony_n2_g2': g2_n2,
        'Prony_n2_tau2_ms': tau2_n2_ms,
        'Prony_n2_g_inf': g_inf_n2,
        'Prony_n2_R2': r2_prony_n2,
        'Delta_R2_Visco': r2_prony_n2 - r2_prony_n1
    })

    # Graphiques
    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "font.size": 11,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "axes.grid": False,
        "lines.linewidth": 1.5
    })

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    step_c = max(1, len(eps_c) // 150)
    step = 80
    if mat['L0'] == 26.03:
        step = 10
        
    # Plot Hyperélastique
    axes[0].plot(eps_c[::step_c], P_exp[::step_c], 'k-', label='Experimental')
    axes[0].plot(eps_c[::step], P_pred_nh[::step], 'x', color='k', markersize=5, label='Neo-Hookean')
    axes[0].plot(eps_c[::step], P_pred_ogden[::step], 'o', markerfacecolor='none', color='k', markersize=6, label='Ogden')
    axes[0].set_xlabel(r'Déformation $\epsilon$ [%]')
    axes[0].set_ylabel(r'Contrainte $\sigma$ [MPa]')
    
    text_ogden = f"Néo-Hookéen:\n$\\mu_1 = {mu1_nh:.4f}$ MPa\n$R^2 = {r2_nh:.4f}$\n\nOgden:\n$\\mu_1 = {mu1_ogden:.4f}$ MPa\n$\\alpha_1 = {alpha1_ogden:.2f}$\n$R^2 = {r2_ogden:.4f}$"
    axes[0].text(0.05, 0.95, text_ogden, transform=axes[0].transAxes, fontsize=9,
                 verticalalignment='top', bbox=dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none'))
    axes[0].legend(loc='lower right', frameon=False, fontsize=9)

    # Plot Viscoélastique
    mask_charge = t_all_ms <= t_ramp_ms
    mask_relax = t_all_ms > t_ramp_ms

    t_charge = t_all_ms[mask_charge]
    F_pred_charge_n2 = (F_pred_n2 / F_max_exp)[mask_charge]
    F_pred_relax_n2 = (F_pred_n2 / F_max_exp)[mask_relax]
    
    F_pred_charge_n1 = (F_pred_n1 / F_max_exp)[mask_charge]
    F_pred_relax_n1 = (F_pred_n1 / F_max_exp)[mask_relax]

    t_relax = t_all_ms[mask_relax]
    step_charge = max(1, len(t_charge) // 7)
    step_relax = 500

    t_prony_pts = np.concatenate([t_charge[::step_charge], t_relax[::step_relax]])
    F_prony_pts_n2 = np.concatenate([F_pred_charge_n2[::step_charge], F_pred_relax_n2[::step_relax]])
    F_prony_pts_n1 = np.concatenate([F_pred_charge_n1[::step_charge], F_pred_relax_n1[::step_relax]])

    axes[1].plot(t_all_ms / 1000, F_all / F_max_exp, 'k-', linewidth=1, alpha=0.9, label='Experimental')          
    axes[1].plot(t_prony_pts / 1000, F_prony_pts_n1, 'x', color='k', markersize=5, label='Prony series (n=1)')
    axes[1].plot(t_prony_pts / 1000, F_prony_pts_n2, 'o', markerfacecolor='none', color='k', markersize=6, label='Prony series (n=2)')
    axes[1].set_xlabel(r'Time [s]')
    axes[1].set_ylabel(r'$F / F_{max}$')
    
    text_prony = f"$\\tau_1 = {tau1_n2_ms/1000:.3f}$ s, $g_1 = {g1_n2:.3f}$\n$\\tau_2 = {tau2_n2_ms/1000:.3f}$ s, $g_2 = {g2_n2:.3f}$\n$R^2 (n=2) = {r2_prony_n2:.4f}$\n$R^2 (n=1) = {r2_prony_n1:.4f}$"
    axes[1].text(0.95, 0.05, text_prony, transform=axes[1].transAxes, fontsize=10,
                 horizontalalignment='right', verticalalignment='bottom', 
                 bbox=dict(boxstyle='square', facecolor='white', alpha=0.9, edgecolor='none'))
    axes[1].legend(loc='center right', frameon=False)

    fig.suptitle(f'{mat["name"]}', fontsize=13, y=1.02)
    plt.tight_layout()
    
    filename = f"Fit_{mat['name'].replace(' ', '_').replace('(', '').replace(')', '')}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

# Export du .csv
results_df = pd.DataFrame(results)
results_df.to_csv('results_parameters.csv', index=False, sep=';', decimal=',')
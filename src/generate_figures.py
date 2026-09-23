import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os

# Academic aesthetic configuration
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'lines.linewidth': 2,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})

def load_simulation_data(run_dir):
    """Loads the generated data."""
    with open(os.path.join(run_dir, 'config.json'), 'r') as f:
        config = json.load(f)
    with open(os.path.join(run_dir, 'run_manifest.json'), 'r') as f:
        manifest = json.load(f)
    timeseries = pd.read_csv(os.path.join(run_dir, 'timeseries.csv'))
    snapshots = np.load(os.path.join(run_dir, 'snapshots.npz'))
    return config, manifest, timeseries, snapshots

def plot_figure3_spatial_profiles(run_dir, output_name="Figure3_Revised.png"):
    """
    Generates Figure 3 (4x3 spatial profiles) ensuring strict vertical order:
    1. Radius (R/R0)
    2. Pressure (p/Pin)
    3. Concentration (phi/phi_soil)
    4. Shear Stress (|tau_b|/|tau_b0|)
    """
    config, manifest, ts, snaps = load_simulation_data(run_dir)
    
    # Extract parameters
    L = config['L']
    R0 = config['R0']
    Pin = config['Pin']
    phi_soil = config['phi_soil']
    
    # Snapshots saved at t/ter = 0, 0.5, 1.0, 2.0, 2.5
    # Select 3 representative instances: t/ter = 0, 1.0, 2.5
    # ter is not stored directly — recalculate from t_final / t_final_over_ter
    t_final = manifest.get('t_final', None)
    t_final_over_ter = manifest.get('t_final_over_ter', None)
    if t_final is not None and t_final_over_ter:
        ter = t_final / t_final_over_ter
    else:
        ter = manifest.get('ter', 949.4)
    target_times = [0.0 * ter, 1.0 * ter, 2.5 * ter]
    
    fig, axes = plt.subplots(4, 3, figsize=(12, 10), sharex=True)
    x = np.linspace(0, L, config['Nx']) / L
    
    t_star = manifest.get('t_star_1D', None)

    for col_idx, target_t in enumerate(target_times):
        # Find the snapshot closest to the target time
        snap_indices = []
        i = 0
        while f"s{i:02d}_t" in snaps:
            snap_indices.append(i)
            i += 1
            
        best_i = snap_indices[0]
        min_diff = abs(snaps[f"s{best_i:02d}_t"][0] - target_t)
        for i in snap_indices[1:]:
            diff = abs(snaps[f"s{i:02d}_t"][0] - target_t)
            if diff < min_diff:
                min_diff = diff
                best_i = i
                
        t_actual = snaps[f"s{best_i:02d}_t"][0]
        R = snaps[f"s{best_i:02d}_R"]
        p = snaps[f"s{best_i:02d}_p"]
        phi = snaps[f"s{best_i:02d}_phi"]
        tau_b = snaps[f"s{best_i:02d}_tau_b"]
        
        tau_b0 = np.abs(tau_b[0]) if np.abs(tau_b[0]) > 0 else 1.0

        # Row 0: Radius
        axes[0, col_idx].plot(x, R / R0, color='tab:blue')
        if col_idx == 0: axes[0, col_idx].set_ylabel('$R / R_0$')
        axes[0, col_idx].set_title(f'$t/t_{{er}} \\approx {t_actual/ter:.1f}$')

        # Row 1: Pressure
        axes[1, col_idx].plot(x, p / Pin, color='tab:orange')
        if col_idx == 0: axes[1, col_idx].set_ylabel('$p / P_{in}$')

        # Row 2: Concentration
        axes[2, col_idx].plot(x, phi / phi_soil, color='tab:green')
        if col_idx == 0: axes[2, col_idx].set_ylabel(r'$\phi / \phi_{soil}$')

        # Row 3: Shear Stress
        axes[3, col_idx].plot(x, np.abs(tau_b) / tau_b0, color='tab:red')
        if col_idx == 0: axes[3, col_idx].set_ylabel('$|\\tau_b| / |\\tau_{b0}|$')
        axes[3, col_idx].set_xlabel('$x / L$')
        
        # Add visual diagnostic if 1D domain is breached
        if t_star is not None and t_actual > t_star:
            for row in range(4):
                axes[row, col_idx].set_facecolor('#ffe6e6') # Pale red background to indicate loss of validity
                if row == 0:
                    axes[row, col_idx].text(0.5, 0.5, '1D Limit\nExceeded', 
                                            transform=axes[row, col_idx].transAxes,
                                            ha='center', va='center', color='red', alpha=0.5, fontsize=14)

    plt.tight_layout()
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"[OK] {output_name} successfully generated.")

def plot_validity_domain(run_dir_Kout0, output_name="Figure_Validity_1D.png"):
    """
    Generates a new figure showing the temporal evolution of validity indicators
    (R_out/L and max|dR/dx|) to rigorously demonstrate the identification of t*.
    """
    config, manifest, ts, _ = load_simulation_data(run_dir_Kout0)
    t_final = manifest.get('t_final', None)
    t_final_over_ter = manifest.get('t_final_over_ter', None)
    if t_final is not None and t_final_over_ter:
        ter = t_final / t_final_over_ter
    else:
        ter = manifest.get('ter', 949.4)
    t_star = manifest.get('t_star_1D', None)
    
    fig, ax1 = plt.subplots(figsize=(8, 5))
    
    t_norm = ts['t'] / ter
    
    # Axis 1: Slenderness ratio
    color1 = 'tab:blue'
    ax1.set_xlabel('Dimensionless time $t / t_{er}$')
    ax1.set_ylabel('Slenderness ratio $R_{out} / L$', color=color1)
    ax1.plot(t_norm, ts['slenderness_ratio'], color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.axhline(y=0.1, color=color1, linestyle=':', label='Validity threshold (0.1)')
    
    # Axis 2: Local slope
    ax2 = ax1.twinx()
    color2 = 'tab:purple'
    ax2.set_ylabel(r'Max local taper $\max|\partial R / \partial x|$', color=color2)
    ax2.plot(t_norm, ts['max_dR_dx'], color=color2, linestyle='--')
    ax2.tick_params(axis='y', labelcolor=color2)
    
    # Vertical line t*
    if t_star is not None:
        t_star_norm = t_star / ter
        ax1.axvline(x=t_star_norm, color='red', linewidth=2, label=f'$t^* \\approx {t_star_norm:.2f}\\,t_{{er}}$')
        ax1.fill_between(t_norm, 0, 1, where=(t_norm > t_star_norm), color='red', alpha=0.1,
                         transform=ax1.get_xaxis_transform())
        # Annotation placed in the valid zone (left of t*), at the top of the graph
        x_annot = max(t_star_norm * 0.5, t_norm.min() + 0.01)
        ax1.text(x_annot, 0.92, '1D assumption\nbreaks down →', color='red', fontsize=9,
                 ha='center', va='top', transform=ax1.get_xaxis_transform())
        
    fig.tight_layout()
    fig.legend(loc='upper left', bbox_to_anchor=(0.15, 0.9))
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"[OK] {output_name} successfully generated.")

if __name__ == "__main__":
    # Update these paths with the exact names of the generated folders
    # Example: 'output/20260802_231742_A_Kout10'
    
    # Figure 3 (Case K_out = 10, which remains 1D valid)
    run_A_Kout10 = "output_etape3/20260803_203802_A_Kout10" 
    if os.path.exists(run_A_Kout10):
        plot_figure3_spatial_profiles(run_A_Kout10, "Figure3_Kout10_Revised.png")
    else:
        print(f"[WARN] Folder {run_A_Kout10} not found.")

    # New Validity figure (Case K_out = 0, which crosses t*)
    run_A_Kout0 = "output_etape3/20260803_183817_A_Kout0"
    if os.path.exists(run_A_Kout0):
        plot_validity_domain(run_A_Kout0, "Figure_Validity_Kout0.png")
    else:
        print(f"[WARN] Folder {run_A_Kout0} not found.")

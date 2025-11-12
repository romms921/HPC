#!/usr/bin/env python
import os
import traceback
import numpy as np
import pandas as pd
import re
import sys
import shutil
import uuid
from itertools import product, islice

# IMPORTANT: Ensure the 'glafic' library is in your Python path
import glafic

# --- File Paths & Parameters ---
model_output_base = '/home/rommulus/Projects/itng_lensing/Simulations/Output'
base_results_path = '/home/rommulus/Projects/itng_lensing/Simulations/Input/System_2'
obs_point_file = os.path.join(base_results_path, 'pos+flux+td_point.dat')
SCRATCH_DIR = os.getenv('SCRATCH_DIR', '/tmp')

# --- Simulation Parameters ---
sim_name = 'Sim 8'
model = 'POW' # Model name for output files
m = np.linspace(1.5, 2.2, 100)
n = np.linspace(0.0, 0.9, 100)
o = np.linspace(0, 360, 100)

m_lens, m_param = 1, 8
n_lens, n_param = 1, 5
o_lens, o_param = 1, 6

constraint_file = os.path.join(base_results_path, 'pos_point.dat')
prior_file = os.path.join(base_results_path, 'prior.dat')

# --- Set these flags to True to enable the corresponding calculations ---
time_delay = False
h0 = False
critical_curve = False

# --- Directory setup ---
final_results_dir = os.path.join(model_output_base, 'individual_results')
progress_status_dir = os.path.join(model_output_base, 'progress_status')
os.makedirs(final_results_dir, exist_ok=True)
os.makedirs(progress_status_dir, exist_ok=True)

# --- Performance Tuning ---
PROGRESS_UPDATE_INTERVAL = 100
CSV_BATCH_SIZE = 100

# --- Parameter definitions (used for CSV headers) ---
pow_params = ['$z_{s,fid}$', 'x', 'y', 'e', '$θ_{e}$', '$r_{Ein}$', '$\gamma$ (PWI)']
sie_params = ['$\sigma$', 'x', 'y', 'e', '$θ_{e}$', '$r_{core}$', 'NaN']
nfw_params = ['M', 'x', 'y', 'e', '$θ_{e}$', 'c or $r_{s}$', 'NaN']
pert_params = ['$z_{s,fid}$', 'x', 'y', '$\gamma$', '$θ_{\gamma}$', 'NaN', '$\kappa$']
model_params = {'POW': pow_params, 'SIE': sie_params, 'ANFW': nfw_params, 'PERT': pert_params}


# --- Helper Functions ---
def calculate_distance(row1, row2):
    return np.sqrt((row1['x'] - row2['x'])**2 + (row1['y'] - row2['y'])**2)

def assign_image_names(df):
    if df.empty or 'Img' in df.columns and pd.notna(df['Img']).all():
        return df
    names = ['A', 'B', 'C', 'D']
    if len(df) > len(names): names.extend([f"Img{i+1}" for i in range(len(names), len(df))])
    
    start_index = df['mag'].idxmax()
    centroid_x, centroid_y = df['x'].mean(), df['y'].mean()
    
    angles = {i: np.arctan2(row['y'] - centroid_y, row['x'] - centroid_x) for i, row in df.iterrows()}
    sorted_indices = sorted(angles, key=angles.get, reverse=True)
    
    start_pos = sorted_indices.index(start_index)
    
    final_names = {}
    for i in range(len(df)):
        original_index = sorted_indices[(start_pos + i) % len(df)]
        final_names[original_index] = names[i]

    df['Img'] = df.index.map(final_names)
    return df

def rms_extract(model_ver, model_path, obs_point_df):
    results = {}
    
    opt_result_file = os.path.join(model_path, f'{model_ver}_optresult.dat')
    with open(opt_result_file, 'r') as file: opt_result = file.readlines()
    last_optimize_index = next((i for i, line in reversed(list(enumerate(opt_result))) if 'optimize' in line), None)
    if last_optimize_index is None: raise ValueError("No 'optimize' line found.")
    opt_result = opt_result[last_optimize_index + 1:]
    
    lens_params_dict = {re.split(r'\s+', line.strip())[1]: [float(x) for x in re.split(r'\s+', line.strip())[3:]] for line in opt_result if line.startswith('lens')}
    source_params = [[float(x) for x in re.split(r'\s+', line.strip())[1:]] for line in opt_result if line.startswith('point')]
    chi2_line = next((line for line in opt_result if 'chi^2' in line), None)
    
    results['chi2'] = float(chi2_line.split('=')[-1].strip().split()[0]) if chi2_line else np.nan
    results['source_x'] = source_params[0][1] if source_params else np.nan
    results['source_y'] = source_params[0][2] if source_params else np.nan
    results['lens_params'] = lens_params_dict

    if h0:
        hubble_line = next((line for line in opt_result if 'hubble =' in line), None)
        results['h0'] = float(hubble_line.split('=')[-1].strip().split()[0]) if hubble_line else np.nan

    obs_point = assign_image_names(obs_point_df.copy())
    out_point_file = os.path.join(model_path, f'{model_ver}_point.dat')
    out_point = pd.read_csv(out_point_file, sep='\s+', header=None, skiprows=1, names=['x', 'y', 'mag', 'td', 'c5', 'c6', 'c7', 'c8'])

    # Match predicted to observed images robustly
    merged_df = pd.merge(obs_point.add_prefix('obs_'), out_point.add_prefix('pred_'), how='cross')
    merged_df['dist'] = np.sqrt((merged_df['obs_x'] - merged_df['pred_x'])**2 + (merged_df['obs_y'] - merged_df['pred_y'])**2)
    merged_df = merged_df.loc[merged_df.groupby('obs_Img')['dist'].idxmin()]

    # --- Aggregate and Per-Image Metrics ---
    results['pos_rms'] = np.sqrt(np.mean(merged_df['dist']**2))
    results['mag_rms'] = np.sqrt(np.mean((np.abs(merged_df['obs_mag']) - np.abs(merged_df['pred_mag']))**2))
    
    # --- PERCENTAGE ERROR CALCULATIONS (RESTORED) ---
    # Calculate percentage error for magnification
    mag_abs_obs = np.abs(merged_df['obs_mag'])
    mag_abs_pred = np.abs(merged_df['pred_mag'])
    # Avoid division by zero if an observed magnitude is zero
    mag_perc_err = np.divide(np.abs(mag_abs_pred - mag_abs_obs), mag_abs_obs, out=np.full_like(mag_abs_obs, np.nan), where=mag_abs_obs!=0) * 100
    results['avg_mag_per'] = np.nanmean(mag_perc_err)
    
    # Store per-image results
    for i, row in merged_df.iterrows():
        img_name = row['obs_Img']
        results[f'pos_{img_name}'] = row['dist']
        results[f'mag_{img_name}'] = row['pred_mag']
        results[f'mag_err_{img_name}'] = row['pred_mag'] - row['obs_mag']
        # --- Store per-image percentage error (RESTORED) ---
        # Find the corresponding percentage error from our calculation above
        original_index = merged_df.index.get_loc(i)
        results[f'mag_per_{img_name}'] = mag_perc_err.iloc[original_index]

    if time_delay:
        td_err = merged_df['pred_td'] - merged_df['obs_td']
        results['td_rms'] = np.sqrt(np.mean(td_err**2))

        # --- TD PERCENTAGE ERROR CALCULATIONS (RESTORED) ---
        obs_td = merged_df['obs_td']
        # Avoid division by zero if an observed time delay is zero
        td_perc_err = np.divide(np.abs(td_err), obs_td, out=np.full_like(obs_td, np.nan), where=obs_td!=0) * 100
        results['avg_td_per'] = np.nanmean(td_perc_err)

        for i, row in merged_df.iterrows():
            img_name = row['obs_Img']
            results[f'td_{img_name}'] = row['pred_td']
            results[f'td_err_{img_name}'] = row['pred_td'] - row['obs_td']
            # --- Store per-image TD percentage error (RESTORED) ---
            original_index = merged_df.index.get_loc(i)
            results[f'td_per_{img_name}'] = td_perc_err.iloc[original_index]
            
    return results


def run_glafic_calculation(params, model_name, worker_temp_dir):
    m_val, n_val, o_val = params
    output_path = os.path.join(worker_temp_dir, model_name)
    
    base_lens_params = [0.261343256161012, 1.0, 0.0, 0.0, 0.107, 23.38, 0.41, 2.0]
    
    current_lens_params = list(base_lens_params)
    current_lens_params[n_param - 1] = n_val
    current_lens_params[o_param - 1] = o_val
    current_lens_params[m_param - 1] = m_val

    glafic.init(0.3, 0.7, -1.0, 0.7, output_path, -3.0, -3.0, 3.0, 3.0, 0.01, 0.01, 1, verb=0)
    glafic.set_secondary('chi2_splane 1', verb=0)
    glafic.set_secondary('chi2_checknimg 0', verb=0)
    glafic.set_secondary('chi2_restart   -1', verb=0)
    glafic.set_secondary('chi2_usemag    1', verb=0)
    glafic.set_secondary('hvary          0' if not h0 else '1', verb=0)
    glafic.set_secondary('ran_seed -122000', verb=0)
    glafic.startup_setnum(1, 0, 1)
    glafic.set_lens(1, 'pow', *current_lens_params)
    glafic.set_point(1, 1.0, 0.0, 0.0)
    glafic.setopt_lens(1, 0, 0, 1, 1, 1, 1, 1, 1)
    glafic.setopt_point(1, 0, 1, 1)
    glafic.model_init(verb=0)
    glafic.readobs_point(constraint_file)
    if prior_file and os.path.exists(prior_file):
        glafic.parprior(prior_file)
    glafic.optimize()
    glafic.findimg()
    
    if critical_curve:
        glafic.writecrit(1.0)
        
    glafic.quit()


def run_single_model(params, worker_temp_dir, obs_point_df):
    m_val, n_val, o_val = round(params[0], 5), round(params[1], 5), round(params[2], 5)
    model_name = f'{model}_{m_val}_{n_val}_{o_val}'

    try:
        run_glafic_calculation((m_val, n_val, o_val), model_name, worker_temp_dir)
        
        extracted_results = rms_extract(model_name, worker_temp_dir, obs_point_df)
        
        out_point_file = os.path.join(worker_temp_dir, f'{model_name}_point.dat')
        num_images = sum(1 for line in open(out_point_file) if line.strip()) - 1 if os.path.exists(out_point_file) else 0
        
        final_result_dict = {'m': m_val, 'n': n_val, 'o': o_val, 'num_images': num_images}
        final_result_dict.update(extracted_results)
        
        lens_params = final_result_dict.pop('lens_params')
        for lens_name, lens_p in lens_params.items():
            model_type = lens_name.upper()
            if model_type not in model_params: continue
            param_names = model_params[model_type]
            for i, param_val in enumerate(lens_p):
                if i < len(param_names) and param_names[i] != 'NaN':
                    col_name = f"{model_type}_{param_names[i]}"
                    final_result_dict[col_name] = param_val

        return final_result_dict
    
    except Exception:
        print(f"!!! WORKER FAILED on model {model_name}.")
        traceback.print_exc()
        sys.stdout.flush()
        return None
    finally:
        for suffix in ['_optresult.dat', '_point.dat', '_crit.dat']:
             f = os.path.join(worker_temp_dir, f"{model_name}{suffix}")
             if os.path.exists(f):
                 try: os.remove(f)
                 except OSError: pass

def write_batch_to_csv(batch, csv_file):
    if not batch: return
    df = pd.DataFrame(batch)
    header = not os.path.exists(csv_file)
    df.to_csv(csv_file, mode='a', header=header, index=False)

def main():
    rank = int(os.getenv('SLURM_PROCID', 0))
    size = int(os.getenv('SLURM_NTASKS', 1))
    
    all_tasks_generator = product(m, n, o)
    tasks_for_this_worker = islice(all_tasks_generator, rank, None, size)
    
    total_combinations = len(m) * len(n) * len(o)
    total_for_this_worker = (total_combinations - rank + size - 1) // size
    
    output_csv_file = os.path.join(final_results_dir, f"results_rank_{rank}.csv")
    status_file = os.path.join(progress_status_dir, f"worker_{rank}.status")
    
    obs_point_df = pd.read_csv(obs_point_file, sep='\s+', header=None, skiprows=1, names=['x', 'y', 'mag', 'pos_err', 'mag_err', 'td', 'td_err', 'parity'])
    
    worker_temp_dir = os.path.join(SCRATCH_DIR, f'worker_temp_{rank}_{uuid.uuid4()}')
    os.makedirs(worker_temp_dir, exist_ok=True)
    
    results_batch = []
    num_completed = 0
    
    try:
        for i, task_params in enumerate(tasks_for_this_worker):
            result = run_single_model(task_params, worker_temp_dir, obs_point_df)
            
            if result:
                results_batch.append(result)
                num_completed += 1
            
            if len(results_batch) >= CSV_BATCH_SIZE or ((i + 1) == total_for_this_worker and results_batch):
                write_batch_to_csv(results_batch, output_csv_file)
                results_batch = [] 

            if (i + 1) % PROGRESS_UPDATE_INTERVAL == 0 or (i + 1) == total_for_this_worker:
                with open(status_file, 'w') as f:
                    f.write(str(i + 1))
    finally:
        if os.path.exists(worker_temp_dir):
            shutil.rmtree(worker_temp_dir)

    print(f"Worker {rank} finished. Total models attempted: {total_for_this_worker}. Successful: {num_completed}.")
    sys.stdout.flush()

if __name__ == "__main__":
    main()